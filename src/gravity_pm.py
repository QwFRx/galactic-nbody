"""Particle-Mesh gravity with CIC assignment and FFT Poisson solves.

Two boundary conditions are supported:

* periodic (cosmic web): standard FFT Poisson solver, k = 0 removed;
* isolated (galaxy collision): Hockney-Eastwood zero padding to 2N and the
  exact isolated Green's function G(r) = -G/r, regularised at r = 0.

The CIC window is deconvolved once (divide by the sinc^2-based window) and a
Gaussian exp(-k^2 r_s^2) filter suppresses grid noise.  Forces are obtained
from a 4-point finite-difference gradient and interpolated back with the
same CIC kernel, which conserves momentum.
"""

from __future__ import annotations

import os
import warnings

import numpy as np

import config

try:
    from scipy import fft as _scipy_fft

    def _rfftn(array: np.ndarray) -> np.ndarray:
        return _scipy_fft.rfftn(array, workers=-1)

    def _irfftn(array: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
        return _scipy_fft.irfftn(array, s=shape, workers=-1)

    HAVE_SCIPY_FFT = True
except Exception:  # pragma: no cover - scipy is a declared dependency
    def _rfftn(array: np.ndarray) -> np.ndarray:
        return np.fft.rfftn(array)

    def _irfftn(array: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
        return np.fft.irfftn(array, s=shape)

    HAVE_SCIPY_FFT = False

try:
    from numba import (njit, prange, get_num_threads, get_thread_id,
                       set_num_threads)
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda fn: fn

    prange = range

    def get_num_threads() -> int:
        return 1

    def get_thread_id() -> int:
        return 0

    def set_num_threads(n: int) -> None:
        return None


# ---------------------------------------------------------------------------
# Numba kernels
# ---------------------------------------------------------------------------
@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _deposit_periodic_numba(pos: np.ndarray, mass: np.ndarray,
                            ox: float, oy: float, oz: float, inv_h: float,
                            n: int, bufs: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        tid = get_thread_id()
        gx = (pos[p, 0] - ox) * inv_h - 0.5
        gy = (pos[p, 1] - oy) * inv_h - 0.5
        gz = (pos[p, 2] - oz) * inv_h - 0.5
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gy))
        k0 = int(np.floor(gz))
        fx = gx - i0
        fy = gy - j0
        fz = gz - k0
        i0 %= n
        j0 %= n
        k0 %= n
        i1 = (i0 + 1) % n
        j1 = (j0 + 1) % n
        k1 = (k0 + 1) % n
        w = mass[p]
        w00 = w * (1.0 - fx) * (1.0 - fy)
        w10 = w * fx * (1.0 - fy)
        w01 = w * (1.0 - fx) * fy
        w11 = w * fx * fy
        bufs[tid, i0, j0, k0] += w00 * (1.0 - fz)
        bufs[tid, i1, j0, k0] += w10 * (1.0 - fz)
        bufs[tid, i0, j1, k0] += w01 * (1.0 - fz)
        bufs[tid, i1, j1, k0] += w11 * (1.0 - fz)
        bufs[tid, i0, j0, k1] += w00 * fz
        bufs[tid, i1, j0, k1] += w10 * fz
        bufs[tid, i0, j1, k1] += w01 * fz
        bufs[tid, i1, j1, k1] += w11 * fz


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _deposit_isolated_numba(pos: np.ndarray, mass: np.ndarray,
                            ox: float, oy: float, oz: float, inv_h: float,
                            n: int, offset: int, bufs: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        gx = (pos[p, 0] - ox) * inv_h - 0.5
        gy = (pos[p, 1] - oy) * inv_h - 0.5
        gz = (pos[p, 2] - oz) * inv_h - 0.5
        if gx < -0.5 or gx >= n - 0.5:
            continue
        if gy < -0.5 or gy >= n - 0.5:
            continue
        if gz < -0.5 or gz >= n - 0.5:
            continue
        tid = get_thread_id()
        i0 = int(np.floor(gx)) + offset
        j0 = int(np.floor(gy)) + offset
        k0 = int(np.floor(gz)) + offset
        fx = gx - np.floor(gx)
        fy = gy - np.floor(gy)
        fz = gz - np.floor(gz)
        i1 = i0 + 1
        j1 = j0 + 1
        k1 = k0 + 1
        w = mass[p]
        w00 = w * (1.0 - fx) * (1.0 - fy)
        w10 = w * fx * (1.0 - fy)
        w01 = w * (1.0 - fx) * fy
        w11 = w * fx * fy
        bufs[tid, i0, j0, k0] += w00 * (1.0 - fz)
        bufs[tid, i1, j0, k0] += w10 * (1.0 - fz)
        bufs[tid, i0, j1, k0] += w01 * (1.0 - fz)
        bufs[tid, i1, j1, k0] += w11 * (1.0 - fz)
        bufs[tid, i0, j0, k1] += w00 * fz
        bufs[tid, i1, j0, k1] += w10 * fz
        bufs[tid, i0, j1, k1] += w01 * fz
        bufs[tid, i1, j1, k1] += w11 * fz


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _gradient_periodic_numba(phi: np.ndarray, inv_12h: float,
                             out: np.ndarray) -> None:
    n = phi.shape[0]
    for i in prange(n):
        im2 = (i - 2) % n
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        ip2 = (i + 2) % n
        for j in range(n):
            jm2 = (j - 2) % n
            jm1 = (j - 1) % n
            jp1 = (j + 1) % n
            jp2 = (j + 2) % n
            for k in range(n):
                km2 = (k - 2) % n
                km1 = (k - 1) % n
                kp1 = (k + 1) % n
                kp2 = (k + 2) % n
                out[i, j, k, 0] = -(phi[im2, j, k] - 8.0 * phi[im1, j, k]
                                    + 8.0 * phi[ip1, j, k] - phi[ip2, j, k]) * inv_12h
                out[i, j, k, 1] = -(phi[i, jm2, k] - 8.0 * phi[i, jm1, k]
                                    + 8.0 * phi[i, jp1, k] - phi[i, jp2, k]) * inv_12h
                out[i, j, k, 2] = -(phi[i, j, km2] - 8.0 * phi[i, j, km1]
                                    + 8.0 * phi[i, j, kp1] - phi[i, j, kp2]) * inv_12h


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _gradient_isolated_numba(phi: np.ndarray, inv_h: float,
                             inv_12h: float, out: np.ndarray) -> None:
    m = phi.shape[0]
    for i in prange(m):
        for j in range(m):
            for k in range(m):
                if i >= 2 and i < m - 2:
                    out[i, j, k, 0] = -(phi[i - 2, j, k] - 8.0 * phi[i - 1, j, k]
                                        + 8.0 * phi[i + 1, j, k]
                                        - phi[i + 2, j, k]) * inv_12h
                elif i == 0:
                    out[i, j, k, 0] = (phi[0, j, k] - phi[1, j, k]) * inv_h
                else:
                    out[i, j, k, 0] = (phi[i - 1, j, k] - phi[i + 1, j, k]) * 0.5 * inv_h
                if j >= 2 and j < m - 2:
                    out[i, j, k, 1] = -(phi[i, j - 2, k] - 8.0 * phi[i, j - 1, k]
                                        + 8.0 * phi[i, j + 1, k]
                                        - phi[i, j + 2, k]) * inv_12h
                elif j == 0:
                    out[i, j, k, 1] = (phi[i, 0, k] - phi[i, 1, k]) * inv_h
                else:
                    out[i, j, k, 1] = (phi[i, j - 1, k] - phi[i, j + 1, k]) * 0.5 * inv_h
                if k >= 2 and k < m - 2:
                    out[i, j, k, 2] = -(phi[i, j, k - 2] - 8.0 * phi[i, j, k - 1]
                                        + 8.0 * phi[i, j, k + 1]
                                        - phi[i, j, k + 2]) * inv_12h
                elif k == 0:
                    out[i, j, k, 2] = (phi[i, j, 0] - phi[i, j, 1]) * inv_h
                else:
                    out[i, j, k, 2] = (phi[i, j, k - 1] - phi[i, j, k + 1]) * 0.5 * inv_h


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _interp_periodic_numba(acc: np.ndarray, pos: np.ndarray,
                           ox: float, oy: float, oz: float, inv_h: float,
                           n: int, out: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        gx = (pos[p, 0] - ox) * inv_h - 0.5
        gy = (pos[p, 1] - oy) * inv_h - 0.5
        gz = (pos[p, 2] - oz) * inv_h - 0.5
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gy))
        k0 = int(np.floor(gz))
        fx = gx - i0
        fy = gy - j0
        fz = gz - k0
        i0 %= n
        j0 %= n
        k0 %= n
        i1 = (i0 + 1) % n
        j1 = (j0 + 1) % n
        k1 = (k0 + 1) % n
        for c in range(3):
            out[p, c] = (
                acc[i0, j0, k0, c] * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
                + acc[i1, j0, k0, c] * fx * (1.0 - fy) * (1.0 - fz)
                + acc[i0, j1, k0, c] * (1.0 - fx) * fy * (1.0 - fz)
                + acc[i1, j1, k0, c] * fx * fy * (1.0 - fz)
                + acc[i0, j0, k1, c] * (1.0 - fx) * (1.0 - fy) * fz
                + acc[i1, j0, k1, c] * fx * (1.0 - fy) * fz
                + acc[i0, j1, k1, c] * (1.0 - fx) * fy * fz
                + acc[i1, j1, k1, c] * fx * fy * fz)


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _interp_isolated_numba(acc: np.ndarray, pos: np.ndarray,
                           ox: float, oy: float, oz: float, inv_h: float,
                           m: int, vx0: float, vy0: float, vz0: float,
                           vx1: float, vy1: float, vz1: float,
                           out: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        px = pos[p, 0]
        py = pos[p, 1]
        pz = pos[p, 2]
        if px < vx0 or px >= vx1 or py < vy0 or py >= vy1 or pz < vz0 or pz >= vz1:
            out[p, 0] = 0.0
            out[p, 1] = 0.0
            out[p, 2] = 0.0
            continue
        gx = (px - ox) * inv_h - 0.5
        gy = (py - oy) * inv_h - 0.5
        gz = (pz - oz) * inv_h - 0.5
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gy))
        k0 = int(np.floor(gz))
        fx = gx - i0
        fy = gy - j0
        fz = gz - k0
        i1 = i0 + 1
        j1 = j0 + 1
        k1 = k0 + 1
        if i0 < 0 or j0 < 0 or k0 < 0 or i1 >= m or j1 >= m or k1 >= m:
            out[p, 0] = 0.0
            out[p, 1] = 0.0
            out[p, 2] = 0.0
            continue
        for c in range(3):
            out[p, c] = (
                acc[i0, j0, k0, c] * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
                + acc[i1, j0, k0, c] * fx * (1.0 - fy) * (1.0 - fz)
                + acc[i0, j1, k0, c] * (1.0 - fx) * fy * (1.0 - fz)
                + acc[i1, j1, k0, c] * fx * fy * (1.0 - fz)
                + acc[i0, j0, k1, c] * (1.0 - fx) * (1.0 - fy) * fz
                + acc[i1, j0, k1, c] * fx * (1.0 - fy) * fz
                + acc[i0, j1, k1, c] * (1.0 - fx) * fy * fz
                + acc[i1, j1, k1, c] * fx * fy * fz)


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _interp_scalar_periodic_numba(field: np.ndarray, pos: np.ndarray,
                                  ox: float, oy: float, oz: float,
                                  inv_h: float, n: int,
                                  out: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        gx = (pos[p, 0] - ox) * inv_h - 0.5
        gy = (pos[p, 1] - oy) * inv_h - 0.5
        gz = (pos[p, 2] - oz) * inv_h - 0.5
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gy))
        k0 = int(np.floor(gz))
        fx = gx - i0
        fy = gy - j0
        fz = gz - k0
        i0 %= n
        j0 %= n
        k0 %= n
        i1 = (i0 + 1) % n
        j1 = (j0 + 1) % n
        k1 = (k0 + 1) % n
        out[p] = (
            field[i0, j0, k0] * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
            + field[i1, j0, k0] * fx * (1.0 - fy) * (1.0 - fz)
            + field[i0, j1, k0] * (1.0 - fx) * fy * (1.0 - fz)
            + field[i1, j1, k0] * fx * fy * (1.0 - fz)
            + field[i0, j0, k1] * (1.0 - fx) * (1.0 - fy) * fz
            + field[i1, j0, k1] * fx * (1.0 - fy) * fz
            + field[i0, j1, k1] * (1.0 - fx) * fy * fz
            + field[i1, j1, k1] * fx * fy * fz)


@njit(parallel=True, fastmath=True, nogil=True, cache=True)
def _interp_scalar_isolated_numba(field: np.ndarray, pos: np.ndarray,
                                  ox: float, oy: float, oz: float,
                                  inv_h: float, m: int,
                                  vx0: float, vy0: float, vz0: float,
                                  vx1: float, vy1: float, vz1: float,
                                  out: np.ndarray) -> None:
    count = pos.shape[0]
    for p in prange(count):
        px = pos[p, 0]
        py = pos[p, 1]
        pz = pos[p, 2]
        if px < vx0 or px >= vx1 or py < vy0 or py >= vy1 or pz < vz0 or pz >= vz1:
            out[p] = 0.0
            continue
        gx = (px - ox) * inv_h - 0.5
        gy = (py - oy) * inv_h - 0.5
        gz = (pz - oz) * inv_h - 0.5
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gy))
        k0 = int(np.floor(gz))
        fx = gx - i0
        fy = gy - j0
        fz = gz - k0
        i1 = i0 + 1
        j1 = j0 + 1
        k1 = k0 + 1
        if i0 < 0 or j0 < 0 or k0 < 0 or i1 >= m or j1 >= m or k1 >= m:
            out[p] = 0.0
            continue
        out[p] = (
            field[i0, j0, k0] * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
            + field[i1, j0, k0] * fx * (1.0 - fy) * (1.0 - fz)
            + field[i0, j1, k0] * (1.0 - fx) * fy * (1.0 - fz)
            + field[i1, j1, k0] * fx * fy * (1.0 - fz)
            + field[i0, j0, k1] * (1.0 - fx) * (1.0 - fy) * fz
            + field[i1, j0, k1] * fx * (1.0 - fy) * fz
            + field[i0, j1, k1] * (1.0 - fx) * fy * fz
            + field[i1, j1, k1] * fx * fy * fz)


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
class PMGravity:
    """Particle-Mesh Poisson solver shared by both scenes."""

    def __init__(self, n: int, box: float, periodic: bool, G: float,
                 center: np.ndarray | None = None,
                 smoothing_cells: float = config.SMOOTHING_CELLS,
                 nthreads: int | None = None) -> None:
        self.n = int(n)
        self.box = float(box)
        self.periodic = bool(periodic)
        self.G = float(G)
        self.h = self.box / self.n
        self.offset = 0 if self.periodic else self.n // 2
        self.nfft = self.n if self.periodic else 2 * self.n
        if center is None:
            center = np.full(3, self.box / 2.0, dtype=np.float64)
        self.center = np.asarray(center, dtype=np.float64)
        self.origin = (np.zeros(3, dtype=np.float64) if self.periodic
                       else self.center - 0.5 * self.box)
        self.smoothing = float(smoothing_cells) * self.h
        self._setup_frequencies()
        self._setup_threads(nthreads)
        self._setup_kernel()
        self.last_density: np.ndarray | None = None
        self.last_phi_ext: np.ndarray | None = None
        self.last_phi: np.ndarray | None = None
        self.last_origin_ext = self.origin - 2.0 * self.h
        self.outside = 0

    # -- setup -------------------------------------------------------------
    def _setup_frequencies(self) -> None:
        f = self.nfft
        fx = np.fft.fftfreq(f, d=1.0 / f)
        fz = np.fft.rfftfreq(f, d=1.0 / f)
        # CIC assignment window: product of sinc^2(k_i h / 2).
        window = (np.sinc(fx[:, None, None] / f)
                  * np.sinc(fx[None, :, None] / f)
                  * np.sinc(fz[None, None, :] / f))
        self.window = window.astype(np.float32)
        kx = 2.0 * np.pi * fx / self.box
        kz = 2.0 * np.pi * fz / self.box
        k2 = (kx[:, None, None] ** 2 + kx[None, :, None] ** 2
              + kz[None, None, :] ** 2)
        self.k2 = k2
        self.gauss = np.exp(-k2 * self.smoothing ** 2).astype(np.float32)

    def _setup_threads(self, nthreads: int | None) -> None:
        cpu = os.cpu_count() or 1
        volume = self.nfft ** 3
        per_thread = 4 * volume
        budget = 1_100_000_000
        auto = max(1, min(cpu, int(budget / max(per_thread, 1))))
        self.nthreads = int(nthreads) if nthreads else auto
        self.nthreads = max(1, min(self.nthreads, cpu))

    def _setup_kernel(self) -> None:
        f = self.nfft
        if self.periodic:
            base = np.zeros_like(self.k2, dtype=np.float32)
            mask = self.k2 > 0.0
            base[mask] = (-4.0 * np.pi * self.gauss[mask]
                          / (self.k2[mask] * self.window[mask]
                             * self.h ** 3)).astype(np.float32)
            self.kernel = base[..., :f // 2 + 1]
        else:
            idx = np.arange(f)
            dx = np.where(idx < f // 2, idx, idx - f) * self.h
            r2 = (dx[:, None, None] ** 2 + dx[None, :, None] ** 2
                  + dx[None, None, :] ** 2)
            soft = 0.5 * self.h
            green = -self.G / np.sqrt(r2 + soft ** 2)
            green_k = _rfftn(green.astype(np.float32))
            self.kernel = ((green_k * self.gauss) / self.window).astype(np.complex64)

    def set_G(self, G: float) -> None:
        """Update the gravitational constant (used by the cosmic solver)."""
        self.G = float(G)
        if not self.periodic:
            self._setup_kernel()

    # -- mass assignment ---------------------------------------------------
    def _deposit_numpy(self, pos: np.ndarray, mass: np.ndarray) -> np.ndarray:
        n = self.n
        h = self.h
        f = self.nfft
        ox, oy, oz = self.origin
        gx = (pos[:, 0] - ox) / h - 0.5
        gy = (pos[:, 1] - oy) / h - 0.5
        gz = (pos[:, 2] - oz) / h - 0.5
        if self.periodic:
            valid = np.ones(pos.shape[0], dtype=bool)
        else:
            valid = ((gx >= -0.5) & (gx < n - 0.5) & (gy >= -0.5)
                     & (gy < n - 0.5) & (gz >= -0.5) & (gz < n - 0.5))
        ix = np.floor(gx).astype(np.int64)
        iy = np.floor(gy).astype(np.int64)
        iz = np.floor(gz).astype(np.int64)
        fx = gx - ix
        fy = gy - iy
        fz = gz - iz
        if self.periodic:
            ix %= n
            iy %= n
            iz %= n
            ix1 = (ix + 1) % n
            iy1 = (iy + 1) % n
            iz1 = (iz + 1) % n
        else:
            ix = ix + self.offset
            iy = iy + self.offset
            iz = iz + self.offset
            ix1 = ix + 1
            iy1 = iy + 1
            iz1 = iz + 1
        flat_size = f ** 3
        out = np.zeros(flat_size, dtype=np.float64)
        w = mass.astype(np.float64)[valid]
        for dx, dy, dz, ww in (
            (ix, iy, iz, (1 - fx) * (1 - fy) * (1 - fz)),
            (ix1, iy, iz, fx * (1 - fy) * (1 - fz)),
            (ix, iy1, iz, (1 - fx) * fy * (1 - fz)),
            (ix1, iy1, iz, fx * fy * (1 - fz)),
            (ix, iy, iz1, (1 - fx) * (1 - fy) * fz),
            (ix1, iy, iz1, fx * (1 - fy) * fz),
            (ix, iy1, iz1, (1 - fx) * fy * fz),
            (ix1, iy1, iz1, fx * fy * fz),
        ):
            idx = (dx[valid] * f + dy[valid]) * f + dz[valid]
            out += np.bincount(idx, weights=w * ww[valid], minlength=flat_size)
        return out.reshape((f, f, f)).astype(np.float32)

    def deposit(self, pos: np.ndarray, mass: np.ndarray) -> np.ndarray:
        """Cloud-in-Cell deposit; returns the full (possibly padded) grid."""
        pos = np.ascontiguousarray(pos, dtype=np.float32)
        mass = np.ascontiguousarray(mass, dtype=np.float32)
        if HAVE_NUMBA:
            set_num_threads(self.nthreads)
            bufs = np.zeros((self.nthreads, self.nfft, self.nfft, self.nfft),
                            dtype=np.float32)
            ox, oy, oz = (float(self.origin[0]), float(self.origin[1]),
                          float(self.origin[2]))
            if self.periodic:
                _deposit_periodic_numba(pos, mass, ox, oy, oz,
                                        1.0 / self.h, self.n, bufs)
            else:
                _deposit_isolated_numba(pos, mass, ox, oy, oz, 1.0 / self.h,
                                        self.n, self.offset, bufs)
            density = bufs.sum(axis=0)
            return density
        return self._deposit_numpy(pos, mass)

    # -- force evaluation --------------------------------------------------
    def _gradient(self, phi: np.ndarray) -> np.ndarray:
        if self.periodic:
            out = np.empty((self.n, self.n, self.n, 3), dtype=np.float32)
            _gradient_periodic_numba(phi, 1.0 / (12.0 * self.h), out)
            return out
        m = phi.shape[0]
        out = np.empty((m, m, m, 3), dtype=np.float32)
        _gradient_isolated_numba(phi, 1.0 / self.h, 1.0 / (12.0 * self.h), out)
        return out

    def accelerations(self, pos: np.ndarray, mass: np.ndarray) -> np.ndarray:
        """Return the PM acceleration for every particle (N, 3) float32."""
        density = self.deposit(pos, mass)
        rho_k = _rfftn(density)
        if self.periodic:
            phi = _irfftn(rho_k * self.kernel * self.G,
                          (self.n, self.n, self.n)).astype(np.float32)
            acc = self._gradient(phi)
            self.last_phi = phi
            self.last_density = density
            return self._interp_periodic(acc, pos)
        phi = _irfftn(rho_k * self.kernel,
                      (self.nfft, self.nfft, self.nfft)).astype(np.float32)
        ext = self.n + 4
        lo = self.offset - 2
        phi_ext = np.ascontiguousarray(phi[lo:lo + ext, lo:lo + ext, lo:lo + ext])
        self.last_phi_ext = phi_ext
        self.last_origin_ext = self.origin - 2.0 * self.h
        central = np.ascontiguousarray(
            density[self.offset:self.offset + self.n,
                    self.offset:self.offset + self.n,
                    self.offset:self.offset + self.n])
        self.last_density = central
        acc_ext = self._gradient(phi_ext)
        return self._interp_isolated(acc_ext, pos)

    def _interp_periodic(self, acc: np.ndarray, pos: np.ndarray) -> np.ndarray:
        out = np.zeros((pos.shape[0], 3), dtype=np.float32)
        if HAVE_NUMBA:
            set_num_threads(self.nthreads)
            ox, oy, oz = map(float, self.origin)
            _interp_periodic_numba(acc, np.ascontiguousarray(pos, np.float32),
                                   ox, oy, oz, 1.0 / self.h, self.n, out)
            return out
        n = self.n
        gx = (pos[:, 0] - self.origin[0]) / self.h - 0.5
        gy = (pos[:, 1] - self.origin[1]) / self.h - 0.5
        gz = (pos[:, 2] - self.origin[2]) / self.h - 0.5
        ix = np.floor(gx).astype(np.int64)
        iy = np.floor(gy).astype(np.int64)
        iz = np.floor(gz).astype(np.int64)
        fx = gx - ix
        fy = gy - iy
        fz = gz - iz
        ix %= n
        iy %= n
        iz %= n
        ix1 = (ix + 1) % n
        iy1 = (iy + 1) % n
        iz1 = (iz + 1) % n
        for c in range(3):
            out[:, c] = (
                acc[ix, iy, iz, c] * (1 - fx) * (1 - fy) * (1 - fz)
                + acc[ix1, iy, iz, c] * fx * (1 - fy) * (1 - fz)
                + acc[ix, iy1, iz, c] * (1 - fx) * fy * (1 - fz)
                + acc[ix1, iy1, iz, c] * fx * fy * (1 - fz)
                + acc[ix, iy, iz1, c] * (1 - fx) * (1 - fy) * fz
                + acc[ix1, iy, iz1, c] * fx * (1 - fy) * fz
                + acc[ix, iy1, iz1, c] * (1 - fx) * fy * fz
                + acc[ix1, iy1, iz1, c] * fx * fy * fz)
        return out.astype(np.float32)

    def _interp_isolated(self, acc: np.ndarray, pos: np.ndarray) -> np.ndarray:
        out = np.zeros((pos.shape[0], 3), dtype=np.float32)
        m = acc.shape[0]
        ox, oy, oz = map(float, self.last_origin_ext)
        if HAVE_NUMBA:
            set_num_threads(self.nthreads)
            vx0, vy0, vz0 = map(float, self.origin)
            vx1 = vy1 = vz1 = float(self.origin[0] + self.box)
            vy1 = float(self.origin[1] + self.box)
            vz1 = float(self.origin[2] + self.box)
            _interp_isolated_numba(acc, np.ascontiguousarray(pos, np.float32),
                                   ox, oy, oz, 1.0 / self.h, m,
                                   vx0, vy0, vz0, vx1, vy1, vz1, out)
            return out
        # NumPy fallback: same kernel with fancy indexing.
        p = np.asarray(pos)
        valid = ((p[:, 0] >= self.origin[0]) & (p[:, 0] < self.origin[0] + self.box)
                 & (p[:, 1] >= self.origin[1]) & (p[:, 1] < self.origin[1] + self.box)
                 & (p[:, 2] >= self.origin[2]) & (p[:, 2] < self.origin[2] + self.box))
        pv = p[valid]
        gx = (pv[:, 0] - ox) / self.h - 0.5
        gy = (pv[:, 1] - oy) / self.h - 0.5
        gz = (pv[:, 2] - oz) / self.h - 0.5
        ix = np.floor(gx).astype(np.int64)
        iy = np.floor(gy).astype(np.int64)
        iz = np.floor(gz).astype(np.int64)
        fx = gx - ix
        fy = gy - iy
        fz = gz - iz
        ix = np.clip(ix, 0, m - 2)
        iy = np.clip(iy, 0, m - 2)
        iz = np.clip(iz, 0, m - 2)
        ix1 = ix + 1
        iy1 = iy + 1
        iz1 = iz + 1
        vals = np.zeros((pv.shape[0], 3), dtype=np.float32)
        for c in range(3):
            vals[:, c] = (
                acc[ix, iy, iz, c] * (1 - fx) * (1 - fy) * (1 - fz)
                + acc[ix1, iy, iz, c] * fx * (1 - fy) * (1 - fz)
                + acc[ix, iy1, iz, c] * (1 - fx) * fy * (1 - fz)
                + acc[ix1, iy1, iz, c] * fx * fy * (1 - fz)
                + acc[ix, iy, iz1, c] * (1 - fx) * (1 - fy) * fz
                + acc[ix1, iy, iz1, c] * fx * (1 - fy) * fz
                + acc[ix, iy1, iz1, c] * (1 - fx) * fy * fz
                + acc[ix1, iy1, iz1, c] * fx * fy * fz)
        out[valid] = vals
        return out

    # -- potential ---------------------------------------------------------
    def potential_at(self, pos: np.ndarray) -> np.ndarray:
        """Interpolate the last computed potential to particle positions."""
        out = np.zeros(pos.shape[0], dtype=np.float32)
        if self.periodic:
            if self.last_phi is None:
                return out
            if HAVE_NUMBA:
                set_num_threads(self.nthreads)
                ox, oy, oz = map(float, self.origin)
                _interp_scalar_periodic_numba(
                    self.last_phi, np.ascontiguousarray(pos, np.float32),
                    ox, oy, oz, 1.0 / self.h, self.n, out)
                return out
            gx = (pos[:, 0] - self.origin[0]) / self.h - 0.5
            gy = (pos[:, 1] - self.origin[1]) / self.h - 0.5
            gz = (pos[:, 2] - self.origin[2]) / self.h - 0.5
            ix = np.floor(gx).astype(np.int64) % self.n
            iy = np.floor(gy).astype(np.int64) % self.n
            iz = np.floor(gz).astype(np.int64) % self.n
            fx = gx - np.floor(gx)
            fy = gy - np.floor(gy)
            fz = gz - np.floor(gz)
            ix1 = (ix + 1) % self.n
            iy1 = (iy + 1) % self.n
            iz1 = (iz + 1) % self.n
            f = self.last_phi
            out = (f[ix, iy, iz] * (1 - fx) * (1 - fy) * (1 - fz)
                   + f[ix1, iy, iz] * fx * (1 - fy) * (1 - fz)
                   + f[ix, iy1, iz] * (1 - fx) * fy * (1 - fz)
                   + f[ix1, iy1, iz] * fx * fy * (1 - fz)
                   + f[ix, iy, iz1] * (1 - fx) * (1 - fy) * fz
                   + f[ix1, iy, iz1] * fx * (1 - fy) * fz
                   + f[ix, iy1, iz1] * (1 - fx) * fy * fz
                   + f[ix1, iy1, iz1] * fx * fy * fz).astype(np.float32)
            return out
        if self.last_phi_ext is None:
            return out
        m = self.last_phi_ext.shape[0]
        ox, oy, oz = map(float, self.last_origin_ext)
        vx0, vy0, vz0 = map(float, self.origin)
        vx1 = float(self.origin[0] + self.box)
        vy1 = float(self.origin[1] + self.box)
        vz1 = float(self.origin[2] + self.box)
        if HAVE_NUMBA:
            set_num_threads(self.nthreads)
            _interp_scalar_isolated_numba(
                self.last_phi_ext, np.ascontiguousarray(pos, np.float32),
                ox, oy, oz, 1.0 / self.h, m,
                vx0, vy0, vz0, vx1, vy1, vz1, out)
            return out
        p = np.asarray(pos)
        valid = ((p[:, 0] >= vx0) & (p[:, 0] < vx1)
                 & (p[:, 1] >= vy0) & (p[:, 1] < vy1)
                 & (p[:, 2] >= vz0) & (p[:, 2] < vz1))
        pv = p[valid]
        gx = (pv[:, 0] - ox) / self.h - 0.5
        gy = (pv[:, 1] - oy) / self.h - 0.5
        gz = (pv[:, 2] - oz) / self.h - 0.5
        ix = np.clip(np.floor(gx).astype(np.int64), 0, m - 2)
        iy = np.clip(np.floor(gy).astype(np.int64), 0, m - 2)
        iz = np.clip(np.floor(gz).astype(np.int64), 0, m - 2)
        fx = gx - np.floor(gx)
        fy = gy - np.floor(gy)
        fz = gz - np.floor(gz)
        f = self.last_phi_ext
        ix1, iy1, iz1 = ix + 1, iy + 1, iz + 1
        vals = (f[ix, iy, iz] * (1 - fx) * (1 - fy) * (1 - fz)
                + f[ix1, iy, iz] * fx * (1 - fy) * (1 - fz)
                + f[ix, iy1, iz] * (1 - fx) * fy * (1 - fz)
                + f[ix1, iy1, iz] * fx * fy * (1 - fz)
                + f[ix, iy, iz1] * (1 - fx) * (1 - fy) * fz
                + f[ix1, iy, iz1] * fx * (1 - fy) * fz
                + f[ix, iy1, iz1] * (1 - fx) * fy * fz
                + f[ix1, iy1, iz1] * fx * fy * fz)
        out[valid] = vals.astype(np.float32)
        return out

    def total_mass_on_grid(self, density: np.ndarray) -> float:
        """Convenience helper for tests and diagnostics."""
        return float(np.sum(density, dtype=np.float64))


def make_gravity(n: int, box: float, periodic: bool, G: float,
                 center: np.ndarray | None = None,
                 smoothing_cells: float = config.SMOOTHING_CELLS
                 ) -> PMGravity:
    """Factory used by both scenes; prints a warning when numba is missing."""
    if not HAVE_NUMBA:
        warnings.warn(
            "numba is not available: using the pure-NumPy fallback. "
            "Particle counts and grid sizes are reduced automatically.",
            RuntimeWarning, stacklevel=2)
    return PMGravity(n, box, periodic, G, center, smoothing_cells)

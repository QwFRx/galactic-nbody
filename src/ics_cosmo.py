"""Cosmological initial conditions: linear power spectrum and Zel'dovich field."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

import config
from . import units

try:  # numba is optional but strongly recommended
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised on installations without numba
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda fn: fn

    prange = range


# ---------------------------------------------------------------------------
# Linear theory
# ---------------------------------------------------------------------------
def eisenstein_hu_power(k: np.ndarray, omega_m: float = config.OMEGA_M,
                        omega_b: float = 0.045,
                        h: float = config.H_DIMENSIONLESS,
                        ns: float = config.NS) -> np.ndarray:
    """Eisenstein & Hu (1998) no-wiggle transfer function.

    P(k) is returned in arbitrary units; the caller normalises to sigma_8.
    k is in h/Mpc.  The shape parameter is
    Gamma = Omega_m h exp(-Omega_b (1 + sqrt(2h)/Omega_m)).
    """
    kk = np.maximum(np.asarray(k, dtype=np.float64), 1.0e-8)
    gamma = omega_m * h * np.exp(-omega_b * (1.0 + np.sqrt(2.0 * h) / omega_m))
    q = kk / gamma
    l0 = np.log(2.0 * np.e + 1.8 * q)
    c0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    transfer = l0 / (l0 + c0 * q * q)
    return kk ** ns * transfer ** 2


def sigma_r(power_fn, radius: float, kmin: float = 1.0e-4,
            kmax: float = 1.0e2, nk: int = 4096) -> float:
    """RMS density contrast in a top-hat sphere of radius R (Mpc/h)."""
    lk = np.linspace(np.log(kmin), np.log(kmax), nk)
    k = np.exp(lk)
    p = power_fn(k)
    x = k * radius
    window = np.where(x > 1.0e-4,
                      3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3,
                      1.0 - x * x / 10.0)
    delta2 = k ** 3 * p / (2.0 * np.pi ** 2) * window ** 2
    return float(np.sqrt(np.trapezoid(delta2, lk)))


def linear_power(k: np.ndarray, omega_m: float = config.OMEGA_M,
                 omega_b: float = 0.045, h: float = config.H_DIMENSIONLESS,
                 ns: float = config.NS,
                 sigma8: float = config.SIGMA8) -> np.ndarray:
    """Eisenstein-Hu power spectrum normalised to sigma_8.

    The returned P(k) has units (Mpc/h)^3 with k in h/Mpc.
    """
    raw = lambda kk: eisenstein_hu_power(kk, omega_m, omega_b, h, ns)
    s8 = sigma_r(raw, 8.0, kmin=1.0e-4, kmax=1.0e2)
    norm = (sigma8 / s8) ** 2
    return raw(k) * norm


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
@dataclass
class CosmicIC:
    pos: np.ndarray        # (Np, 3) float32, comoving Mpc/h
    pvel: np.ndarray       # (Np, 3) float32, a^2 dx/dt in 100 km/s
    mass: np.ndarray       # (Np,) float32, Omega_m / Np
    a_start: float
    delta: np.ndarray      # (G, G, G) float32 density contrast on the PM grid
    box: float
    grid: int
    particle_grid: int


@njit(parallel=True, cache=True, fastmath=True, nogil=True)
def _sample_field_numba(field: np.ndarray, qpos: np.ndarray,
                        inv_h: float, n: int, out: np.ndarray) -> None:
    """Periodic trilinear interpolation of a 3-component grid field."""
    m = qpos.shape[0]
    for p in prange(m):
        gx = qpos[p, 0] * inv_h - 0.5
        gy = qpos[p, 1] * inv_h - 0.5
        gz = qpos[p, 2] * inv_h - 0.5
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
                field[c, i0, j0, k0] * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
                + field[c, i1, j0, k0] * fx * (1.0 - fy) * (1.0 - fz)
                + field[c, i0, j1, k0] * (1.0 - fx) * fy * (1.0 - fz)
                + field[c, i1, j1, k0] * fx * fy * (1.0 - fz)
                + field[c, i0, j0, k1] * (1.0 - fx) * (1.0 - fy) * fz
                + field[c, i1, j0, k1] * fx * (1.0 - fy) * fz
                + field[c, i0, j1, k1] * (1.0 - fx) * fy * fz
                + field[c, i1, j1, k1] * fx * fy * fz)


def _sample_field_numpy(field: np.ndarray, qpos: np.ndarray,
                        inv_h: float, n: int) -> np.ndarray:
    from scipy.ndimage import map_coordinates
    coords = np.empty((3, qpos.shape[0]), dtype=np.float32)
    for c in range(3):
        coords[c] = qpos[:, c] * inv_h - 0.5
    return np.stack([map_coordinates(field[c], coords, order=1, mode="wrap")
                     for c in range(3)], axis=1).astype(np.float32)


def make_cosmic_ics(particles: int, grid: int, box: float,
                    a_start: float = 0.02, seed: int = 12345,
                    omega_m: float = config.OMEGA_M,
                    omega_l: float = config.OMEGA_L,
                    h: float = config.H_DIMENSIONLESS,
                    ns: float = config.NS,
                    sigma8: float = config.SIGMA8) -> CosmicIC:
    """Build a Zel'dovich realisation of a flat LambdaCDM density field.

    The linear field is generated on the PM grid and normalised to sigma_8.
    Particles live on an approximately cubic lattice Np^3 <= particles and
    the displacement field is trilinearly sampled at those positions.
    """
    g = int(grid)
    lbox = float(box)
    a0 = float(a_start)
    if not (0.0 < a0 <= 1.0):
        raise ValueError("a_start must lie in (0, 1]")
    rng = np.random.default_rng(int(seed))

    white = rng.standard_normal((g, g, g), dtype=np.float64)
    white_k = np.fft.fftn(white)
    del white
    k1 = 2.0 * np.pi * np.fft.fftfreq(g, d=lbox / g)
    k2 = (k1[:, None, None] ** 2 + k1[None, :, None] ** 2
          + k1[None, None, :] ** 2)
    kk = np.sqrt(k2)
    power = linear_power(kk, omega_m, 0.045, h, ns, sigma8)
    nbar_grid = g ** 3 / lbox ** 3
    delta_k = white_k * np.sqrt(power * nbar_grid)
    del white_k, power
    delta_k[0, 0, 0] = 0.0 + 0.0j

    # Zel'dovich displacement: Psi_k = i k delta_k / k^2.
    safe_k2 = np.where(k2 > 0.0, k2, 1.0)
    psi = np.empty((3, g, g, g), dtype=np.float32)
    for axis in range(3):
        comp = np.zeros((g, g, g), dtype=np.complex128)
        if axis == 0:
            comp = 1j * k1[:, None, None] * delta_k / safe_k2
        elif axis == 1:
            comp = 1j * k1[None, :, None] * delta_k / safe_k2
        else:
            comp = 1j * k1[None, None, :] * delta_k / safe_k2
        comp[0, 0, 0] = 0.0
        psi[axis] = np.fft.ifftn(comp).real.astype(np.float32)
        del comp

    delta = np.fft.ifftn(delta_k).real.astype(np.float32)
    del delta_k, k2, kk

    np_grid = max(8, int(round(float(particles) ** (1.0 / 3.0))))
    np_grid = min(np_grid, 512)
    count = np_grid ** 3
    spacing = lbox / np_grid
    axis = (np.arange(np_grid, dtype=np.float64) + 0.5) * spacing
    qpos = np.empty((count, 3), dtype=np.float32)
    qpos[:, 0] = np.repeat(axis, np_grid * np_grid)
    qpos[:, 1] = np.tile(np.repeat(axis, np_grid), np_grid)
    qpos[:, 2] = np.tile(axis, np_grid * np_grid)

    if HAVE_NUMBA:
        psi_sample = np.empty((count, 3), dtype=np.float32)
        _sample_field_numba(psi, qpos, 1.0 / spacing, g, psi_sample)
    else:
        psi_sample = _sample_field_numpy(psi, qpos, 1.0 / spacing, g)

    d_start = float(units.growth_factor(a0, omega_m, omega_l))
    f_start = float(units.growth_rate(a0, omega_m, omega_l))
    e_start = float(units.hubble_e(a0, omega_m, omega_l))

    pos = qpos + np.float32(d_start) * psi_sample
    pos %= lbox
    pos = pos.astype(np.float32)

    # p = a^2 dx/dt = a^2 H f D Psi in code units (H0 = 1).
    pvel = (a0 ** 2 * e_start * f_start * d_start * psi_sample).astype(np.float32)

    mass = np.full(count, omega_m / count, dtype=np.float32)
    return CosmicIC(pos=pos, pvel=pvel, mass=mass, a_start=a0, delta=delta,
                    box=lbox, grid=g, particle_grid=np_grid)


def growth_from_integral(a: float, omega_m: float = config.OMEGA_M,
                         omega_l: float = config.OMEGA_L) -> float:
    """Unnormalised growth factor (helper kept here for the tests)."""
    return float(units.growth_factor(a, omega_m, omega_l, normalized=False))

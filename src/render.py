"""Software point-splatter renderer with bloom, tone mapping and a 3D box.

The hot loops are numba kernels.  Each thread accumulates into its own
image buffer and the buffers are merged afterwards, so there are no data
races.  A pure-NumPy fallback is provided for installations without numba.
"""

from __future__ import annotations

import os

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter

import config

try:
    from numba import njit, prange, get_num_threads, set_num_threads, get_thread_id
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda fn: fn

    prange = range

    def get_num_threads() -> int:
        return 1

    def set_num_threads(n: int) -> None:
        return None

    def get_thread_id() -> int:
        return 0


def _build_lut(anchors: np.ndarray, size: int = 256) -> np.ndarray:
    t = np.linspace(0.0, 1.0, size)
    lut = np.empty((size, 3), dtype=np.float32)
    for channel in range(3):
        lut[:, channel] = np.interp(t, anchors[:, 0], anchors[:, channel + 1])
    return lut


INFERNO_LUT = _build_lut(config.INFERNO_ANCHORS)
WEB_LUT = _build_lut(config.WEB_ANCHORS)


@njit(parallel=True, cache=True, fastmath=True, nogil=True)
def _splat_points_numba(pos: np.ndarray, rgb: np.ndarray, weight: np.ndarray,
                        eye: np.ndarray, fwd: np.ndarray, rgt: np.ndarray,
                        up: np.ndarray, focal: float, width: int, height: int,
                        bufs: np.ndarray) -> None:
    count = pos.shape[0]
    for i in prange(count):
        dx = pos[i, 0] - eye[0]
        dy = pos[i, 1] - eye[1]
        dz = pos[i, 2] - eye[2]
        depth = dx * fwd[0] + dy * fwd[1] + dz * fwd[2]
        if depth <= 1.0e-3:
            continue
        inv = 1.0 / depth
        sx = 0.5 * width + (dx * rgt[0] + dy * rgt[1] + dz * rgt[2]) * inv * focal
        sy = 0.5 * height - (dx * up[0] + dy * up[1] + dz * up[2]) * inv * focal
        if sx < 0.0 or sx >= width - 1.0 or sy < 0.0 or sy >= height - 1.0:
            continue
        ix = int(sx)
        iy = int(sy)
        fx = sx - ix
        fy = sy - iy
        tid = get_thread_id()
        base = iy * width + ix
        w = weight[i]
        w00 = w * (1.0 - fx) * (1.0 - fy)
        w10 = w * fx * (1.0 - fy)
        w01 = w * (1.0 - fx) * fy
        w11 = w * fx * fy
        r = rgb[i, 0]
        g = rgb[i, 1]
        b = rgb[i, 2]
        bufs[tid, base, 0] += r * w00
        bufs[tid, base, 1] += g * w00
        bufs[tid, base, 2] += b * w00
        bufs[tid, base + 1, 0] += r * w10
        bufs[tid, base + 1, 1] += g * w10
        bufs[tid, base + 1, 2] += b * w10
        bufs[tid, base + width, 0] += r * w01
        bufs[tid, base + width, 1] += g * w01
        bufs[tid, base + width, 2] += b * w01
        bufs[tid, base + width + 1, 0] += r * w11
        bufs[tid, base + width + 1, 1] += g * w11
        bufs[tid, base + width + 1, 2] += b * w11


@njit(parallel=True, cache=True, fastmath=True, nogil=True)
def _splat_scalar_numba(pos: np.ndarray, weight: np.ndarray,
                        eye: np.ndarray, fwd: np.ndarray, rgt: np.ndarray,
                        up: np.ndarray, focal: float, width: int, height: int,
                        depth_center: float, half_slab: float,
                        bufs: np.ndarray) -> None:
    count = pos.shape[0]
    for i in prange(count):
        dx = pos[i, 0] - eye[0]
        dy = pos[i, 1] - eye[1]
        dz = pos[i, 2] - eye[2]
        depth = dx * fwd[0] + dy * fwd[1] + dz * fwd[2]
        if depth <= 1.0e-3 or abs(depth - depth_center) > half_slab:
            continue
        inv = 1.0 / depth
        sx = 0.5 * width + (dx * rgt[0] + dy * rgt[1] + dz * rgt[2]) * inv * focal
        sy = 0.5 * height - (dx * up[0] + dy * up[1] + dz * up[2]) * inv * focal
        if sx < 0.0 or sx >= width - 1.0 or sy < 0.0 or sy >= height - 1.0:
            continue
        ix = int(sx)
        iy = int(sy)
        fx = sx - ix
        fy = sy - iy
        tid = get_thread_id()
        w = weight[i]
        base = iy * width + ix
        bufs[tid, base] += w * (1.0 - fx) * (1.0 - fy)
        bufs[tid, base + 1] += w * fx * (1.0 - fy)
        bufs[tid, base + width] += w * (1.0 - fx) * fy
        bufs[tid, base + width + 1] += w * fx * fy


class Renderer:
    """Render snapshots to RGB arrays; also handles PNG export."""

    def __init__(self, width: int, height: int) -> None:
        self.width = int(width)
        self.height = int(height)
        cpu = os.cpu_count() or 1
        per_thread = self.width * self.height * 3 * 4
        self.max_threads = max(1, min(cpu, int(1_200_000_000 / per_thread)))
        if not HAVE_NUMBA:
            self.max_threads = 1

    # -- public API --------------------------------------------------------
    def render(self, snap: dict, camera, color_mode: int = 0,
               show_dm: bool = False, exposure: float = 4.0,
               slab: float = 0.35, supersample: int = 1) -> np.ndarray:
        """Render a snapshot and return an (H, W, 3) uint8 RGB array."""
        if supersample > 1:
            big = Renderer(self.width * supersample, self.height * supersample)
            image = big.render(snap, camera, color_mode, show_dm, exposure,
                               slab, supersample=1)
            image = image.reshape(self.height, supersample,
                                  self.width, supersample, 3).mean(axis=(1, 3))
            return image.astype(np.uint8)
        if snap["scene"] == "cosmic":
            return self._render_cosmic(snap, camera, exposure, slab)
        return self._render_collision(snap, camera, color_mode, show_dm, exposure)

    # -- collision ---------------------------------------------------------
    def _render_collision(self, snap: dict, camera, color_mode: int,
                          show_dm: bool, exposure: float) -> np.ndarray:
        pos = np.asarray(snap["pos"], dtype=np.float32)
        comp = snap["comp"]
        gal = snap["gal"]
        speed = np.asarray(snap["speed"], dtype=np.float32)
        if not show_dm:
            keep = comp != 0
            if not np.all(keep):
                pos = pos[keep]
                comp = comp[keep]
                gal = gal[keep]
                speed = speed[keep]
        rgb, weight = self._collision_style(snap, pos, comp, gal, speed,
                                            color_mode)
        img = self._splat_points(pos, rgb, weight, camera)
        return self._finish_points(img, exposure)

    def _collision_style(self, snap: dict, pos: np.ndarray, comp: np.ndarray,
                         gal: np.ndarray, speed: np.ndarray, mode: int
                         ) -> tuple[np.ndarray, np.ndarray]:
        count = comp.shape[0]
        component_weight = np.array([0.22, 1.0, 0.85], dtype=np.float32)
        weight = component_weight[comp] * np.float32(0.055)
        if mode == 0:
            rgb = config.GALAXY_COLORS[gal].copy()
        elif mode == 1:
            rgb = config.COMPONENT_COLORS[comp].copy()
        elif mode == 2:
            dens = snap.get("density")
            if dens is None:
                t = np.zeros(count, dtype=np.float32)
            else:
                mean = float(np.mean(dens)) + 1.0e-30
                sampled = self._sample_density(dens, pos, snap["center"],
                                               snap["box"], periodic=False)
                ratio = sampled / mean
                t = np.log10(1.0 + np.maximum(ratio, 0.0)) / np.log10(51.0)
                t = np.clip(t, 0.0, 1.0).astype(np.float32)
            rgb = INFERNO_LUT[(t * 255).astype(np.int32)]
            weight *= (0.35 + 1.4 * t)
        else:
            scale = float(np.percentile(speed, 99.0)) + 1.0e-6
            t = np.clip(speed / scale, 0.0, 1.0)
            rgb = INFERNO_LUT[(t * 255).astype(np.int32)]
        return np.ascontiguousarray(rgb, dtype=np.float32), weight.astype(np.float32)

    # -- cosmic ------------------------------------------------------------
    def _render_cosmic(self, snap: dict, camera, exposure: float,
                       slab: float) -> np.ndarray:
        box = float(snap["box"])
        center = np.asarray(snap["center"], dtype=np.float64)
        camera.target = center
        depth_center = float(np.dot(center - camera.eye(), camera.basis()[0]))
        half_slab = max(0.02, slab) * 0.5 * box
        col = self._splat_scalar(np.asarray(snap["pos"], dtype=np.float32),
                                 np.asarray(snap["mass"], dtype=np.float32),
                                 camera, depth_center, half_slab)
        # A small blur suppresses the shot noise of individual particles and
        # turns the 2x2 CIC footprint into a continuous density field.
        col = gaussian_filter(col, sigma=1.6)
        total = float(col.sum())
        occupied = max(int(np.count_nonzero(col > 0.0)), 1)
        mean = total / occupied
        ratio = col / max(mean, 1.0e-30)
        t = (np.log10(np.maximum(ratio, 1.0e-6)) + 1.5) / 3.0
        t = np.clip(np.power(np.clip(t, 0.0, 1.0), 1.3)
                    * (exposure / 5.0), 0.0, 1.0)
        rgb = WEB_LUT[(t * 255).astype(np.int32)]
        rgb *= (0.25 + 0.75 * t)[..., None]
        bright = np.maximum(rgb - 0.55, 0.0)
        bloom = np.stack([gaussian_filter(bright[..., c], sigma=2.5)
                          for c in range(3)], axis=-1)
        rgb = np.clip(rgb + 0.7 * bloom, 0.0, 1.0)
        rgb = self._vignette(rgb)
        frame = (rgb * 255.0).astype(np.uint8)
        frame = self._draw_box(frame, camera, center, box)
        return frame

    # -- low-level splatting ----------------------------------------------
    def _splat_points(self, pos: np.ndarray, rgb: np.ndarray,
                      weight: np.ndarray, camera) -> np.ndarray:
        w, h = self.width, self.height
        eye, fwd, rgt, focal = camera.projection(w, h)
        eye = np.ascontiguousarray(eye, dtype=np.float32)
        fwd = np.ascontiguousarray(fwd, dtype=np.float32)
        rgt = np.ascontiguousarray(rgt, dtype=np.float32)
        up = np.ascontiguousarray(np.cross(rgt.astype(np.float64),
                                           fwd.astype(np.float64)), dtype=np.float32)
        if HAVE_NUMBA:
            nt = self.max_threads
            set_num_threads(nt)
            bufs = np.zeros((nt, w * h, 3), dtype=np.float32)
            _splat_points_numba(np.ascontiguousarray(pos, dtype=np.float32),
                                np.ascontiguousarray(rgb, dtype=np.float32),
                                np.ascontiguousarray(weight, dtype=np.float32),
                                eye, fwd, rgt, up, float(focal), w, h, bufs)
            return bufs.sum(axis=0).reshape((h, w, 3))
        return self._splat_points_numpy(pos, rgb, weight, camera)

    def _splat_points_numpy(self, pos: np.ndarray, rgb: np.ndarray,
                            weight: np.ndarray, camera) -> np.ndarray:
        w, h = self.width, self.height
        eye, fwd, rgt, focal = camera.projection(w, h)
        up = np.cross(rgt, fwd)
        rel = np.asarray(pos, dtype=np.float64) - eye
        depth = rel @ fwd
        valid = depth > 1.0e-3
        inv = 1.0 / np.maximum(depth, 1.0e-9)
        sx = w * 0.5 + (rel @ rgt) * inv * focal
        sy = h * 0.5 - (rel @ up) * inv * focal
        ix = np.floor(sx).astype(np.int64)
        iy = np.floor(sy).astype(np.int64)
        valid &= (ix >= 0) & (ix < w - 1) & (iy >= 0) & (iy < h - 1) & (depth > 0)
        ix = ix[valid]
        iy = iy[valid]
        fx = sx[valid] - ix
        fy = sy[valid] - iy
        src = np.asarray(rgb, dtype=np.float64)[valid] * weight[valid, None]
        out = np.zeros((h * w, 3), dtype=np.float64)
        for ox, oy, ww in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)),
                           (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
            idx = (iy + oy) * w + (ix + ox)
            for c in range(3):
                out[:, c] += np.bincount(idx, weights=src[:, c] * ww,
                                         minlength=h * w)
        return out.reshape((h, w, 3)).astype(np.float32)

    def _splat_scalar(self, pos: np.ndarray, mass: np.ndarray, camera,
                      depth_center: float, half_slab: float) -> np.ndarray:
        w, h = self.width, self.height
        eye, fwd, rgt, focal = camera.projection(w, h)
        eye = np.ascontiguousarray(eye, dtype=np.float32)
        fwd = np.ascontiguousarray(fwd, dtype=np.float32)
        rgt = np.ascontiguousarray(rgt, dtype=np.float32)
        up = np.ascontiguousarray(np.cross(rgt.astype(np.float64),
                                           fwd.astype(np.float64)), dtype=np.float32)
        if HAVE_NUMBA:
            nt = self.max_threads
            set_num_threads(nt)
            bufs = np.zeros((nt, w * h), dtype=np.float32)
            _splat_scalar_numba(np.ascontiguousarray(pos, dtype=np.float32),
                                np.ascontiguousarray(mass, dtype=np.float32),
                                eye, fwd, rgt, up, float(focal), w, h,
                                float(depth_center), float(half_slab), bufs)
            return bufs.sum(axis=0).reshape((h, w))
        rel = np.asarray(pos, dtype=np.float64) - eye
        depth = rel @ fwd
        valid = (depth > 1.0e-3) & (np.abs(depth - depth_center) <= half_slab)
        inv = 1.0 / np.maximum(depth, 1.0e-9)
        sx = w * 0.5 + (rel @ rgt) * inv * focal
        sy = h * 0.5 - (rel @ up) * inv * focal
        ix = np.floor(sx).astype(np.int64)
        iy = np.floor(sy).astype(np.int64)
        valid &= (ix >= 0) & (ix < w - 1) & (iy >= 0) & (iy < h - 1)
        ix = ix[valid]
        iy = iy[valid]
        fx = sx[valid] - ix
        fy = sy[valid] - iy
        m = np.asarray(mass, dtype=np.float64)[valid]
        out = np.zeros(h * w, dtype=np.float64)
        for ox, oy, ww in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)),
                           (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
            idx = (iy + oy) * w + (ix + ox)
            out += np.bincount(idx, weights=m * ww, minlength=h * w)
        return out.reshape((h, w)).astype(np.float32)

    # -- helper fields -----------------------------------------------------
    @staticmethod
    def _sample_density(density: np.ndarray, pos: np.ndarray,
                        center: np.ndarray, box: float,
                        periodic: bool) -> np.ndarray:
        n = density.shape[0]
        h = box / n
        origin = np.asarray(center, dtype=np.float32) - np.float32(0.5 * box)
        rel = (np.asarray(pos, dtype=np.float32) - origin) / np.float32(h)
        idx = np.floor(rel).astype(np.int64)
        if periodic:
            idx %= n
        else:
            inside = np.all((idx >= 0) & (idx < n), axis=1)
            idx = np.clip(idx, 0, n - 1)
        vals = density[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.float32)
        if not periodic:
            vals[~inside] = 0.0
        return vals

    @staticmethod
    def _vignette(rgb: np.ndarray) -> np.ndarray:
        h, w = rgb.shape[:2]
        yy, xx = np.mgrid[0:h, 0:w]
        nx = (xx - 0.5 * w) / (0.5 * w)
        ny = (yy - 0.5 * h) / (0.5 * h)
        r2 = nx * nx + ny * ny
        factor = np.clip(1.0 - 0.28 * r2, 0.0, 1.0)[..., None]
        return rgb * factor

    def _finish_points(self, img: np.ndarray, exposure: float) -> np.ndarray:
        # A sub-pixel blur turns isolated particles into a continuous glow;
        # without it a 250k preview looks like a field of coloured dots.
        img = np.stack([gaussian_filter(img[..., c], sigma=0.8)
                        for c in range(3)], axis=-1)
        positive = img[img > 0.0]
        if positive.size:
            scale = float(np.percentile(positive, 99.5))
        else:
            scale = 1.0
        img = img / (scale + 1.0e-9)
        img = img * float(exposure)
        bright = np.maximum(img - 2.5, 0.0)
        bloom = np.stack([gaussian_filter(bright[..., c], sigma=2.0)
                          for c in range(3)], axis=-1)
        img = img + 0.55 * bloom
        img = np.arcsinh(img) / max(np.arcsinh(exposure), 1.0e-6)
        img = np.clip(img, 0.0, 1.0)
        img = self._vignette(img)
        return (img * 255.0).astype(np.uint8)

    def _draw_box(self, frame: np.ndarray, camera, center: np.ndarray,
                  box: float) -> np.ndarray:
        h, w = frame.shape[:2]
        half = 0.5 * box
        corners = np.array([[x, y, z] for x in (-half, half)
                            for y in (-half, half) for z in (-half, half)])
        corners = corners + center
        xy, depth = camera.project(corners, w, h)
        image = Image.fromarray(frame)
        draw = ImageDraw.Draw(image)
        edges = ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                 (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7))
        for a, b in edges:
            if depth[a] > 0 and depth[b] > 0:
                draw.line([tuple(xy[a]), tuple(xy[b])], fill=(80, 100, 135), width=1)
        origin = center - half
        axis_len = 0.16 * box
        ends = np.array([origin + np.array([axis_len, 0, 0]),
                         origin + np.array([0, axis_len, 0]),
                         origin + np.array([0, 0, axis_len])])
        oxy, _ = camera.project(np.vstack((origin, ends)), w, h)
        colors = ((235, 80, 80), (90, 220, 110), (90, 150, 245))
        labels = ("x", "y", "z")
        for i in range(3):
            draw.line([tuple(oxy[0]), tuple(oxy[i + 1])], fill=colors[i], width=2)
            draw.text((float(oxy[i + 1][0]) + 3.0, float(oxy[i + 1][1]) - 4.0),
                      labels[i], fill=colors[i])
        return np.asarray(image)


def save_png(frame: np.ndarray, path: str) -> None:
    """Write an RGB uint8 array to a PNG file, creating parent directories."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    Image.fromarray(frame).save(path)

"""Central configuration and physical constants for the galactic N-body project.

Code units
----------
Collision scene
    length        : kpc
    velocity      : km/s
    mass          : Msun
    time          : kpc / (km/s) = 0.9778 Gyr
    gravitational : G = 4.30091e-6 kpc (km/s)^2 / Msun

Cosmic scene
    length        : Mpc/h (comoving)
    velocity      : 100 km/s
    mass          : arbitrary (particles carry Omega_m / N each)
    time          : (Mpc/h) / (100 km/s) = 14.379 Gyr for h = 0.68
    expansion     : H0 = 100 (km/s)/(Mpc/h) = 1 in code units

The PM force resolution is roughly two cells.  In the default collision
setup this is a few kpc, and in the cosmic scene it is one to two Mpc/h.
There is no gas, star formation or feedback: both scenes are purely
collisionless.  The cosmic scene therefore has no halo-scale resolution
below the PM cell size.
"""

from __future__ import annotations

import importlib.util
import warnings
from dataclasses import dataclass, field

import numpy as np

# ---------------------------------------------------------------------------
# Collision-scene constants
# ---------------------------------------------------------------------------
G_KPC: float = 4.30091e-6            # kpc (km/s)^2 / Msun
KPC_KM: float = 3.0856775814913673e16  # km per kpc
GYR_S: float = 3.15576e16              # s per Gyr
KPC_KMPS_GYR: float = KPC_KM / GYR_S   # 0.97779 Gyr per kpc/(km/s)

MSUN_KG: float = 1.98847e30
KMS_TO_KPC_PER_GYR: float = KPC_KMPS_GYR

# ---------------------------------------------------------------------------
# Cosmic-scene constants
# ---------------------------------------------------------------------------
OMEGA_M: float = 0.3
OMEGA_L: float = 0.7
H_DIMENSIONLESS: float = 0.68
SIGMA8: float = 0.8
NS: float = 0.96
BOX_MPCH: float = 100.0
Z_START: float = 49.0
MPC_KM: float = 3.0856775814913673e19
# One code time unit is (Mpc/h)/(100 km/s).  Because 1 Mpc/h = 3.0857e19/h km,
# this is 9.778 / h Gyr (14.379 Gyr for h = 0.68).
COSMIC_TIME_GYR: float = MPC_KM / (H_DIMENSIONLESS * 100.0 * GYR_S)
COSMIC_VEL_KMS: float = 100.0

# ---------------------------------------------------------------------------
# Galaxy model defaults (Milky-Way-like primary)
# ---------------------------------------------------------------------------
M_HALO: float = 1.0e12
M_DISK: float = 5.0e10
M_BULGE: float = 1.0e10
R_DISK: float = 3.5
Z_DISK: float = 0.35
A_BULGE: float = 0.7
A_HALO: float = 30.0
R_HALO_MAX: float = 100.0    # truncation radius; M(<R_HALO_MAX) is used for the orbit
R_BULGE_MAX: float = 8.0
R_DISK_MAX_FACTOR: float = 5.0
TOOMRE_Q: float = 1.5
ETA: float = 0.2             # leapfrog timestep safety factor
SMOOTHING_CELLS: float = 0.30  # Gaussian PM smoothing length in cells
MIN_DISK_SHARE: float = 0.25  # minimum fraction of particles in the stellar disk

# Discrete collision-box side lengths in kpc.
BOX_LEVELS: tuple[float, ...] = (64.0, 96.0, 128.0, 192.0, 256.0, 384.0,
                                 512.0, 768.0, 1024.0)

# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------
QUALITY: dict[str, dict[str, object]] = {
    "preview": {"particles": 250_000, "grid": 64, "render": (640, 360),
                "window": (1280, 720)},
    "high": {"particles": 2_000_000, "grid": 128, "render": (960, 540),
             "window": (1280, 720)},
    "ultra": {"particles": 8_000_000, "grid": 256, "render": (1280, 720),
              "window": (1600, 900)},
}

MAX_PARTICLES: int = 8_000_000
MAX_GRID: int = 512

# Default number of headless steps when the user does not pass --steps.
# The collision run stops earlier when it reaches ~0.55 Gyr; the cosmic run
# stops exactly at a = 1 (z = 0).  These are safety caps, not fixed counts.
DEFAULT_STEPS: dict[str, int] = {"collision": 800, "cosmic": 400}


@dataclass
class SimConfig:
    """Fully resolved run configuration passed to both scenes and the UI."""

    scene: str = "collision"
    quality: str = "preview"
    particles: int | None = None
    grid: int | None = None
    width: int | None = None
    height: int | None = None
    seed: int = 12345
    box: float | None = None
    z_start: float = Z_START
    mass_ratio: float = 1.0
    pericenter: float = 15.0
    eccentricity: float = 1.0
    separation: float = 80.0
    incl1: float = 0.0
    incl2: float = 60.0
    steps: int | None = None
    render: str | None = None
    frames: str | None = None
    frame_count: int = 0
    analyze: bool = False
    speed: float = 1.0

    # Resolved values (filled by resolve()).
    n_particles: int = field(init=False, default=0)
    n_grid: int = field(init=False, default=0)
    win_w: int = field(init=False, default=1280)
    win_h: int = field(init=False, default=720)
    render_w: int = field(init=False, default=640)
    render_h: int = field(init=False, default=360)

    def resolve(self) -> "SimConfig":
        """Resolve preset values and validate the combination of arguments."""
        if self.quality not in QUALITY:
            raise ValueError(f"unknown quality {self.quality!r}")
        preset = QUALITY[self.quality]
        n = int(self.particles if self.particles is not None
                else preset["particles"])  # type: ignore[arg-type]
        if n < 1000:
            raise ValueError("--particles must be at least 1000")
        if n > MAX_PARTICLES:
            raise ValueError(f"--particles must not exceed {MAX_PARTICLES:,}")
        g = int(self.grid if self.grid is not None else preset["grid"])  # type: ignore[arg-type]
        if g < 8 or g > MAX_GRID:
            raise ValueError(f"--grid must be between 8 and {MAX_GRID}")
        if g & (g - 1):
            # A power of two is not strictly required, but keeps FFT sizes fast.
            pass
        rw, rh = preset["render"]  # type: ignore[misc]
        ww, wh = preset["window"]  # type: ignore[misc]
        # The pure-NumPy fallback is intentionally conservative: it is slow
        # and allocates much more temporary memory than the JIT kernels.
        if importlib.util.find_spec("numba") is None:
            warnings.warn(
                "numba is not installed; using the pure-NumPy fallback and "
                "reducing the default particle count and grid to 50,000 / 32.",
                RuntimeWarning, stacklevel=2)
            n = min(n, 50_000)
            g = min(g, 32)
            rw = min(int(rw), 480)
            rh = min(int(rh), 270)
        self.n_particles = n
        self.n_grid = g
        self.render_w = int(rw)
        self.render_h = int(rh)
        self.win_w = int(self.width if self.width is not None else ww)
        self.win_h = int(self.height if self.height is not None else wh)
        if self.win_w < 320 or self.win_h < 200:
            raise ValueError("window must be at least 320x200")
        if self.scene not in ("collision", "cosmic"):
            raise ValueError("scene must be 'collision' or 'cosmic'")
        if self.scene == "collision":
            if self.mass_ratio <= 0:
                raise ValueError("--mass-ratio must be positive")
            if self.pericenter <= 0 or self.separation <= self.pericenter:
                raise ValueError("--separation must exceed --pericenter > 0")
            if self.eccentricity < 0:
                raise ValueError("--eccentricity must be non-negative")
        else:
            if self.z_start < 0:
                raise ValueError("--z-start must be non-negative")
        return self

    def default_steps(self) -> int:
        return DEFAULT_STEPS[self.scene]

    def summary(self) -> str:
        return (f"scene={self.scene} quality={self.quality} "
                f"N={self.n_particles:,} grid={self.n_grid}³ "
                f"window={self.win_w}x{self.win_h}")


# ---------------------------------------------------------------------------
# Memory guard
# ---------------------------------------------------------------------------
def available_memory_bytes() -> int:
    """Return available physical memory in bytes (Linux sysconf fallback)."""
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        return int(os_sysconf_pages() * 4096)
    except Exception:  # pragma: no cover - extremely defensive
        return 4 * 1024 ** 3


def os_sysconf_pages() -> int:
    import os
    return os.sysconf("SC_AVPHYS_PAGES")


def estimate_memory_bytes(cfg: SimConfig) -> int:
    """Conservative peak-memory estimate for particles, grids and FFTs."""
    n = cfg.n_particles
    g = cfg.n_grid
    particle = n * (3 * 4 * 2 + 4 * 5)          # pos, vel, mass, ids, speed
    snaps = n * 3 * 4 * 2                        # double-buffered snapshots
    if cfg.scene == "collision":
        nfft = 2 * g
        fft = nfft ** 3 * 4 * 3                  # density, potential, kernel
        deposit = min(1_200_000_000, nfft ** 3 * 4 * 8)
        grid = nfft ** 3 * 4
    else:
        fft = g ** 3 * 4 * 3
        deposit = min(1_200_000_000, g ** 3 * 4 * 8)
        grid = g ** 3 * 4
    render = 1_200_000_000                        # 4K thread-local RGB buffers
    return int(particle + snaps + fft + deposit + grid + render + 300_000_000)


def check_memory(cfg: SimConfig) -> None:
    """Raise a clear MemoryError when the requested run does not fit."""
    need = estimate_memory_bytes(cfg)
    have = available_memory_bytes()
    if need > have * 0.9:
        raise MemoryError(
            f"not enough memory for {cfg.n_particles:,} particles on a "
            f"{cfg.n_grid}³ grid: need about {need / 1e9:.1f} GB, "
            f"available {have / 1e9:.1f} GB. Reduce --particles/--grid or "
            f"use a lower --quality preset.")


# ---------------------------------------------------------------------------
# Palette definitions used by the software renderer
# ---------------------------------------------------------------------------
GALAXY_COLORS: np.ndarray = np.array(
    [[0.15, 0.75, 1.00], [1.00, 0.50, 0.12]], dtype=np.float32)
COMPONENT_COLORS: np.ndarray = np.array(
    [[0.22, 0.22, 0.30], [1.00, 0.62, 0.22], [0.62, 0.80, 1.00]],
    dtype=np.float32)

# Inferno-like control points (t, r, g, b).
INFERNO_ANCHORS = np.array([
    [0.00, 0.001, 0.000, 0.014],
    [0.13, 0.120, 0.045, 0.320],
    [0.28, 0.330, 0.070, 0.430],
    [0.45, 0.560, 0.130, 0.360],
    [0.62, 0.780, 0.280, 0.180],
    [0.78, 0.930, 0.540, 0.060],
    [0.90, 0.990, 0.780, 0.180],
    [1.00, 0.990, 0.980, 0.750],
], dtype=np.float32)

# Cosmic density colormap: black -> indigo -> violet -> magenta -> gold -> white.
WEB_ANCHORS = np.array([
    [0.00, 0.000, 0.000, 0.000],
    [0.18, 0.020, 0.015, 0.090],
    [0.36, 0.080, 0.040, 0.260],
    [0.54, 0.240, 0.090, 0.470],
    [0.70, 0.530, 0.220, 0.620],
    [0.84, 0.880, 0.520, 0.400],
    [0.93, 1.000, 0.780, 0.320],
    [1.00, 1.000, 0.980, 0.900],
], dtype=np.float32)

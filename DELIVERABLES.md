# Source listing
## `README.md`

```markdown
[IMAGE PROMPT: A cinematic three-quarter view of two spiral galaxies colliding in deep space, the larger blue-white disk edged with pink star-forming knots, the smaller orange-gold disk stretched into a tidal tail, a faint violet cosmic-web backdrop, volumetric glow, ultra-detailed astrophotography, 16:9]

# Galactic N-Body

A compact, fully self-contained 3D N-body laboratory: two colliding disk galaxies
and a ΛCDM cosmic web, both evolved with the same particle-mesh gravity solver,
rendered with a software splatterer, and driven by a pygame interface or from the
command line. It supports up to **8,000,000 particles** and runs end-to-end on a
CPU-only machine.

No OpenGL, no GPU, no audio, no network. Python 3.10+ with numpy, scipy, numba,
pygame, pillow and matplotlib.

---

## Screenshots

### 1. Galaxy collision

> **Screenshot placeholder 1 of 3** — save a live capture to
> `images/screenshot_collision.png`.

![Galaxy collision](images/screenshot_collision.png)

### 2. Cosmic web

> **Screenshot placeholder 2 of 3** — save a live capture to
> `images/screenshot_cosmic.png`.

![Cosmic web](images/screenshot_cosmic.png)

### 3. Analysis output

> **Screenshot placeholder 3 of 3** — run `python main.py --analyze` and keep
> `images/analysis_cosmic.png` (or `images/analysis_collision.png`).

![Analysis plots](images/screenshot_analysis.png)

---

## Features

- Two scenes: a **galaxy collision** with a Hernquist halo + bulge and an
  exponential/sech² disk, and a **cosmic web** grown from Zel'dovich initial
  conditions in a flat ΛCDM box.
- Particle-mesh gravity with cloud-in-cell (CIC) deposition and interpolation,
  Fourier Poisson solve, CIC window deconvolution, Gaussian smoothing and a
  4-point gradient. The collision scene uses an isolated **Hockney–Eastwood 2×
  padded** solve with a `-G/r` Green's function; the cosmic scene uses a
  periodic solve with the `k = 0` mode removed.
- Kick–drift–kick leapfrog with a **global adaptive timestep**.
- Self-consistent galaxy initial conditions: Jeans dispersions, Toomre
  `Q = 1.5`, asymmetric drift, a Kepler pair orbit, and two independently
  inclined disks.
- Adaptive isolated box with 10% hysteresis, refreshed every 50 steps.
- Live pygame UI: double-buffered, lock-protected snapshots, orbit camera,
  colour modes, dark-matter toggle, cosmic slab control, PNG export, 1920×1080
  capture, frame recording, restart and a help overlay.
- Headless modes: single PNG, PNG sequence, and analysis plots. **Pygame is
  never imported in headless mode.**
- Pure-NumPy fallback when numba is unavailable: it warns and automatically
  reduces the defaults.

---

## Installation

```bash
git clone <your-repo-url> galactic-nbody
cd galactic-nbody
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Python 3.10 or newer is required. On a fresh machine the first numba call
compiles the kernels and takes a few seconds; later runs use the on-disk cache.

---

## Quick Start

```bash
# Interactive galaxy collision (250k particles, 64^3 grid)
python main.py

# Interactive cosmic web
python main.py --scene cosmic

# Headless single renders
python main.py --render images/collision.png
python main.py --scene cosmic --render images/cosmic.png

# A PNG sequence, then analysis
python main.py --frames images/frames --frame-count 90
python main.py --scene cosmic --analyze
```

---

## Command-Line Options

| Option | Meaning | Default |
| --- | --- | --- |
| `--scene {collision,cosmic}` | Simulation scene | `collision` |
| `--quality {preview,high,ultra}` | Particle/grid/resolution preset | `preview` |
| `--particles N` | Particle count, `1000 … 8,000,000` | preset |
| `--grid G` | PM grid size per dimension, `8 … 512` | preset |
| `--width W`, `--height H` | Interactive window size | preset |
| `--seed S` | Random seed | `12345` |
| `--box L` | Cosmic box side in Mpc/h | `100` |
| `--z-start Z` | Cosmic starting redshift | `49` |
| `--mass-ratio q` | Secondary / primary mass | `1.0` |
| `--pericenter R` | Pair pericentre in kpc | `15` |
| `--eccentricity e` | Pair eccentricity | `1.0` |
| `--separation R` | Initial pair separation in kpc | `80` |
| `--incl1 DEG`, `--incl2 DEG` | Disk inclinations | `0`, `60` |
| `--steps N` | Fixed number of steps (headless) | scene endpoint |
| `--render PATH` | Write one PNG and exit | — |
| `--frames DIR` | Write a PNG sequence and exit | — |
| `--frame-count N` | Number of frames | `90` |
| `--analyze` | Run and save a diagnostic plot | — |

Without `--steps`, the collision scene advances to about **0.55 Gyr** (safety cap
800 steps) and the cosmic scene stops exactly at **a = 1** (safety cap 400
steps).

### Quality presets

| Preset | Particles | Grid | Render size | Typical use |
| --- | --- | --- | --- | --- |
| `preview` | 250,000 | 64³ | 640×360 | live interaction |
| `high` | 2,000,000 | 128³ | 960×540 | live / short movies |
| `ultra` | 8,000,000 | 256³ | 1280×720 | offline renders |

`ultra` is deliberately an offline target: the UI remains responsive because
simulation and rendering run in separate threads, but a live 8M-particle view is
limited by the software renderer.

---

## Interactive Controls

| Input | Action |
| --- | --- |
| Left mouse drag | Orbit the camera |
| Mouse wheel | Zoom |
| Arrow keys | Orbit / change inclination |
| `SPACE` | Pause / resume |
| `C` | Cycle colour mode (galaxy, component, density, speed) |
| `D` | Show / hide dark matter (collision) |
| `+` / `-` | Increase / decrease simulation speed |
| `[` / `]` | Decrease / increase cosmic slab thickness |
| `S` | Save the current preview PNG to `images/` |
| `F` | Save a 1920×1080 PNG |
| `V` | Start / stop frame recording |
| `R` | Reset the camera |
| `T` | Restart with a new seed |
| `H` | Toggle the help overlay |
| `ESC` | Quit |

---

## Physics

### Units

The collision scene uses **kpc, km/s, M☉**, so `G = 4.30091 × 10⁻⁶
kpc (km/s)² M☉⁻¹` and one time unit is `0.9778 Gyr`. The cosmic scene uses
comoving **Mpc/h** and `H₀ = 100 km/s/(Mpc/h)`, i.e. `h = 0.68`,
`Ω_m = 0.30`, `Ω_Λ = 0.70`, `σ₈ = 0.80`, `n_s = 0.96`.

### Particle-mesh gravity

Particles are assigned to the grid with CIC weights, the density is transformed
with `scipy.fft.rfftn(..., workers=-1)` (or `numpy.fft` if scipy is absent),
multiplied by the Poisson kernel, and differentiated with a 4-point central
gradient. Each Fourier mode is divided by the CIC assignment window so that the
force is not artificially smoothed by the mass-assignment kernel. A Gaussian
smoothing of `0.30` cells suppresses grid-scale noise.

- **Periodic (cosmic):** `k = 0` is set to zero; the kernel uses `1/h³` so that
  particle masses are dimensionless overdensity weights.
- **Isolated (collision):** Hockney–Eastwood 2× zero padding and the analytic
  `-G/r` Green's function give vacuum boundary conditions. Physical cell masses
  are used directly.

The CIC deposit conserves mass exactly, and the periodic force conserves total
momentum to round-off.

### Galaxy initial conditions

The primary has `M_halo = 10¹² M☉` (Hernquist, `a = 30 kpc`, truncated at
100 kpc), `M_bulge = 10¹⁰ M☉` (Hernquist, `a = 0.7 kpc`), and
`M_disk = 5 × 10¹⁰ M☉` (exponential, `R_d = 3.5 kpc`, `z_d = 0.35 kpc`,
sech² vertical profile). The secondary is scaled by `--mass-ratio`. The stellar
disk always receives at least 25% of the particles; the rest follow the mass
fractions.

Radial velocity dispersions come from the spherical Jeans equation,

```
σ²(r) = (1/ρ) ∫_r^R ρ(r') G M(<r') / r'² dr',
```

and the disk is given a Toomre `Q = 1.5` radial dispersion. Asymmetric drift is
included through

```
v_φ² = v_c² − σ_R² [ d ln Σ/d ln R + d ln σ_R²/d ln R
                     + 1 − κ²/(4Ω²) ].
```

The pair is placed on a Kepler orbit from `--separation`, `--pericenter` and
`--eccentricity`, with independent disk inclinations. The two galaxies are then
recentred to zero total momentum.

### Cosmic initial conditions

The Eisenstein–Hu no-wiggle transfer function is normalized to `σ₈ = 0.8`, and
a deterministic Gaussian random field is generated in Fourier space from the
same seed. Particles start at `z = 49` on a uniform grid with Zel'dovich
displacements and velocities. The linear growth factor is obtained by exact
integration of

```
D(a) ∝ H(a) ∫_0^a da' / (a' H(a'))³,
```

normalized to `D(1) = 1`. For an Einstein–de Sitter universe this reduces to
`D = a`; for the default cosmology `D(1) = 0.779` before normalization.

The cosmic evolution is a comoving leapfrog in `ln a`:

```
dx/d ln a = p / (a² H),
dp/d ln a = −∇φ / H,
∇²φ = (3/2) Ω_m δ / a.
```

### Integrator

Both scenes use kick–drift–kick leapfrog. The collision timestep is

```
dt = min( η √(ε/a_max), 0.5 · cell / v_max ),
```

with `η = 0.2`. The cosmic timestep is limited by the maximum displacement per
`ln a` and by the remaining distance to `a = 1`. The collision box is refreshed
every 50 steps only when the 99th-percentile radius moved by more than 10%
(10% hysteresis).

---

## Output

- `--render PATH` writes one PNG after advancing the scene.
- `--frames DIR` writes `frame_0000.png`, `frame_0001.png`, … evenly spaced
  over the requested run.
- `--analyze` writes `images/analysis_cosmic.png` (measured versus linear
  `P(k)`) or `images/analysis_collision.png` (dark-matter density and enclosed
  mass versus Hernquist).
- The live UI writes captures to `images/render_*.png`, `images/frame_1920_*.png`
  and, while recording, `images/record_*/frame_*.png`.

---

## Testing

```bash
python -m pytest -q
```

The suite checks:

- CIC mass conservation in periodic and isolated boxes;
- the isolated point-mass force against `G M / r²` beyond 3–4 cells (5%);
- periodic momentum conservation;
- leapfrog energy drift over 10 orbits (< 10⁻³) and time reversibility;
- Hernquist enclosed mass and circular velocity against numerical integrals
  (5%), plus a sampled-halo virial ratio in `0.7 … 1.3`;
- `D(1) ≈ 0.779` unnormalised, `D(1) = 1` normalised, and `D = a` for EdS;
- a cosmic age of 13.5–13.9 Gyr and the `σ₈` normalization.

With warm numba caches the suite finishes in about one second; a cold run is
still well under a minute.

---

## Performance Notes

- `--quality preview` targets interactive use on 8 CPU cores; `high` is a good
  compromise for short captures; `ultra` is intended for `--render`/`--frames`.
- The renderer accumulates into one RGB buffer per thread and merges the
  buffers afterwards, so there are no races and no global locks in the splat
  loop.
- The simulation runs in a background thread and publishes immutable
  double-buffered snapshots. The UI never reads a particle array while it is
  being mutated.
- The memory guard estimates particle arrays, double buffers, FFT workspace,
  deposition temporaries and render buffers. If the request exceeds 90% of
  available RAM, the program exits with a clear error before allocating.
- On a machine without numba the program warns, forces the pure-NumPy path and
  caps the defaults at 50,000 particles and a 32³ grid.

---

## Limitations

This is an educational, collisionless PM code, not a research cosmology solver.
Be aware of the following:

- **No gas, star formation, feedback, black holes or magnetic fields.** Both
  scenes are pure N-body.
- **PM resolution is about two cells.** Structures smaller than that are
  smoothed away. In the collision scene the dark-matter cusp inside roughly two
  cells is flattened by the grid; the analysis plot therefore starts at two
  cells. In the cosmic scene a cell is 1–2 Mpc/h at the default settings, so
  individual halos and galaxies are not resolved.
- **The collision profile is not converged below the softening scale.** The
  solver uses cell-scale Gaussian smoothing rather than per-particle softening;
  very close encounters are softened by the grid.
- **Cosmic energy is not conserved.** The comoving equations include explicit
  time dependence through `a` and `H(a)`; the displayed cosmic energy is a
  diagnostic, not a conserved quantity. The UI hides its drift for that reason.
- **Global timestep.** All particles share one `dt`, chosen from the global
  maximum acceleration and speed. A single tightly bound particle can slow the
  whole run.
- **No MPI, GPU or tree gravity.** Everything is single-node, CPU-only, and the
  FFT cost scales as `N_grid³ log N_grid`. An 8M-particle `ultra` run with a
  256³ grid is an offline job.
- **The renderer is an artistic splatterer, not radiative transfer.** Colour
  modes, exposure, bloom and the cosmic slab are visual choices; they are not
  calibrated observables.
- **Initial conditions are approximate.** The Jeans/Toomre setup produces a
  reasonably cool disk, but it is not a full distribution-function model;
  transients at startup are possible, especially at low particle counts.
- **Pygame is required only for the live UI.** Headless rendering and tests
  never create a window.

---

## References

1. Hockney, R. W. & Eastwood, J. W. (1988). *Computer Simulation Using
   Particles*. Adam Hilger.
2. Hernquist, L. (1990). An analytical model for spherical galaxies and
   bulges. *ApJ*, 356, 359.
3. Toomre, A. (1964). On the gravitational stability of a disk of stars.
   *ApJ*, 139, 1217.
4. Binney, J. & Tremaine, S. (2008). *Galactic Dynamics* (2nd ed.). Princeton
   University Press.
5. Springel, V. (2005). The cosmological simulation code GADGET-2. *MNRAS*,
   364, 1105.
6. Eisenstein, D. J. & Hu, W. (1998). Baryonic features in the matter transfer
   function. *ApJ*, 496, 605.
7. Zel'dovich, Ya. B. (1970). Gravitational instability: an approximate theory
   for large density perturbations. *A&A*, 5, 84.
8. Barnes, J. & Hut, P. (1986). A hierarchical O(N log N) force-calculation
   algorithm. *Nature*, 324, 446.
9. Press, W. H., Teukolsky, S. A., Vetterling, W. T. & Flannery, B. P. (2007).
   *Numerical Recipes* (3rd ed.). Cambridge University Press.
10. Planck Collaboration (2020). Planck 2018 results. VI. Cosmological
    parameters. *A&A*, 641, A6.

---

## License

MIT. See [LICENSE](LICENSE).

```
## `requirements.txt`

```text
numpy>=1.24
scipy>=1.10
numba>=0.57
pygame>=2.5
pillow>=10.0
matplotlib>=3.7
pytest>=7.4

```
## `LICENSE`

```text
MIT License

Copyright (c) 2026 qwfrx

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

```
## `.gitignore`

```text
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
dist/
*.egg-info/
.eggs/

# Virtual environments
.venv/
venv/
env/

# Tooling
.pytest_cache/
.mypy_cache/
.ruff_cache/
.coverage
htmlcov/

# Generated renders and analysis output
images/*.png
images/*.jpg
images/recording/
images/frames/
!images/.gitkeep

# Numba cache
__pycache__/

# OS/editor
.DS_Store
Thumbs.db
.vscode/
.idea/

```
## `main.py`

```python
"""Command-line entry point for the galactic N-body simulation.

Examples
--------
    python main.py
    python main.py --scene cosmic --quality preview --steps 120
    python main.py --render images/collision.png --steps 220
    python main.py --frames images/frames --frame-count 90 --steps 200
    python main.py --analyze --scene cosmic --steps 200

Pygame is imported only in the interactive branch, so every headless mode runs
on machines without a display server.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Any

import numpy as np

import config
from src.analysis import analyze
from src.camera import Camera
from src.integrator import CollisionSim, CosmicSim, make_simulation, run_steps
from src.render import Renderer, save_png


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    p = argparse.ArgumentParser(
        prog="main.py",
        description="3D N-body galaxy collision and cosmic-web simulation.")
    p.add_argument("--scene", choices=("collision", "cosmic"),
                   default="collision", help="simulation scene")
    p.add_argument("--quality", choices=tuple(config.QUALITY),
                   default="preview", help="quality preset")
    p.add_argument("--particles", type=int, default=None,
                   help="particle count (up to 8,000,000)")
    p.add_argument("--grid", type=int, default=None,
                   help="PM grid size per dimension")
    p.add_argument("--width", type=int, default=None,
                   help="interactive window width")
    p.add_argument("--height", type=int, default=None,
                   help="interactive window height")
    p.add_argument("--seed", type=int, default=12345,
                   help="random seed for the initial conditions")
    p.add_argument("--box", type=float, default=None,
                   help="cosmic box side in Mpc/h")
    p.add_argument("--z-start", type=float, default=config.Z_START,
                   help="starting redshift of the cosmic scene")
    p.add_argument("--mass-ratio", type=float, default=1.0,
                   help="secondary / primary mass ratio")
    p.add_argument("--pericenter", type=float, default=15.0,
                   help="pericentre of the Kepler pair orbit in kpc")
    p.add_argument("--eccentricity", type=float, default=1.0,
                   help="eccentricity of the pair orbit")
    p.add_argument("--separation", type=float, default=80.0,
                   help="initial separation in kpc")
    p.add_argument("--incl1", type=float, default=0.0,
                   help="primary disk inclination in degrees")
    p.add_argument("--incl2", type=float, default=60.0,
                   help="secondary disk inclination in degrees")
    p.add_argument("--steps", type=int, default=None,
                   help="number of steps for headless runs")
    p.add_argument("--render", metavar="PATH", default=None,
                   help="write a single PNG and exit")
    p.add_argument("--frames", metavar="DIR", default=None,
                   help="write a PNG sequence to DIR and exit")
    p.add_argument("--frame-count", type=int, default=90,
                   help="number of frames for --frames (default 90)")
    p.add_argument("--analyze", action="store_true",
                   help="run the scene and save a diagnostic plot")
    return p


def make_config(args: argparse.Namespace) -> config.SimConfig:
    """Turn parsed arguments into a validated :class:`config.SimConfig`."""
    cfg = config.SimConfig(
        scene=args.scene,
        quality=args.quality,
        particles=args.particles,
        grid=args.grid,
        width=args.width,
        height=args.height,
        seed=args.seed,
        box=args.box,
        z_start=args.z_start,
        mass_ratio=args.mass_ratio,
        pericenter=args.pericenter,
        eccentricity=args.eccentricity,
        separation=args.separation,
        incl1=args.incl1,
        incl2=args.incl2,
        steps=args.steps,
        render=args.render,
        frames=args.frames,
        frame_count=args.frame_count,
        analyze=args.analyze,
    )
    if cfg.steps is not None and cfg.steps < 0:
        raise ValueError("--steps must be non-negative")
    return cfg.resolve()


def advance(sim: CollisionSim | CosmicSim, cfg: config.SimConfig) -> int:
    """Run the requested scene to its default or explicit endpoint."""
    if cfg.steps is not None:
        run_steps(sim, cfg.steps, progress=True)
        return int(cfg.steps)
    cap = cfg.default_steps()
    count = 0
    if isinstance(sim, CollisionSim):
        while sim.stats["time_gyr"] < 0.55 and count < cap:
            sim.step()
            count += 1
            if count % max(1, cap // 10) == 0:
                print(f"  step {count}/{cap}", flush=True)
    else:
        while sim.a < 1.0 - 1.0e-10 and count < cap:
            sim.step()
            count += 1
            if count % max(1, cap // 10) == 0:
                print(f"  step {count}/{cap} (a={sim.a:.3f})", flush=True)
    return count


def frame_camera(snap: dict[str, Any]) -> Camera:
    """Choose a reproducible camera for a stand-alone render or movie frame."""
    camera = Camera()
    box = float(snap["box"])
    extent = (0.82 if snap["scene"] == "cosmic" else 0.28) * box
    camera.reset(np.asarray(snap["center"], dtype=np.float64), extent)
    return camera


def render_still(sim: CollisionSim | CosmicSim, cfg: config.SimConfig,
                 path: str) -> str:
    """Render the current state at the preset resolution and save it."""
    snap = sim.snapshot()
    camera = frame_camera(snap)
    renderer = Renderer(cfg.render_w, cfg.render_h)
    frame = renderer.render(snap, camera, color_mode=0, show_dm=True,
                            exposure=5.0, slab=0.35)
    save_png(frame, path)
    return path


def export_frames(sim: CollisionSim | CosmicSim, cfg: config.SimConfig,
                  directory: str, count: int) -> int:
    """Render ``count`` frames spread over the requested run and save them."""
    if count < 1:
        raise ValueError("--frame-count must be at least 1")
    os.makedirs(directory, exist_ok=True)
    renderer = Renderer(cfg.render_w, cfg.render_h)
    written = 0

    if cfg.steps is not None:
        total = max(1, int(cfg.steps))
        for step in range(1, total + 1):
            sim.step()
            target = max(1, round(total * (written + 1) / count))
            if step >= target:
                snap = sim.snapshot()
                camera = frame_camera(snap)
                frame = renderer.render(snap, camera, 0, True, 5.0, 0.35)
                path = os.path.join(directory, f"frame_{written:04d}.png")
                save_png(frame, path)
                written += 1
                if written % 10 == 0:
                    print(f"  frame {written}/{count}", flush=True)
                if written >= count:
                    break
    elif isinstance(sim, CollisionSim):
        targets = np.linspace(0.55 / count, 0.55, count)
        next_target = 0
        while next_target < count and sim.steps < cfg.default_steps():
            sim.step()
            while (next_target < count
                   and sim.stats["time_gyr"] >= targets[next_target]):
                snap = sim.snapshot()
                camera = frame_camera(snap)
                frame = renderer.render(snap, camera, 0, True, 5.0, 0.35)
                path = os.path.join(directory, f"frame_{written:04d}.png")
                save_png(frame, path)
                written += 1
                next_target += 1
        # Duplicate the last frame if the target time was reached early.
        while written < count:
            snap = sim.snapshot()
            camera = frame_camera(snap)
            frame = renderer.render(snap, camera, 0, True, 5.0, 0.35)
            save_png(frame, os.path.join(directory, f"frame_{written:04d}.png"))
            written += 1
    else:
        a0 = max(sim.a, 1.0e-6)
        targets = a0 * (1.0 / a0) ** (np.arange(1, count + 1) / count)
        next_target = 0
        while next_target < count and sim.a < 1.0 - 1.0e-10:
            sim.step()
            while (next_target < count
                   and sim.a >= targets[next_target] - 1.0e-12):
                snap = sim.snapshot()
                camera = frame_camera(snap)
                frame = renderer.render(snap, camera, 0, True, 5.0, 0.35)
                path = os.path.join(directory, f"frame_{written:04d}.png")
                save_png(frame, path)
                written += 1
                next_target += 1
        while written < count:
            snap = sim.snapshot()
            camera = frame_camera(snap)
            frame = renderer.render(snap, camera, 0, True, 5.0, 0.35)
            save_png(frame, os.path.join(directory, f"frame_{written:04d}.png"))
            written += 1
    return written


def print_summary(sim: CollisionSim | CosmicSim, cfg: config.SimConfig) -> None:
    """Print a compact run report."""
    snap = sim.snapshot()
    st = snap["stats"]
    print(f"  {cfg.summary()}")
    if snap["scene"] == "collision":
        print(f"  t = {st['time_gyr']:.3f} Gyr, separation = "
              f"{st['separation_kpc']:.1f} kpc, virial = {st['virial']:.3f}, "
              f"energy drift = {100.0 * st['energy_drift']:+.3f}%")
    else:
        print(f"  a = {st['a']:.4f}, z = {st['z']:.3f}, "
              f"age = {st['age_gyr']:.3f} Gyr, delta_rms = "
              f"{st['delta_rms']:.3f}, virial = {st['virial']:.3f}")


def main(argv: list[str] | None = None) -> int:
    """Program entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg = make_config(args)
        config.check_memory(cfg)
    except (ValueError, MemoryError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"Galactic N-body: {cfg.summary()}")
    start = time.time()
    headless = bool(args.render or args.frames or args.analyze)
    if not headless:
        from src.ui import App
        App(cfg).run()
        return 0

    sim = make_simulation(cfg)
    if args.frames:
        count = export_frames(sim, cfg, args.frames, args.frame_count)
        print(f"  wrote {count} frames to {args.frames}")
        print_summary(sim, cfg)
    else:
        steps = advance(sim, cfg)
        print(f"  advanced {steps} steps in {time.time() - start:.1f} s")
        print_summary(sim, cfg)

    if args.render:
        path = render_still(sim, cfg, args.render)
        print(f"  saved {path}")

    if args.analyze:
        analyze(cfg, sim)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

```
## `config.py`

```python
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

```
## `src/__init__.py`

```python

```
## `src/units.py`

```python
"""Unit conversions and cosmological background functions.

Collision units are kpc, km/s and Msun with G = 4.30091e-6.  Cosmic units
are Mpc/h, 100 km/s and code time (14.379 Gyr for h = 0.68), with H0 = 1.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

import config


def collision_time_gyr(t_code: np.ndarray | float) -> np.ndarray | float:
    """Convert collision code time to Gyr."""
    return np.asarray(t_code) * config.KPC_KMPS_GYR


def dynamical_time_gyr(density_msun_kpc3: np.ndarray | float) -> np.ndarray | float:
    """Spherical dynamical time t_dyn = sqrt(3 pi / (32 G rho))."""
    rho = np.maximum(np.asarray(density_msun_kpc3, dtype=np.float64), 1e-30)
    return np.sqrt(3.0 * np.pi / (32.0 * config.G_KPC * rho)) * config.KPC_KMPS_GYR


def hubble_e(a: np.ndarray | float, omega_m: float = config.OMEGA_M,
              omega_l: float = config.OMEGA_L) -> np.ndarray | float:
    """Dimensionless Friedmann expansion rate H(a)/H0.

    H(a) = H0 sqrt(Omega_m a^-3 + Omega_L) for a flat universe.
    """
    a_arr = np.asarray(a, dtype=np.float64)
    e = np.sqrt(omega_m * a_arr ** -3 + omega_l)
    return e if e.ndim else float(e)


def cosmic_age_gyr(a: np.ndarray | float, omega_m: float = config.OMEGA_M,
                   omega_l: float = config.OMEGA_L,
                   h: float = config.H_DIMENSIONLESS) -> np.ndarray | float:
    """Analytic age of a flat LambdaCDM universe.

    t(a) = 2 / (3 H0 sqrt(Omega_L)) * asinh(sqrt(Omega_L/Omega_m) a^(3/2)).
    The result is returned in Gyr; H0 is 100 h km/s/Mpc.
    """
    a_arr = np.asarray(a, dtype=np.float64)
    if omega_l <= 0.0:
        factor = 2.0 / (3.0 * np.sqrt(omega_m)) * a_arr ** 1.5
    else:
        factor = 2.0 / (3.0 * np.sqrt(omega_l)) * np.arcsinh(
            np.sqrt(omega_l / omega_m) * a_arr ** 1.5)
    unit = 977.792221 / (100.0 * h)  # Gyr per unit of 1/H0
    out = factor * unit
    return out if out.ndim else float(out)


def _growth_integrand(x: float, omega_m: float, omega_l: float) -> float:
    e = np.sqrt(omega_m * x ** -3 + omega_l)
    return 1.0 / (x * e) ** 3


def growth_factor(a: np.ndarray | float, omega_m: float = config.OMEGA_M,
                  omega_l: float = config.OMEGA_L,
                  normalized: bool = True) -> np.ndarray | float:
    """Linear growth factor from the exact integral solution.

    D(a) = (5 Omega_m / 2) H(a)/H0 * integral_0^a da' / (a' H(a')/H0)^3.

    With ``normalized=True`` D(1) = 1.  With ``normalized=False`` the result
    tends to a as a -> 0, so D(1) ~ 0.78 for the default cosmology and
    D(a) = a exactly for Einstein-de Sitter.  The integral is evaluated with
    adaptive quadrature and is deterministic.
    """
    def one(x: float) -> float:
        if x <= 0.0:
            return 0.0
        e = np.sqrt(omega_m * x ** -3 + omega_l)
        integral, _ = quad(_growth_integrand, 0.0, x,
                           args=(omega_m, omega_l), limit=200,
                           epsabs=1e-12, epsrel=1e-10)
        return 2.5 * omega_m * e * integral

    a_arr = np.asarray(a, dtype=np.float64)
    raw = np.vectorize(one, otypes=[np.float64])(a_arr)
    if not normalized:
        return float(raw) if np.ndim(raw) == 0 else raw
    norm = one(1.0)
    out = raw / norm
    return float(out) if np.ndim(out) == 0 else out


def growth_rate(a: float, omega_m: float = config.OMEGA_M,
                omega_l: float = config.OMEGA_L) -> float:
    """Logarithmic growth rate f = d ln D / d ln a (finite difference)."""
    h = 1.0e-3
    d1 = float(growth_factor(a * np.exp(-h), omega_m, omega_l))
    d2 = float(growth_factor(a * np.exp(h), omega_m, omega_l))
    return (np.log(d2) - np.log(d1)) / (2.0 * h)


def mean_matter_density(a: float, omega_m: float = config.OMEGA_M,
                        h: float = config.H_DIMENSIONLESS) -> float:
    """Mean matter density in Msun per (Mpc/h)^3 (for the info panel)."""
    rho_c = 2.77536627e11 * h ** 2  # Msun / Mpc^3
    return omega_m * rho_c / a ** 3 * h ** 3  # per (Mpc/h)^3


def sky_angle(x: float, y: float) -> float:
    """Helper used by the info panel for a rough angular scale."""
    return float(np.degrees(np.arctan2(y, x)))

```
## `src/ics_galaxy.py`

```python
"""Galaxy models and initial conditions for the collision scene.

The primary galaxy is Milky-Way-like.  The secondary is obtained by scaling
all component masses by the mass ratio q.  Particle counts follow component
mass fractions, except that the stellar disk always receives at least
MIN_DISK_SHARE of the particles so that tidal tails are well sampled.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import config
from . import units


@dataclass
class Galaxy:
    pos: np.ndarray      # (N, 3) float32, kpc, centred on zero
    vel: np.ndarray      # (N, 3) float32, km/s, zero net momentum
    mass: np.ndarray     # (N,) float32, Msun
    comp: np.ndarray     # (N,) int8: 0 halo, 1 bulge, 2 disk
    gal: np.ndarray      # (N,) int8: 0 or 1
    mass_total: float


@dataclass
class CollisionIC:
    pos: np.ndarray
    vel: np.ndarray
    mass: np.ndarray
    comp: np.ndarray
    gal: np.ndarray
    box: float
    center: np.ndarray
    mass_total: float
    mass_primary: float


# ---------------------------------------------------------------------------
# Analytical profiles
# ---------------------------------------------------------------------------
def hernquist_mass(r: np.ndarray | float, mass: float, scale: float
                   ) -> np.ndarray | float:
    """Hernquist enclosed mass M(<r) = M r^2 / (r+a)^2."""
    rr = np.asarray(r, dtype=np.float64)
    out = mass * rr ** 2 / (rr + scale) ** 2
    return out if out.ndim else float(out)


def hernquist_density(r: np.ndarray | float, mass: float, scale: float
                      ) -> np.ndarray | float:
    """Hernquist density rho = M a / (2 pi r (r+a)^3)."""
    rr = np.maximum(np.asarray(r, dtype=np.float64), 1.0e-8)
    out = mass * scale / (2.0 * np.pi * rr * (rr + scale) ** 3)
    return out if out.ndim else float(out)


def hernquist_potential(r: np.ndarray | float, mass: float, scale: float
                        ) -> np.ndarray | float:
    """Hernquist potential Phi = -G M / (r+a)."""
    rr = np.asarray(r, dtype=np.float64)
    out = -config.G_KPC * mass / (rr + scale)
    return out if out.ndim else float(out)


def miyamoto_nagai_potential(radius: np.ndarray, z: np.ndarray, mass: float,
                             scale: float, thickness: float) -> np.ndarray:
    """Miyamoto-Nagai disk potential Phi = -G M / sqrt(R^2+(a+sqrt(z^2+b^2))^2)."""
    return -config.G_KPC * mass / np.sqrt(
        radius ** 2 + (scale + np.sqrt(z ** 2 + thickness ** 2)) ** 2)


def _total_potential(radius: np.ndarray, z: np.ndarray,
                     mh: float, mb: float, md: float,
                     ah: float, ab: float, ad: float, bd: float) -> np.ndarray:
    r = np.sqrt(radius ** 2 + z ** 2)
    return (hernquist_potential(r, mh, ah)
            + hernquist_potential(r, mb, ab)
            + miyamoto_nagai_potential(radius, z, md, ad, bd))


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
def sample_hernquist(n: int, mass: float, scale: float, r_max: float,
                     rng: np.random.Generator) -> np.ndarray:
    """Draw radii from a truncated Hernquist profile and random directions."""
    if n <= 0:
        return np.zeros((0, 3), dtype=np.float64)
    frac = float(hernquist_mass(r_max, mass, scale)) / mass
    u = rng.random(n) * frac
    s = np.sqrt(np.maximum(u, 0.0))
    r = scale * s / np.maximum(1.0 - s, 1.0e-8)
    r = np.minimum(r, r_max)
    cos_t = 2.0 * rng.random(n) - 1.0
    phi = 2.0 * np.pi * rng.random(n)
    sin_t = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0))
    return np.column_stack((r * sin_t * np.cos(phi),
                            r * sin_t * np.sin(phi),
                            r * cos_t))


def sample_exponential_disk(n: int, scale: float, thickness: float,
                            r_max: float, rng: np.random.Generator
                            ) -> np.ndarray:
    """Draw an exponential disk with a sech^2 vertical profile."""
    if n <= 0:
        return np.zeros((0, 3), dtype=np.float64)
    radius = np.empty(n, dtype=np.float64)
    filled = 0
    while filled < n:
        draw = rng.gamma(2.0, scale, n - filled)
        draw = draw[draw < r_max]
        radius[filled:filled + draw.size] = draw
        filled += draw.size
    u = rng.random(n)
    z = thickness * np.arctanh(2.0 * u - 1.0)
    z = np.clip(z, -5.0 * thickness, 5.0 * thickness)
    phi = 2.0 * np.pi * rng.random(n)
    return np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))


def jeans_sigma(radii: np.ndarray, mass: float, scale: float, r_max: float,
                n_table: int = 400) -> np.ndarray:
    """Isotropic spherical Jeans dispersion for a truncated Hernquist model.

    sigma^2(r) = (1/rho) integral_r^{r_max} rho(r') G M(<r') / r'^2 dr'.
    The integral is tabulated and linearly interpolated.
    """
    r_tab = np.geomspace(max(0.005 * scale, 1.0e-3), r_max, n_table)
    rho = hernquist_density(r_tab, mass, scale)
    menc = hernquist_mass(r_tab, mass, scale)
    integrand = rho * config.G_KPC * menc / r_tab ** 2
    cumulative = np.zeros(n_table, dtype=np.float64)
    dr = np.diff(r_tab)
    for i in range(n_table - 2, -1, -1):
        cumulative[i] = cumulative[i + 1] + 0.5 * (
            integrand[i] + integrand[i + 1]) * dr[i]
    sigma2 = cumulative / np.maximum(rho, 1.0e-30)
    sigma = np.sqrt(np.maximum(sigma2, 0.0))
    sigma = np.maximum(sigma, 1.0)
    rr = np.clip(radii, r_tab[0], r_tab[-1])
    return np.interp(rr, r_tab, sigma)


def _rotation_y(angle_deg: float) -> np.ndarray:
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _rotation_x(angle_deg: float) -> np.ndarray:
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rotation_z(angle_deg: float) -> np.ndarray:
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ---------------------------------------------------------------------------
# Galaxy construction
# ---------------------------------------------------------------------------
def _circular_velocity(radius: np.ndarray, z: np.ndarray,
                       mh: float, mb: float, md: float,
                       ah: float, ab: float, ad: float, bd: float,
                       r_halo_max: float, r_bulge_max: float) -> np.ndarray:
    """v_c from the enclosed-mass model plus the MN disk derivative.

    v_c^2 = G (M_halo(<R) + M_bulge(<R)) / R + R dPhi_MN/dR.
    The spheroidal components are truncated at their sampling radii.
    """
    r_h = np.minimum(radius, r_halo_max)
    r_b = np.minimum(radius, r_bulge_max)
    m_enc = hernquist_mass(r_h, mh, ah) + hernquist_mass(r_b, mb, ab)
    vc2_sph = config.G_KPC * m_enc / np.maximum(radius, 1.0e-6)
    dr = np.maximum(1.0e-3 * radius, 2.0e-3)
    phi_p = miyamoto_nagai_potential(radius + dr, z, md, ad, bd)
    phi_m = miyamoto_nagai_potential(radius - dr, z, md, ad, bd)
    vc2_disk = np.maximum(radius * (phi_p - phi_m) / (2.0 * dr), 0.0)
    return np.sqrt(np.maximum(vc2_sph + vc2_disk, 0.0))


def build_galaxy(n: int, mass_factor: float, galaxy_id: int,
                 rng: np.random.Generator,
                 mh: float = config.M_HALO, md: float = config.M_DISK,
                 mb: float = config.M_BULGE, rd: float = config.R_DISK,
                 zd: float = config.Z_DISK, ah: float = config.A_HALO,
                 ab: float = config.A_BULGE,
                 r_halo_max: float = config.R_HALO_MAX) -> Galaxy:
    """Construct a three-component galaxy with self-consistent kinematics."""
    mh *= mass_factor
    md *= mass_factor
    mb *= mass_factor
    m_total = mh + md + mb

    disk_frac = max(config.MIN_DISK_SHARE, md / m_total)
    rest = 1.0 - disk_frac
    bulge_frac = rest * mb / (mh + mb)
    halo_frac = rest - bulge_frac

    n_disk = max(1, int(round(n * disk_frac)))
    n_bulge = max(1, int(round(n * bulge_frac)))
    n_halo = max(1, n - n_disk - n_bulge)
    if n_halo + n_bulge + n_disk != n:
        n_halo = n - n_disk - n_bulge

    # --- halo -------------------------------------------------------------
    pos_h = sample_hernquist(n_halo, mh, ah, r_halo_max, rng)
    r_h = np.linalg.norm(pos_h, axis=1)
    sig_h = jeans_sigma(r_h, mh, ah, r_halo_max)
    vel_h = rng.normal(0.0, 1.0, (n_halo, 3)) * sig_h[:, None]
    mass_h = np.full(n_halo, float(hernquist_mass(r_halo_max, mh, ah)) / n_halo,
                     dtype=np.float64)
    comp_h = np.zeros(n_halo, dtype=np.int8)

    # --- bulge ------------------------------------------------------------
    r_bulge_max = config.R_BULGE_MAX
    pos_b = sample_hernquist(n_bulge, mb, ab, r_bulge_max, rng)
    r_b = np.linalg.norm(pos_b, axis=1)
    sig_b = jeans_sigma(r_b, mb, ab, r_bulge_max)
    vel_b = rng.normal(0.0, 1.0, (n_bulge, 3)) * sig_b[:, None]
    mass_b = np.full(n_bulge, float(hernquist_mass(r_bulge_max, mb, ab)) / n_bulge,
                     dtype=np.float64)
    comp_b = np.ones(n_bulge, dtype=np.int8)

    # --- exponential disk -------------------------------------------------
    r_disk_max = config.R_DISK_MAX_FACTOR * rd
    pos_d = sample_exponential_disk(n_disk, rd, zd, r_disk_max, rng)
    radius = np.hypot(pos_d[:, 0], pos_d[:, 1])
    zpos = pos_d[:, 2]

    # Circular velocity from the total potential (halo + bulge + MN disk).
    vc = _circular_velocity(radius, zpos, mh, mb, md, ah, ab, rd, zd,
                            r_halo_max, r_bulge_max)
    vc2 = vc ** 2
    dr2 = np.maximum(1.0e-3 * radius, 2.0e-3)
    vc_p = _circular_velocity(radius + dr2, zpos, mh, mb, md, ah, ab, rd, zd,
                              r_halo_max, r_bulge_max)
    vc_m = _circular_velocity(np.maximum(radius - dr2, 1.0e-4), zpos,
                              mh, mb, md, ah, ab, rd, zd,
                              r_halo_max, r_bulge_max)
    dvc_dr = (vc_p - vc_m) / (2.0 * dr2)

    omega = vc / np.maximum(radius, 1.0e-6)
    kappa2 = np.maximum(2.0 * omega * (omega + dvc_dr), 0.0)
    kappa = np.sqrt(kappa2)

    sigma_surface = md / (2.0 * np.pi * rd ** 2) * np.exp(-radius / rd)
    sigma_r = config.TOOMRE_Q * 3.36 * config.G_KPC * sigma_surface / \
        np.maximum(kappa, 1.0e-8)

    # Asymmetric drift: v_phi^2 = v_c^2 - sigma_R^2 * (dlnSigma/dlnR
    # + dln sigma_R^2/dlnR + 1 - kappa^2/(4 Omega^2)).
    # The dispersion derivative is evaluated on a smooth radial table and
    # interpolated back: the particle radii are unordered, so differentiating
    # them directly would produce enormous spurious velocities.
    with np.errstate(divide="ignore", invalid="ignore"):
        r_tab = np.geomspace(max(0.02, float(radius.min())),
                             max(0.1, float(radius.max())), 256)
        vc_tab = _circular_velocity(r_tab, np.zeros_like(r_tab), mh, mb, md,
                                    ah, ab, rd, zd, r_halo_max, r_bulge_max)
        dr_tab = np.maximum(1.0e-3 * r_tab, 2.0e-3)
        vcp_tab = _circular_velocity(r_tab + dr_tab, np.zeros_like(r_tab), mh,
                                     mb, md, ah, ab, rd, zd, r_halo_max,
                                     r_bulge_max)
        vcm_tab = _circular_velocity(np.maximum(r_tab - dr_tab, 1.0e-4),
                                     np.zeros_like(r_tab), mh, mb, md, ah, ab,
                                     rd, zd, r_halo_max, r_bulge_max)
        dvc_tab = (vcp_tab - vcm_tab) / (2.0 * dr_tab)
        om_tab = vc_tab / np.maximum(r_tab, 1.0e-6)
        k2_tab = np.maximum(2.0 * om_tab * (om_tab + dvc_tab), 0.0)
        sig_surf_tab = md / (2.0 * np.pi * rd ** 2) * np.exp(-r_tab / rd)
        sig_tab = config.TOOMRE_Q * 3.36 * config.G_KPC * sig_surf_tab / \
            np.maximum(np.sqrt(k2_tab), 1.0e-8)
        dln_sigma2_tab = np.gradient(
            np.log(np.maximum(sig_tab, 1.0e-6) ** 2), np.log(r_tab))
        dln_sigma2 = np.interp(radius, r_tab, dln_sigma2_tab)
        dln_sigma = -radius / rd
    vphi2 = vc2 - sigma_r ** 2 * (
        dln_sigma + np.nan_to_num(dln_sigma2, nan=0.0)
        + 1.0 - kappa2 / np.maximum(4.0 * omega ** 2, 1.0e-12))
    vphi2 = np.where(np.isfinite(vphi2), vphi2, vc2)
    vphi2 = np.maximum(vphi2, 0.04 * vc2)
    vphi = np.sqrt(vphi2)
    sigma_phi = sigma_r * kappa / np.maximum(2.0 * omega, 1.0e-8)
    sigma_z = np.sqrt(np.maximum(np.pi * config.G_KPC * sigma_surface * zd, 0.0))

    phi_az = np.arctan2(pos_d[:, 1], pos_d[:, 0])
    vr = rng.normal(0.0, 1.0, n_disk) * sigma_r
    vphi_rand = vphi + rng.normal(0.0, 1.0, n_disk) * sigma_phi
    vz = rng.normal(0.0, 1.0, n_disk) * sigma_z
    vel_d = np.column_stack((
        vr * np.cos(phi_az) - vphi_rand * np.sin(phi_az),
        vr * np.sin(phi_az) + vphi_rand * np.cos(phi_az),
        vz,
    ))
    mass_d = np.full(n_disk, md / n_disk, dtype=np.float64)
    comp_d = np.full(n_disk, 2, dtype=np.int8)

    pos = np.vstack((pos_h, pos_b, pos_d)).astype(np.float64)
    vel = np.vstack((vel_h, vel_b, vel_d)).astype(np.float64)
    mass = np.concatenate((mass_h, mass_b, mass_d))
    comp = np.concatenate((comp_h, comp_b, comp_d))
    gal = np.full(pos.shape[0], galaxy_id, dtype=np.int8)

    # Zero the centre of mass and the net momentum.
    mt = float(mass.sum())
    com = (mass[:, None] * pos).sum(axis=0) / mt
    pos -= com
    vcm = (mass[:, None] * vel).sum(axis=0) / mt
    vel -= vcm

    return Galaxy(pos=pos.astype(np.float32), vel=vel.astype(np.float32),
                  mass=mass.astype(np.float32), comp=comp, gal=gal,
                  mass_total=mt)


def galaxy_mass_inside_truncation(mass_factor: float,
                                  mh: float = config.M_HALO,
                                  md: float = config.M_DISK,
                                  mb: float = config.M_BULGE,
                                  ah: float = config.A_HALO,
                                  ab: float = config.A_BULGE,
                                  r_halo_max: float = config.R_HALO_MAX,
                                  r_bulge_max: float = config.R_BULGE_MAX) -> float:
    """Total mass used for the Kepler pair orbit."""
    return float(mass_factor * (
        hernquist_mass(r_halo_max, mh, ah)
        + hernquist_mass(r_bulge_max, mb, ab) + md))


def _kepler_velocity(r0: float, rp: float, ecc: float, gm: float
                     ) -> tuple[float, float, float]:
    """Return (v_radial, v_tangential, true_anomaly) for the initial state."""
    if ecc < 1.0e-4:
        return 0.0, float(np.sqrt(gm / r0)), 0.0
    cos_nu = ((1.0 + ecc) * rp / r0 - 1.0) / ecc
    cos_nu = float(np.clip(cos_nu, -1.0, 1.0))
    nu = float(np.arccos(cos_nu))
    p = rp * (1.0 + ecc)
    vscale = np.sqrt(gm / p)
    vr = vscale * ecc * np.sin(nu)
    vt = vscale * (1.0 + ecc * cos_nu)
    return float(vr), float(vt), nu


def make_collision(cfg: config.SimConfig, rng: np.random.Generator
                   ) -> CollisionIC:
    """Build both galaxies and place them on the requested Kepler orbit."""
    n_total = int(cfg.n_particles)
    q = float(cfg.mass_ratio)
    n1 = max(2, int(round(n_total / (1.0 + q))))
    n2 = max(2, n_total - n1)

    g1 = build_galaxy(n1, 1.0, 0, rng)
    g2 = build_galaxy(n2, q, 1, rng)

    m1 = g1.mass_total
    m2 = g2.mass_total
    gm = config.G_KPC * (m1 + m2)
    r0 = float(cfg.separation)
    vr, vt, _ = _kepler_velocity(r0, float(cfg.pericenter),
                                 float(cfg.eccentricity), gm)

    # Rotate each disk (and its velocity field) before placement.
    rot1 = _rotation_z(0.0) @ _rotation_x(float(cfg.incl1))
    rot2 = _rotation_z(30.0) @ _rotation_x(float(cfg.incl2))
    pos1 = g1.pos.astype(np.float64) @ rot1.T
    vel1 = g1.vel.astype(np.float64) @ rot1.T
    pos2 = g2.pos.astype(np.float64) @ rot2.T
    vel2 = g2.vel.astype(np.float64) @ rot2.T

    # Relative position of galaxy 2 with respect to galaxy 1 is (-r0, 0, 0).
    # The pair approaches pericentre and orbits prograde (L_z > 0).
    v_rel = np.array([vr, -vt, 0.0], dtype=np.float64)
    v1 = -(m2 / (m1 + m2)) * v_rel
    v2 = (m1 / (m1 + m2)) * v_rel

    pos1 += np.array([0.5 * r0, 0.0, 0.0])
    pos2 += np.array([-0.5 * r0, 0.0, 0.0])
    vel1 += v1
    vel2 += v2

    pos = np.vstack((pos1, pos2)).astype(np.float32)
    vel = np.vstack((vel1, vel2)).astype(np.float32)
    mass = np.concatenate((g1.mass, g2.mass)).astype(np.float32)
    comp = np.concatenate((g1.comp, g2.comp))
    gal = np.concatenate((g1.gal, g2.gal))

    mt = float(mass.sum())
    com = (mass[:, None] * pos).sum(axis=0) / mt
    pos -= com.astype(np.float32)
    vcm = (mass[:, None] * vel).sum(axis=0) / mt
    vel -= vcm.astype(np.float32)

    # Initial box: 2.4 times the 99th-percentile radius, snapped to a level.
    radius = np.linalg.norm(pos, axis=1)
    target = 2.4 * float(np.percentile(radius, 99.0))
    box = min(config.BOX_LEVELS, key=lambda level: abs(level - target))
    box = max(box, 64.0)
    return CollisionIC(pos=pos, vel=vel, mass=mass, comp=comp, gal=gal,
                       box=float(box), center=np.zeros(3, dtype=np.float32),
                       mass_total=mt, mass_primary=m1)

```
## `src/ics_cosmo.py`

```python
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

```
## `src/gravity_pm.py`

```python
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

```
## `src/integrator.py`

```python
"""Time integration and the two simulation scenes.

Both scenes use kick-drift-kick leapfrog with a global adaptive timestep.
The collision scene is isolated and uses physical units (kpc, km/s, Msun).
The cosmic scene is periodic and integrates the comoving equations of motion
in ln a:

    dx/dln a = p / (a^2 H)
    dp/dln a = -grad(phi) / H

with p = a^2 dx/dt and grad^2(phi) = (3/2) Omega_m delta / a.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

import config
from . import ics_cosmo, ics_galaxy, units
from .gravity_pm import PMGravity, HAVE_NUMBA

try:
    from numba import njit, prange
except Exception:  # pragma: no cover
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda fn: fn
    prange = range


# ---------------------------------------------------------------------------
# Direct summation used by the physics tests
# ---------------------------------------------------------------------------
@njit(cache=True, fastmath=True, nogil=True)
def direct_acceleration(pos: np.ndarray, mass: np.ndarray, G: float,
                        soft: float, out: np.ndarray) -> None:
    n = pos.shape[0]
    eps2 = soft * soft
    for i in range(n):
        ax = 0.0
        ay = 0.0
        az = 0.0
        for j in range(n):
            if i == j:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            r2 = dx * dx + dy * dy + dz * dz + eps2
            factor = G * mass[j] / (r2 * np.sqrt(r2))
            ax += dx * factor
            ay += dy * factor
            az += dz * factor
        out[i, 0] = ax
        out[i, 1] = ay
        out[i, 2] = az


@njit(cache=True, fastmath=True, nogil=True)
def direct_acceleration_serial(pos: np.ndarray, mass: np.ndarray, G: float,
                               soft: float, out: np.ndarray) -> None:
    direct_acceleration(pos, mass, G, soft, out)


@njit(cache=True, fastmath=True, nogil=True)
def direct_energy(pos: np.ndarray, vel: np.ndarray, mass: np.ndarray,
                  G: float, soft: float) -> tuple[float, float]:
    n = pos.shape[0]
    kinetic = 0.0
    potential = 0.0
    for i in range(n):
        vx = vel[i, 0]
        vy = vel[i, 1]
        vz = vel[i, 2]
        kinetic += 0.5 * mass[i] * (vx * vx + vy * vy + vz * vz)
        for j in range(i + 1, n):
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz + soft * soft)
            potential -= G * mass[i] * mass[j] / r
    return kinetic, potential


@njit(cache=True, fastmath=True, nogil=True)
def leapfrog_direct_step(pos: np.ndarray, vel: np.ndarray, mass: np.ndarray,
                         G: float, soft: float, dt: float) -> None:
    acc = np.empty_like(pos)
    direct_acceleration(pos, mass, G, soft, acc)
    vel += 0.5 * dt * acc
    pos += dt * vel
    direct_acceleration(pos, mass, G, soft, acc)
    vel += 0.5 * dt * acc


# ---------------------------------------------------------------------------
# Collision scene
# ---------------------------------------------------------------------------
class CollisionSim:
    """Two-galaxy collision evolved with an isolated PM solver."""

    def __init__(self, cfg: config.SimConfig) -> None:
        self.cfg = cfg
        rng = np.random.default_rng(cfg.seed)
        ic = ics_galaxy.make_collision(cfg, rng)
        self.pos = ic.pos.copy()
        self.vel = ic.vel.copy()
        self.mass = ic.mass.copy()
        self.comp = ic.comp.copy()
        self.gal = ic.gal.copy()
        self.box = ic.box
        self.center = ic.center.copy()
        self.n = cfg.n_grid
        self.gravity = PMGravity(self.n, self.box, False, config.G_KPC,
                                 center=self.center)
        self.acc = self.gravity.accelerations(self.pos, self.mass)
        self.time = 0.0
        self.dt = 0.0
        self.steps = 0
        self.outside = 0
        self.min_sep = float("inf")
        self.t_peri = 0.0
        self._last_box_step = 0
        self.speed = np.linalg.norm(self.vel, axis=1).astype(np.float32)
        k0, w0 = self._energies()
        self.energy0 = k0 + w0
        self._drift = 0.0
        self._virial = 2.0 * k0 / max(abs(w0), 1.0e-30)
        self.stats: dict[str, float] = {}
        self._update_stats()

    # -- helpers -----------------------------------------------------------
    def _energies(self) -> tuple[float, float]:
        phi = self.gravity.potential_at(self.pos)
        kinetic = 0.5 * float(np.sum(self.mass * np.sum(self.vel ** 2, axis=1)))
        potential = 0.5 * float(np.sum(self.mass * phi))
        return kinetic, potential

    def _stars(self, galaxy: int) -> np.ndarray:
        return (self.gal == galaxy) & (self.comp > 0)

    def _iterative_com(self, galaxy: int) -> tuple[np.ndarray, np.ndarray]:
        sel = self._stars(galaxy)
        if not np.any(sel):
            return np.zeros(3), np.zeros(3)
        idx = np.flatnonzero(sel)
        radius = 30.0
        com = np.average(self.pos[idx], axis=0, weights=self.mass[idx])
        for _ in range(3):
            d = np.linalg.norm(self.pos[idx] - com, axis=1)
            keep = d < radius
            if np.count_nonzero(keep) < 10:
                break
            idx = idx[keep]
            com = np.average(self.pos[idx], axis=0, weights=self.mass[idx])
            radius = max(15.0, 0.8 * radius)
        vel = np.average(self.vel[idx], axis=0, weights=self.mass[idx])
        return com.astype(np.float64), vel.astype(np.float64)

    def _choose_dt(self) -> float:
        cell = self.box / self.n
        amax = float(np.max(np.linalg.norm(self.acc, axis=1))) + 1.0e-12
        vmax = float(np.max(self.speed)) + 1.0e-6
        dt = min(config.ETA * np.sqrt(cell / amax), 0.5 * cell / vmax)
        if self.dt > 0.0:
            dt = min(dt, 1.5 * self.dt)
        return float(max(dt, 1.0e-6))

    def _adapt_box(self) -> None:
        if self.steps - self._last_box_step < 50:
            return
        self._last_box_step = self.steps
        mt = float(self.mass.sum())
        com = (self.mass[:, None] * self.pos).sum(axis=0) / mt
        radius = np.linalg.norm(self.pos - com, axis=1)
        target = 2.4 * float(np.percentile(radius, 99.0))
        target = max(target, config.BOX_LEVELS[0])
        level = min(config.BOX_LEVELS, key=lambda value: abs(value - target))
        if level > self.box * 1.1 or level < self.box * 0.9:
            self.box = float(level)
            self.center = com.astype(np.float32)
            self.gravity = PMGravity(self.n, self.box, False, config.G_KPC,
                                     center=self.center)
            self.acc = self.gravity.accelerations(self.pos, self.mass)
        self.outside = int(np.count_nonzero(
            np.any(np.abs(self.pos - self.center) > 0.5 * self.box, axis=1)))

    # -- evolution ---------------------------------------------------------
    def step(self) -> float:
        dt = self._choose_dt()
        self.vel += np.float32(0.5 * dt) * self.acc
        self.pos += np.float32(dt) * self.vel
        self.time += dt
        self.acc = self.gravity.accelerations(self.pos, self.mass)
        self.vel += np.float32(0.5 * dt) * self.acc
        self.dt = dt
        self.steps += 1
        self.speed = np.linalg.norm(self.vel, axis=1).astype(np.float32)
        self._adapt_box()
        self._update_stats()
        return dt

    def _update_stats(self) -> None:
        com1, v1 = self._iterative_com(0)
        com2, v2 = self._iterative_com(1)
        sep = float(np.linalg.norm(com2 - com1))
        dv = v2 - v1
        vrel = float(np.linalg.norm(dv))
        if sep < self.min_sep:
            self.min_sep = sep
            self.t_peri = self.time
        since_peri = max(0.0, self.time - self.t_peri) if sep > 1.02 * self.min_sep else 0.0
        m_small = min(float(self.mass[self.gal == 0].sum()),
                      float(self.mass[self.gal == 1].sum()))
        ln_lambda = 3.0
        t_df_code = 1.17 * sep ** 2 * vrel / (
            config.G_KPC * max(m_small, 1.0) * ln_lambda)
        rho = (float(self.mass.sum()) / (4.0 / 3.0 * np.pi * max(sep, 1.0) ** 3))
        t_dyn = float(units.dynamical_time_gyr(rho))
        if self.steps % 5 == 0 or not hasattr(self, "_drift"):
            kinetic, potential = self._energies()
            energy = kinetic + potential
            self._drift = (energy - self.energy0) / max(abs(self.energy0), 1.0e-30)
            self._virial = 2.0 * kinetic / max(abs(potential), 1.0e-30)
        self.stats = {
            "time_gyr": self.time * config.KPC_KMPS_GYR,
            "dt_myr": self.dt * config.KPC_KMPS_GYR * 1000.0,
            "separation_kpc": sep,
            "relative_velocity": vrel,
            "since_pericenter_gyr": since_peri * config.KPC_KMPS_GYR,
            "t_df_gyr": t_df_code * config.KPC_KMPS_GYR,
            "t_dyn_gyr": t_dyn,
            "energy_drift": self._drift,
            "virial": self._virial,
            "outside": self.outside,
            "box_kpc": self.box,
            "cell_kpc": self.box / self.n,
        }

    # -- output ------------------------------------------------------------
    def snapshot(self) -> dict:
        return {
            "scene": "collision",
            "pos": self.pos.copy(),
            "vel": self.vel.copy(),
            "mass": self.mass,
            "comp": self.comp,
            "gal": self.gal,
            "speed": self.speed.copy(),
            "density": (self.gravity.last_density.copy()
                        if self.gravity.last_density is not None else None),
            "box": self.box,
            "center": self.center.copy(),
            "grid": self.n,
            "time": self.time,
            "a": 1.0,
            "z": 0.0,
            "stats": dict(self.stats),
        }


# ---------------------------------------------------------------------------
# Cosmic web scene
# ---------------------------------------------------------------------------
class CosmicSim:
    """Periodic comoving PM simulation of the cosmic web."""

    def __init__(self, cfg: config.SimConfig) -> None:
        self.cfg = cfg
        box = float(cfg.box if cfg.box is not None else config.BOX_MPCH)
        a_start = 1.0 / (1.0 + float(cfg.z_start))
        ic = ics_cosmo.make_cosmic_ics(cfg.n_particles, cfg.n_grid, box,
                                       a_start=a_start, seed=cfg.seed)
        self.ic = ic
        self.pos = ic.pos.copy()
        self.pvel = ic.pvel.copy()
        self.mass = ic.mass.copy()
        self.box = box
        self.n = cfg.n_grid
        self.a = float(ic.a_start)
        self.time = float(units.cosmic_age_gyr(self.a))
        self.steps = 0
        self.gravity = PMGravity(self.n, self.box, True,
                                 self._g_eff(self.a), smoothing_cells=1.0)
        self.acc = self.gravity.accelerations(self.pos, self.mass)
        self.delta = (ic.delta * np.float32(units.growth_factor(self.a))).astype(np.float32)
        self.delta_rms = float(np.std(self.delta))
        self.node_frac = float(np.mean(self.delta > 10.0))
        self.void_frac = float(np.mean(self.delta < -0.8))
        self.speed = (np.linalg.norm(self.pvel, axis=1) / self.a
                      * config.COSMIC_VEL_KMS).astype(np.float32)
        self.energy0 = self._energy()
        self.stats: dict[str, float] = {}
        self._update_stats()

    def _g_eff(self, a: float) -> float:
        # 4 pi G rho_bar = (3/2) Omega_m / a in H0 = 1 units.  Particles carry
        # Omega_m / N each, so G_eff = (3/2) V / (4 pi a).
        return 1.5 * self.box ** 3 / (4.0 * np.pi * a)

    def _choose_dln_a(self) -> float:
        a = self.a
        e = float(units.hubble_e(a))
        vmax = float(np.max(np.linalg.norm(self.pvel, axis=1)))
        rate = vmax / max(a * a * e, 1.0e-30)  # displacement per unit ln a
        cell = self.box / self.n
        limit = 0.02
        if rate > 0.0:
            limit = min(limit, 0.25 * cell / rate)
        if a < 1.0:
            limit = min(limit, float(np.log(1.0 / a)))
        return float(max(min(limit, 0.05), 1.0e-4))

    def _energy(self) -> float:
        phi = self.gravity.potential_at(self.pos)
        vpec = self.pvel / self.a
        kinetic = 0.5 * float(np.sum(self.mass * np.sum(vpec ** 2, axis=1)))
        potential = 0.5 * float(np.sum(self.mass * phi))
        return kinetic + potential

    def step(self) -> float:
        if self.a >= 1.0 - 1.0e-10:
            return 0.0
        a = self.a
        e = float(units.hubble_e(a))
        dln = self._choose_dln_a()
        self.pvel += self.acc * np.float32(0.5 * dln / e)
        self.pos += (self.pvel / np.float32(a * a * e)) * np.float32(dln)
        np.mod(self.pos, self.box, out=self.pos)
        a_new = a * float(np.exp(dln))
        if a_new > 1.0:
            a_new = 1.0
        self.a = a_new
        self.gravity.set_G(self._g_eff(a_new))
        self.acc = self.gravity.accelerations(self.pos, self.mass)
        e_new = float(units.hubble_e(a_new))
        self.pvel += self.acc * np.float32(0.5 * dln / e_new)
        self.steps += 1
        self.time = float(units.cosmic_age_gyr(a_new))
        self.speed = (np.linalg.norm(self.pvel, axis=1) / self.a
                      * config.COSMIC_VEL_KMS).astype(np.float32)
        if self.steps % 5 == 0:
            self._update_density_stats()
        self._update_stats()
        return dln

    def _update_density_stats(self) -> None:
        dens = self.gravity.last_density
        if dens is None:
            return
        mean = float(np.mean(dens))
        if mean <= 0.0:
            return
        delta = dens / mean - 1.0
        self.delta = delta
        self.delta_rms = float(np.std(delta))
        self.node_frac = float(np.mean(delta > 10.0))
        self.void_frac = float(np.mean(delta < -0.8))

    def _update_stats(self) -> None:
        if self.steps <= 1 or self.steps % 5 == 0 or self.a >= 1.0 - 1.0e-10:
            energy = self._energy()
            self._drift = (energy - self.energy0) / max(abs(self.energy0), 1.0e-30)
            phi = self.gravity.potential_at(self.pos)
            vpec = self.pvel / self.a
            kinetic = 0.5 * float(np.sum(self.mass * np.sum(vpec ** 2, axis=1)))
            potential = 0.5 * float(np.sum(self.mass * phi))
            self._virial = 2.0 * kinetic / max(abs(potential), 1.0e-30)
        self.stats = {
            "a": self.a,
            "z": 1.0 / self.a - 1.0,
            "age_gyr": self.time,
            "delta_rms": self.delta_rms,
            "node_frac": self.node_frac,
            "void_frac": self.void_frac,
            "box_mpch": self.box,
            "cell_mpch": self.box / self.n,
            "energy_drift": self._drift,
            "virial": self._virial,
            "dln_a": self._choose_dln_a(),
        }

    def snapshot(self) -> dict:
        return {
            "scene": "cosmic",
            "pos": self.pos.copy(),
            "vel": self.pvel.copy(),
            "mass": self.mass,
            "comp": np.zeros(self.pos.shape[0], dtype=np.int8),
            "gal": np.zeros(self.pos.shape[0], dtype=np.int8),
            "speed": self.speed.copy(),
            "density": (self.gravity.last_density
                        if self.gravity.last_density is not None else None),
            "box": self.box,
            "center": np.full(3, 0.5 * self.box, dtype=np.float32),
            "grid": self.n,
            "time": self.time,
            "a": self.a,
            "z": 1.0 / self.a - 1.0,
            "stats": dict(self.stats),
        }


def make_simulation(cfg: config.SimConfig) -> CollisionSim | CosmicSim:
    """Create the simulation requested by the configuration."""
    if cfg.scene == "collision":
        return CollisionSim(cfg)
    return CosmicSim(cfg)


def run_steps(sim: CollisionSim | CosmicSim, count: int,
              progress: bool = False) -> None:
    """Advance a simulation by a fixed number of steps (headless helper)."""
    for i in range(count):
        sim.step()
        if progress and (i + 1) % max(1, count // 10) == 0:
            print(f"  step {i + 1}/{count}", flush=True)

```
## `src/render.py`

```python
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

```
## `src/camera.py`

```python
"""Orbiting perspective camera shared by the live view and the renderer."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Camera:
    """Simple yaw-pitch-distance camera looking at a target point."""

    yaw: float = 35.0
    pitch: float = 25.0
    distance: float = 400.0
    target: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    fov: float = 50.0
    min_distance_scale: float = 0.15
    max_distance_scale: float = 8.0
    _extent: float = 300.0

    def reset(self, target: np.ndarray, extent: float) -> None:
        """Place the camera at a pleasant three-quarter view of the scene."""
        self.target = np.asarray(target, dtype=np.float64).copy()
        self._extent = float(extent)
        self.distance = 1.35 * float(extent)
        self.yaw = 35.0
        self.pitch = 25.0

    def eye(self) -> np.ndarray:
        cy, sy = np.cos(np.radians(self.yaw)), np.sin(np.radians(self.yaw))
        cp, sp = np.cos(np.radians(self.pitch)), np.sin(np.radians(self.pitch))
        direction = np.array([cp * cy, cp * sy, sp], dtype=np.float64)
        return self.target + self.distance * direction

    def basis(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (forward, right, up) unit vectors."""
        forward = self.target - self.eye()
        norm = np.linalg.norm(forward)
        forward = forward / max(norm, 1.0e-12)
        world_up = np.array([0.0, 0.0, 1.0])
        right = np.cross(forward, world_up)
        if np.linalg.norm(right) < 1.0e-6:
            right = np.array([1.0, 0.0, 0.0])
        right = right / np.linalg.norm(right)
        up = np.cross(right, forward)
        return forward, right, up

    def projection(self, width: int, height: int
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        forward, right, up = self.basis()
        focal = 0.5 * height / np.tan(0.5 * np.radians(self.fov))
        return self.eye(), forward, right, focal

    def orbit(self, dyaw: float, dpitch: float) -> None:
        self.yaw += float(dyaw)
        self.pitch = float(np.clip(self.pitch + dpitch, -89.0, 89.0))

    def zoom(self, factor: float) -> None:
        lo = self.min_distance_scale * self._extent
        hi = self.max_distance_scale * self._extent
        self.distance = float(np.clip(self.distance * factor, lo, hi))

    def project(self, points: np.ndarray, width: int, height: int
                ) -> tuple[np.ndarray, np.ndarray]:
        """Project points to pixel coordinates; returns (xy, depth)."""
        eye, forward, right, focal = self.projection(width, height)
        rel = np.asarray(points, dtype=np.float64) - eye
        depth = rel @ forward
        safe = np.maximum(depth, 1.0e-9)
        sx = width * 0.5 + (rel @ right) / safe * focal
        sy = height * 0.5 - (rel @ np.cross(right, forward)) / safe * focal
        return np.stack((sx, sy), axis=1), depth

```
## `src/sim_thread.py`

```python
"""Background simulation thread with lock-protected double buffering.

The UI thread never touches the live particle arrays.  The worker publishes a
finished snapshot into the inactive buffer and then flips an index under a
short lock, so a render can keep a stable reference while the next step is
already running.  This keeps the interactive loop responsive even for a
multi-million particle run.
"""

from __future__ import annotations

import threading
import time
from typing import Any

import config
from .integrator import CollisionSim, CosmicSim, make_simulation


class SimulationThread(threading.Thread):
    """Run a :class:`CollisionSim` or :class:`CosmicSim` in the background."""

    def __init__(self, cfg: config.SimConfig,
                 steps_per_second: float = 60.0) -> None:
        super().__init__(name="nbody-simulation", daemon=True)
        self.cfg = cfg
        self.base_rate = float(steps_per_second)
        self.speed = 1.0
        self.sim: CollisionSim | CosmicSim | None = None
        self.steps_done = 0
        self.sim_rate = 0.0
        self.version = 0

        self._stop_event = threading.Event()
        self._pause_event = threading.Event()
        self._ready = threading.Event()
        self._lock = threading.Lock()
        self._error: BaseException | None = None

        # Two snapshot slots: only one is written at a time.
        self._slots: list[dict[str, Any] | None] = [None, None]
        self._front = 0

        n = cfg.n_particles
        if n <= 500_000:
            self.snapshot_every = 1
        elif n <= 2_000_000:
            self.snapshot_every = 3
        else:
            self.snapshot_every = 8

    # -- worker ------------------------------------------------------------
    def run(self) -> None:
        """Create the simulation and advance it until :meth:`stop`."""
        try:
            self.sim = make_simulation(self.cfg)
            self._publish(self.sim)
            self._ready.set()
            last = time.perf_counter()
            accumulator = 0.0
            rate_window = last
            rate_steps = 0
            while not self._stop_event.is_set():
                if self._pause_event.is_set():
                    time.sleep(0.03)
                    last = time.perf_counter()
                    continue
                now = time.perf_counter()
                elapsed = min(now - last, 0.10)
                last = now
                accumulator += self.base_rate * self.speed * elapsed
                count = int(accumulator)
                if count <= 0:
                    time.sleep(0.002)
                    continue
                accumulator -= count
                count = min(count, 64)
                for _ in range(count):
                    if self._stop_event.is_set():
                        break
                    assert self.sim is not None
                    self.sim.step()
                    self.steps_done += 1
                    rate_steps += 1
                    if self.steps_done % self.snapshot_every == 0:
                        self._publish(self.sim)
                self._publish(self.sim)
                if (isinstance(self.sim, CosmicSim)
                        and self.sim.a >= 1.0 - 1.0e-10):
                    # The run has reached the present day; stop burning CPU.
                    self._pause_event.set()
                stamp = time.perf_counter()
                if stamp - rate_window >= 0.5:
                    self.sim_rate = rate_steps / (stamp - rate_window)
                    rate_window = stamp
                    rate_steps = 0
        except BaseException as exc:  # keep the message for the UI thread
            self._error = exc
        finally:
            self._ready.set()

    def _publish(self, sim: CollisionSim | CosmicSim) -> None:
        snapshot = sim.snapshot()
        with self._lock:
            back = 1 - self._front
            self._slots[back] = snapshot
            self._front = back
            self.version += 1

    # -- public API --------------------------------------------------------
    def wait_ready(self, timeout: float | None = None) -> bool:
        """Block until the first snapshot exists or construction failed."""
        return self._ready.wait(timeout)

    def error(self) -> BaseException | None:
        """Return the exception raised by the worker, if any."""
        return self._error

    def latest(self) -> dict[str, Any] | None:
        """Return the newest complete snapshot (immutable after publication)."""
        with self._lock:
            return self._slots[self._front]

    def stop(self) -> None:
        """Ask the worker to finish and wait briefly for it."""
        self._stop_event.set()
        if self.is_alive():
            self.join(timeout=2.0)

    def pause(self) -> None:
        self._pause_event.set()

    def resume(self) -> None:
        self._pause_event.clear()

    def toggle_pause(self) -> bool:
        """Flip the pause state and return ``True`` when now paused."""
        if self._pause_event.is_set():
            self._pause_event.clear()
            return False
        self._pause_event.set()
        return True

    def is_paused(self) -> bool:
        return self._pause_event.is_set()

    def set_speed(self, value: float) -> None:
        """Set the multiplicative speed in the range ``[0.1, 32]``."""
        self.speed = float(min(32.0, max(0.1, value)))

    def finished(self) -> bool:
        """True when a cosmic run has reached a = 1 (today)."""
        snap = self.latest()
        return bool(snap is not None and snap.get("scene") == "cosmic"
                    and snap.get("a", 0.0) >= 1.0 - 1.0e-9)

```
## `src/ui.py`

```python
"""Pygame front end: orbiting camera, live panel, help and frame recording.

The simulation runs in :class:`src.sim_thread.SimulationThread`.  This module
only imports pygame when the interactive mode is actually selected, so the
headless CLI never creates a window or initialises an audio/video subsystem.
"""

from __future__ import annotations

import os
import time
from typing import Any

import numpy as np
import pygame

import config
from .camera import Camera
from .render import Renderer, save_png
from .sim_thread import SimulationThread

COLOR_MODES = ("galaxy", "component", "density", "speed")


class App:
    """Interactive viewer for either simulation scene."""

    def __init__(self, cfg: config.SimConfig,
                 thread: SimulationThread | None = None) -> None:
        self.cfg = cfg
        self.thread = thread or SimulationThread(cfg)
        self.width = cfg.win_w
        self.height = cfg.win_h
        pygame.init()
        self.screen = pygame.display.set_mode(
            (self.width, self.height), pygame.RESIZABLE)
        pygame.display.set_caption("Galactic N-Body — collision / cosmic web")
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("dejavusans", 15)
        self.small = pygame.font.SysFont("dejavusans", 13)
        self.mono = pygame.font.SysFont("dejavusansmono", 14)
        self.big = pygame.font.SysFont("dejavusans", 22, bold=True)

        self.renderer = Renderer(cfg.render_w, cfg.render_h)
        self.hires = Renderer(1920, 1080)
        self.camera = Camera()
        self.camera_ready = False
        self.color_mode = 0
        self.show_dm = self.cfg.scene == "collision"
        self.exposure = 5.0
        self.slab = 0.35
        self.show_help = False
        self.running = True
        self.dragging = False
        self.recording = False
        self.record_index = 0
        self.record_tick = 0
        self.record_dir = os.path.join("images", "recording")
        self.last_frame: np.ndarray | None = None
        self.last_snapshot_id: int | None = None
        self.last_signature: tuple[Any, ...] | None = None
        self.status = "Initialising simulation..."
        self.status_until = time.time() + 4.0
        self.thread.start()

    # -- lifecycle ---------------------------------------------------------
    def run(self) -> None:
        """Enter the pygame event/render loop."""
        while self.running:
            self._events()
            self._update_camera()
            self._draw()
            self.clock.tick(60)
        self.thread.stop()
        pygame.quit()

    # -- input -------------------------------------------------------------
    def _events(self) -> None:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False
            elif event.type == pygame.VIDEORESIZE:
                self.width = max(320, event.w)
                self.height = max(200, event.h)
                self.screen = pygame.display.set_mode(
                    (self.width, self.height), pygame.RESIZABLE)
            elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                self.dragging = True
            elif event.type == pygame.MOUSEBUTTONUP and event.button == 1:
                self.dragging = False
            elif event.type == pygame.MOUSEMOTION and self.dragging:
                self.camera.orbit(event.rel[0] * 0.32,
                                  -event.rel[1] * 0.32)
            elif event.type == pygame.MOUSEWHEEL:
                self.camera.zoom(0.90 if event.y > 0 else 1.11)
            elif event.type == pygame.KEYDOWN:
                self._key(event)

    def _key(self, event: pygame.event.Event) -> None:
        key = event.key
        if key == pygame.K_ESCAPE:
            self.running = False
        elif key == pygame.K_SPACE:
            paused = self.thread.toggle_pause()
            self._toast("Paused" if paused else "Running")
        elif key == pygame.K_h:
            self.show_help = not self.show_help
        elif key == pygame.K_c:
            self.color_mode = (self.color_mode + 1) % len(COLOR_MODES)
            self._toast(f"Colour mode: {COLOR_MODES[self.color_mode]}")
        elif key == pygame.K_d and self.cfg.scene == "collision":
            self.show_dm = not self.show_dm
            self._toast(f"Dark matter: {'on' if self.show_dm else 'off'}")
        elif key == pygame.K_r:
            self._reset_camera()
            self._toast("Camera reset")
        elif key in (pygame.K_PLUS, pygame.K_EQUALS, pygame.K_KP_PLUS):
            self.thread.set_speed(self.thread.speed * 1.5)
            self._toast(f"Speed: {self.thread.speed:.2f}x")
        elif key in (pygame.K_MINUS, pygame.K_KP_MINUS):
            self.thread.set_speed(self.thread.speed / 1.5)
            self._toast(f"Speed: {self.thread.speed:.2f}x")
        elif key == pygame.K_LEFT:
            self.camera.orbit(-4.0, 0.0)
        elif key == pygame.K_RIGHT:
            self.camera.orbit(4.0, 0.0)
        elif key == pygame.K_UP:
            self.camera.orbit(0.0, 3.0)
        elif key == pygame.K_DOWN:
            self.camera.orbit(0.0, -3.0)
        elif key == pygame.K_LEFTBRACKET and self.cfg.scene == "cosmic":
            self.slab = max(0.05, self.slab - 0.05)
            self._toast(f"Slab: {self.slab:.2f} of the box")
        elif key == pygame.K_RIGHTBRACKET and self.cfg.scene == "cosmic":
            self.slab = min(1.00, self.slab + 0.05)
            self._toast(f"Slab: {self.slab:.2f} of the box")
        elif key == pygame.K_s:
            self._save_frame()
        elif key == pygame.K_f:
            self._save_hires()
        elif key == pygame.K_v:
            self.recording = not self.recording
            if self.recording:
                os.makedirs(self.record_dir, exist_ok=True)
            self._toast("Recording: " + ("on" if self.recording else "off"))
        elif key == pygame.K_t:
            self._restart()

    def _toast(self, text: str) -> None:
        self.status = text
        self.status_until = time.time() + 3.0

    # -- simulation helpers ------------------------------------------------
    def _snapshot(self) -> dict[str, Any] | None:
        return self.thread.latest()

    def _reset_camera(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        extent = self._extent(snap)
        self.camera.reset(np.asarray(snap["center"], dtype=np.float64), extent)
        self.camera_ready = True

    def _extent(self, snap: dict[str, Any]) -> float:
        box = float(snap["box"])
        if snap["scene"] == "cosmic":
            return 0.82 * box
        return 0.28 * box

    def _update_camera(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        if not self.camera_ready:
            self._reset_camera()
        # The collision box grows adaptively; keep the camera framing while
        # the user is not actively zooming by updating only the extent.
        if snap["scene"] == "collision":
            target = np.asarray(snap["center"], dtype=np.float64)
            if np.linalg.norm(target - self.camera.target) > 1.0e-3:
                self.camera.target = target
        else:
            self.camera.target = np.asarray(snap["center"], dtype=np.float64)

    def _restart(self) -> None:
        self.thread.stop()
        self.cfg.seed += 1
        self.thread = SimulationThread(self.cfg)
        self.thread.start()
        self.camera_ready = False
        self.last_frame = None
        self.last_snapshot_id = None
        self.last_signature = None
        self._toast(f"Restarted with seed {self.cfg.seed}")

    # -- rendering ---------------------------------------------------------
    def _make_frame(self, snap: dict[str, Any]) -> np.ndarray:
        signature = (self.thread.version, self.color_mode, self.show_dm,
                     round(self.exposure, 3), round(self.slab, 3),
                     tuple(np.round(self.camera.target, 3)),
                     round(self.camera.yaw, 2), round(self.camera.pitch, 2),
                     round(self.camera.distance, 3))
        if self.last_frame is not None and signature == self.last_signature:
            return self.last_frame
        frame = self.renderer.render(
            snap, self.camera, self.color_mode, self.show_dm,
            exposure=self.exposure, slab=self.slab)
        self.last_frame = frame
        self.last_signature = signature
        if self.recording:
            self.record_tick += 1
            if self.record_tick % 2 == 0:
                os.makedirs(self.record_dir, exist_ok=True)
                path = os.path.join(self.record_dir,
                                    f"frame_{self.record_index:05d}.png")
                save_png(frame, path)
                self.record_index += 1
        return frame

    def _save_frame(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        frame = self.renderer.render(snap, self.camera, self.color_mode,
                                     self.show_dm, self.exposure, self.slab)
        path = os.path.join("images",
                            time.strftime("render_%Y%m%d_%H%M%S.png"))
        save_png(frame, path)
        self._toast(f"Saved {path}")

    def _save_hires(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        frame = self.hires.render(snap, self.camera, self.color_mode,
                                  self.show_dm, self.exposure, self.slab)
        path = os.path.join("images",
                            time.strftime("frame_1920_%Y%m%d_%H%M%S.png"))
        save_png(frame, path)
        self._toast(f"Saved {path} (1920x1080)")

    # -- drawing -----------------------------------------------------------
    def _draw(self) -> None:
        self.screen.fill((4, 6, 12))
        snap = self._snapshot()
        if snap is not None:
            if self.status == "Initialising simulation...":
                self.status_until = 0.0
            frame = self._make_frame(snap)
            source = pygame.surfarray.make_surface(
                np.transpose(frame, (1, 0, 2)))
            scaled = pygame.transform.smoothscale(
                source, (self.width, self.height))
            self.screen.blit(scaled, (0, 0))
        else:
            self._text("Initialising simulation...", 30, 30, self.big,
                       (220, 230, 255))
            if self.thread.error() is not None:
                self._text(str(self.thread.error()), 30, 70, self.font,
                           (255, 110, 110))

        self._draw_panel(snap)
        self._draw_hint()
        if self.show_help:
            self._draw_help()
        if time.time() < self.status_until:
            self._draw_toast()
        pygame.display.flip()

    def _text(self, text: str, x: float, y: float,
              font: pygame.font.Font, color: tuple[int, int, int]) -> None:
        self.screen.blit(font.render(text, True, color), (int(x), int(y)))

    def _draw_panel(self, snap: dict[str, Any] | None) -> None:
        lines: list[tuple[str, str]] = []
        lines.append(("N particles", f"{self.cfg.n_particles:,}"))
        lines.append(("Grid", f"{self.cfg.n_grid}³"))
        lines.append(("FPS", f"{self.clock.get_fps():5.1f}"))
        lines.append(("Sim steps/s", f"{self.thread.sim_rate:5.1f}"))
        lines.append(("Speed", f"{self.thread.speed:5.2f}x"))
        lines.append(("State", "PAUSED" if self.thread.is_paused() else "running"))
        if snap is not None:
            st = snap["stats"]
            if snap["scene"] == "collision":
                lines.append(("Time", f"{st['time_gyr']:.3f} Gyr"))
                lines.append(("Timestep", f"{st['dt_myr']:.2f} Myr"))
                lines.append(("Separation", f"{st['separation_kpc']:.1f} kpc"))
                lines.append(("Relative v", f"{st['relative_velocity']:.0f} km/s"))
                lines.append(("Since pericentre",
                              f"{st['since_pericenter_gyr']:.3f} Gyr"))
                lines.append(("t_dyn", f"{st['t_dyn_gyr']:.3f} Gyr"))
                lines.append(("t_df", f"{st['t_df_gyr']:.1f} Gyr"))
                lines.append(("Box / cell",
                              f"{st['box_kpc']:.0f} / {st['cell_kpc']:.2f} kpc"))
                lines.append(("Outside box", f"{int(st['outside'])}"))
            else:
                lines.append(("Scale factor a", f"{st['a']:.4f}"))
                lines.append(("Redshift z", f"{st['z']:.3f}"))
                lines.append(("Age", f"{st['age_gyr']:.3f} Gyr"))
                lines.append(("delta_rms", f"{st['delta_rms']:.3f}"))
                lines.append(("Nodes > 10", f"{100.0 * st['node_frac']:.2f} %"))
                lines.append(("Voids < -0.8", f"{100.0 * st['void_frac']:.1f} %"))
                lines.append(("Box / cell",
                              f"{st['box_mpch']:.0f} / {st['cell_mpch']:.2f} Mpc/h"))
                lines.append(("Slab", f"{self.slab:.2f} box"))
            if snap["scene"] == "collision":
                lines.append(("E drift", f"{100.0 * st['energy_drift']:+.3f} %"))
            lines.append(("Virial 2K/|W|", f"{st['virial']:.3f}"))
        if snap is not None and snap["scene"] == "collision":
            lines.append(("Colour", COLOR_MODES[self.color_mode]))
            lines.append(("Dark matter", "on" if self.show_dm else "off"))
        if self.recording:
            lines.append(("Recording", f"{self.record_index} frames"))

        panel_w = 330
        panel_h = 14 + 19 * len(lines)
        surface = pygame.Surface((panel_w, panel_h), pygame.SRCALPHA)
        surface.fill((6, 10, 20, 205))
        pygame.draw.rect(surface, (60, 80, 120, 180),
                         pygame.Rect(0, 0, panel_w, panel_h), 1)
        self.screen.blit(surface, (12, 12))
        self._text(self._scene_title(snap), 22, 18, self.big, (235, 242, 255))
        y = 48
        for label, value in lines:
            self._text(label, 22, y, self.small, (145, 168, 200))
            value_surface = self.small.render(value, True, (225, 235, 250))
            self.screen.blit(value_surface,
                             (panel_w - 18 - value_surface.get_width(), y))
            y += 19

    def _scene_title(self, snap: dict[str, Any] | None) -> str:
        if snap is None:
            return "GALACTIC N-BODY"
        if snap["scene"] == "collision":
            return "GALAXY COLLISION"
        return "COSMIC WEB"

    def _draw_hint(self) -> None:
        hint = ("SPACE pause   C colour   D dark matter   S save PNG   "
                "F 1920x1080   V record   T restart   H help")
        if self.cfg.scene == "cosmic":
            hint += "   [ ] slab"
        self._text(hint, 16, self.height - 26, self.small, (135, 155, 185))

    def _draw_toast(self) -> None:
        width = min(self.width - 40, 620)
        surface = pygame.Surface((width, 36), pygame.SRCALPHA)
        surface.fill((16, 28, 48, 220))
        pygame.draw.rect(surface, (90, 140, 210, 220),
                         pygame.Rect(0, 0, width, 36), 1)
        self.screen.blit(surface, (20, self.height - 72))
        self._text(self.status, 32, self.height - 66, self.small,
                   (215, 230, 250))

    def _draw_help(self) -> None:
        rows = [
            ("Mouse drag", "orbit the camera"),
            ("Mouse wheel", "zoom in / out"),
            ("Arrow keys", "orbit / change inclination"),
            ("SPACE", "pause or resume"),
            ("C", "cycle colour mode"),
            ("D", "show / hide dark matter (collision)"),
            ("+ / -", "simulation speed"),
            ("S", "save the current preview PNG"),
            ("F", "save a 1920x1080 PNG"),
            ("V", "start / stop frame recording"),
            ("[ / ]", "change the cosmic slab thickness"),
            ("R", "reset the camera"),
            ("T", "restart with a new seed"),
            ("H", "close this help"),
            ("ESC", "quit"),
        ]
        width = min(620, self.width - 60)
        height = 30 + 26 * len(rows)
        x = (self.width - width) // 2
        y = max(20, (self.height - height) // 2)
        surface = pygame.Surface((width, height), pygame.SRCALPHA)
        surface.fill((6, 10, 20, 235))
        pygame.draw.rect(surface, (80, 120, 180, 230),
                         pygame.Rect(0, 0, width, height), 1)
        self.screen.blit(surface, (x, y))
        self._text("CONTROLS", x + 20, y + 10, self.big, (235, 242, 255))
        for i, (key, action) in enumerate(rows):
            yy = y + 42 + 26 * i
            self._text(key, x + 24, yy, self.mono, (255, 205, 120))
            self._text(action, x + 190, yy, self.font, (205, 220, 240))

```
## `src/analysis.py`

```python
"""Post-processing: cosmic power spectrum and collision radial profiles.

Both diagnostics are intentionally modest.  The PM grid sets the smallest
available scale, so the power spectrum is only meaningful above a few cells
and the radial profile is compared with the *analytic* Hernquist halo outside
the PM resolution limit rather than with a fit to the stars.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np

import config
from . import ics_cosmo, ics_galaxy, units
from .integrator import CollisionSim, CosmicSim, make_simulation


# ---------------------------------------------------------------------------
# Cosmic power spectrum
# ---------------------------------------------------------------------------
def cosmic_power_spectrum(sim: CosmicSim) -> tuple[np.ndarray, np.ndarray,
                                                   np.ndarray]:
    """Return (k, P_measured, P_linear) in (h/Mpc, (Mpc/h)^3, (Mpc/h)^3)."""
    density = sim.gravity.last_density
    if density is None:
        density = sim.gravity.deposit(sim.pos, sim.mass)
    mean = float(np.mean(density, dtype=np.float64))
    if mean <= 0.0:
        raise RuntimeError("empty density grid")
    delta = density.astype(np.float64) / mean - 1.0
    n = delta.shape[0]
    n_cells = n ** 3
    box = float(sim.box)
    h = box / n
    volume = box ** 3

    from scipy import fft as sfft
    delta_k = sfft.rfftn(delta, workers=-1)
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=h)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n, d=h)
    k2 = (kx[:, None, None] ** 2 + kx[None, :, None] ** 2
          + kz[None, None, :] ** 2)
    k = np.sqrt(k2)
    # CIC assignment window, deconvolved mode by mode.
    window = (np.sinc(kx[:, None, None] * h / (2.0 * np.pi))
              * np.sinc(kx[None, :, None] * h / (2.0 * np.pi))
              * np.sinc(kz[None, None, :] * h / (2.0 * np.pi)))
    # numpy's forward FFT is unnormalised: delta_k = N_cells * delta_tilde.
    # Therefore P(k) = V |delta_tilde|^2 = V |delta_k|^2 / N_cells^2.
    power = (np.abs(delta_k) ** 2 * volume
             / (n_cells ** 2 * np.maximum(window, 1.0e-6) ** 2))

    actual_n = sim.pos.shape[0]
    shot = volume / max(actual_n, 1)
    power = np.maximum(power - shot, 1.0e-30)

    kmin = 2.0 * np.pi / box
    kmax = np.pi * n / box
    edges = np.geomspace(max(kmin, 1.0e-6), kmax, 33)
    k_mid = np.sqrt(edges[:-1] * edges[1:])
    p_bin = np.zeros_like(k_mid)
    for i in range(k_mid.size):
        mask = (k >= edges[i]) & (k < edges[i + 1])
        if np.any(mask):
            p_bin[i] = float(np.mean(power[mask]))
    valid = p_bin > 0.0
    k_out = k_mid[valid]
    p_out = p_bin[valid]
    d = float(units.growth_factor(sim.a))
    p_lin = ics_cosmo.linear_power(k_out) * d * d
    return k_out, p_out, p_lin


def plot_cosmic_power(sim: CosmicSim, path: str) -> str:
    """Write a P(k) comparison plot and return the output path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k, p_sim, p_lin = cosmic_power_spectrum(sim)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=140)
    ax.set_facecolor("#080b14")
    fig.patch.set_facecolor("#080b14")
    ax.loglog(k, p_lin, color="#62b7ff", lw=1.6,
              label=r"linear $\Lambda$CDM $\times D(a)^2$")
    ax.loglog(k, p_sim, color="#ffb15e", lw=1.5, marker="o", ms=2.5,
              label="PM simulation")
    ax.set_xlabel(r"$k\;[h\,\mathrm{Mpc}^{-1}]$", color="#cfd9ea")
    ax.set_ylabel(r"$P(k)\;[(h^{-1}\mathrm{Mpc})^3]$", color="#cfd9ea")
    ax.set_title(f"Cosmic web at a = {sim.a:.3f}  (z = {1.0 / sim.a - 1.0:.2f})",
                 color="#eef4ff")
    ax.grid(True, which="both", color="#1d2a42", lw=0.6)
    ax.tick_params(colors="#adbbd1")
    for spine in ax.spines.values():
        spine.set_color("#2a3a57")
    leg = ax.legend(frameon=False, labelcolor="#d8e4f7")
    for text in leg.get_texts():
        text.set_color("#d8e4f7")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Collision radial profile
# ---------------------------------------------------------------------------
def collision_radial_profile(sim: CollisionSim) -> tuple[np.ndarray, ...]:
    """Measure the primary DM halo profile and compare with Hernquist."""
    sel = (sim.gal == 0) & (sim.comp == 0)
    if not np.any(sel):
        raise RuntimeError("no dark-matter particles in galaxy 0")
    mass = sim.mass[sel].astype(np.float64)
    pos = sim.pos[sel].astype(np.float64)
    com = np.average(pos, axis=0, weights=mass)
    radius = np.linalg.norm(pos - com, axis=1)
    m_sampled = float(mass.sum())

    # The PM force is only resolved beyond roughly two cells; below that the
    # sampled profile is flattened by the grid and is not a physical result.
    r_min = max(2.0, 2.0 * float(sim.box) / float(sim.n))
    edges = np.geomspace(r_min, config.R_HALO_MAX, 28)
    shell_mass, _ = np.histogram(radius, bins=edges, weights=mass)
    volume = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    density = shell_mass / volume
    r_mid = np.sqrt(edges[:-1] * edges[1:])
    # The particles sample a Hernquist profile truncated at R_HALO_MAX, so
    # rescale the sampled mass back to the full (untruncated) halo mass before
    # comparing with the analytic profile.
    truncated_fraction = float(ics_galaxy.hernquist_mass(
        config.R_HALO_MAX, 1.0, config.A_HALO))
    m_full = m_sampled / max(truncated_fraction, 1.0e-12)
    analytic = np.asarray(ics_galaxy.hernquist_density(
        r_mid, m_full, config.A_HALO))
    # The histogram places each shell's mass at the shell midpoint; subtract
    # half a shell to approximate M(<r_mid) rather than M(<r_hi).
    enclosed = np.cumsum(shell_mass) - 0.5 * shell_mass
    analytic_mass = np.asarray(ics_galaxy.hernquist_mass(
        r_mid, m_full, config.A_HALO))
    return r_mid, density, analytic, enclosed, analytic_mass


def plot_collision_profile(sim: CollisionSim, path: str) -> str:
    """Write density and enclosed-mass diagnostics and return the path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r, rho, rho_h, menc, menc_h = collision_radial_profile(sim)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), dpi=140)
    fig.patch.set_facecolor("#080b14")
    for ax in axes:
        ax.set_facecolor("#080b14")
        ax.grid(True, which="both", color="#1d2a42", lw=0.6)
        ax.tick_params(colors="#adbbd1")
        for spine in ax.spines.values():
            spine.set_color("#2a3a57")

    axes[0].loglog(r, rho_h, color="#62b7ff", lw=1.5,
                   label="Hernquist model")
    axes[0].loglog(r, rho, color="#ffb15e", lw=1.2, marker="o", ms=2.5,
                   label="simulation")
    axes[0].set_xlabel(r"$r\;[\mathrm{kpc}]$", color="#cfd9ea")
    axes[0].set_ylabel(r"$\rho\;[M_\odot\,\mathrm{kpc}^{-3}]$",
                       color="#cfd9ea")
    axes[0].set_title("Dark-matter density", color="#eef4ff")

    axes[1].loglog(r, menc_h, color="#62b7ff", lw=1.5, label="Hernquist")
    axes[1].loglog(r, menc, color="#ffb15e", lw=1.2, marker="o", ms=2.5,
                   label="simulation")
    axes[1].set_xlabel(r"$r\;[\mathrm{kpc}]$", color="#cfd9ea")
    axes[1].set_ylabel(r"$M(<r)\;[M_\odot]$", color="#cfd9ea")
    axes[1].set_title("Enclosed mass", color="#eef4ff")

    for ax in axes:
        leg = ax.legend(frameon=False)
        for text in leg.get_texts():
            text.set_color("#d8e4f7")
    fig.suptitle(
        f"Galaxy collision at t = {sim.stats['time_gyr']:.3f} Gyr, "
        f"separation = {sim.stats['separation_kpc']:.1f} kpc",
        color="#eef4ff")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def advance(sim: CollisionSim | CosmicSim, cfg: config.SimConfig) -> int:
    """Advance a simulation using the same defaults as the headless CLI."""
    if cfg.steps is not None:
        count = int(cfg.steps)
        for _ in range(count):
            sim.step()
        return count
    cap = cfg.default_steps()
    count = 0
    if isinstance(sim, CollisionSim):
        while sim.stats["time_gyr"] < 0.55 and count < cap:
            sim.step()
            count += 1
    else:
        while sim.a < 1.0 - 1.0e-10 and count < cap:
            sim.step()
            count += 1
    return count


def analyze(cfg: config.SimConfig,
            sim: CollisionSim | CosmicSim | None = None) -> str:
    """Run the requested simulation (if needed) and save its diagnostic plot."""
    if sim is None:
        sim = make_simulation(cfg)
        advance(sim, cfg)
    if cfg.scene == "cosmic":
        assert isinstance(sim, CosmicSim)
        output = os.path.join("images", "analysis_cosmic.png")
        path = plot_cosmic_power(sim, output)
        k, p_sim, p_lin = cosmic_power_spectrum(sim)
        print(f"  P(k): {k.size} bins, a={sim.a:.4f}, z={1.0 / sim.a - 1.0:.3f}")
    else:
        assert isinstance(sim, CollisionSim)
        output = os.path.join("images", "analysis_collision.png")
        path = plot_collision_profile(sim, output)
        print(f"  profile: t={sim.stats['time_gyr']:.3f} Gyr, "
              f"separation={sim.stats['separation_kpc']:.1f} kpc")
    print(f"  saved {path}")
    return path

```
## `tests/test_physics.py`

```python
"""Physics regression tests for the galactic N-body project.

The suite is intentionally small and fast: it checks the mass assignment,
the PM force normalisation, momentum conservation, the leapfrog integrator,
the Hernquist model, the linear growth factor and the cosmological age.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

import config
from src import ics_cosmo, ics_galaxy, units
from src.gravity_pm import PMGravity
from src.integrator import direct_energy, leapfrog_direct_step


# ---------------------------------------------------------------------------
# Particle-mesh gravity
# ---------------------------------------------------------------------------
def test_cic_mass_conservation() -> None:
    """CIC deposit preserves the total mass in periodic and isolated boxes."""
    rng = np.random.default_rng(11)
    n = 32
    box = 32.0
    pos = rng.random((5000, 3), dtype=np.float32) * np.float32(box)
    mass = rng.random(5000, dtype=np.float32)
    total = float(mass.sum())
    for periodic in (True, False):
        solver = PMGravity(n, box, periodic, 1.0)
        grid = solver.deposit(pos, mass)
        assert float(grid.sum()) == pytest.approx(total, rel=1.0e-6)


def test_point_mass_force_beyond_three_cells() -> None:
    """The isolated PM force matches G M / r^2 to 5% for r >= 4 cells."""
    n = 64
    box = 64.0
    center = np.array([32.0, 32.0, 32.0])
    solver = PMGravity(n, box, False, 1.0, center=center)
    source = np.array([[32.0, 32.0, 32.0]], dtype=np.float32)
    distances = np.array([4.0, 5.0, 6.0, 8.0])
    targets = np.stack((32.0 + distances, np.full(4, 32.0),
                        np.full(4, 32.0)), axis=1).astype(np.float32)
    pos = np.vstack((source, targets))
    mass = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    acc = solver.accelerations(pos, mass)
    measured = np.abs(acc[1:, 0].astype(np.float64))
    expected = 1.0 / distances ** 2
    relative = np.abs(measured - expected) / expected
    assert np.all(relative < 0.05), relative


def test_periodic_momentum_conservation() -> None:
    """CIC deposit and interpolation conserve total momentum in a periodic box."""
    rng = np.random.default_rng(7)
    n = 32
    box = 32.0
    pos = rng.random((20000, 3), dtype=np.float32) * np.float32(box)
    mass = np.ones(20000, dtype=np.float32)
    solver = PMGravity(n, box, True, 1.0)
    acc = solver.accelerations(pos, mass)
    net = np.sum(mass[:, None] * acc, axis=0, dtype=np.float64)
    scale = float(np.sum(mass * np.linalg.norm(acc, axis=1)))
    assert np.linalg.norm(net) / scale < 1.0e-5


# ---------------------------------------------------------------------------
# Integrator
# ---------------------------------------------------------------------------
def test_leapfrog_energy_and_reversibility() -> None:
    """KDK leapfrog conserves energy over ten orbits and is time reversible."""
    g = 1.0
    mass = np.array([0.5, 0.5])
    pos0 = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    vel0 = np.array([[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]])
    dt = 2.0 * np.pi / 4000.0

    pos = pos0.copy()
    vel = vel0.copy()
    k0, w0 = direct_energy(pos, vel, mass, g, 0.0)
    e0 = k0 + w0
    for _ in range(40000):
        leapfrog_direct_step(pos, vel, mass, g, 0.0, dt)
    k1, w1 = direct_energy(pos, vel, mass, g, 0.0)
    assert abs((k1 + w1 - e0) / e0) < 1.0e-3

    pos = pos0.copy()
    vel = vel0.copy()
    for _ in range(500):
        leapfrog_direct_step(pos, vel, mass, g, 0.0, dt)
    vel = -vel
    for _ in range(500):
        leapfrog_direct_step(pos, vel, mass, g, 0.0, dt)
    assert np.max(np.abs(pos - pos0)) < 1.0e-8
    assert np.max(np.abs(vel + vel0)) < 1.0e-8


# ---------------------------------------------------------------------------
# Hernquist model
# ---------------------------------------------------------------------------
def test_hernquist_mass_velocity_and_virial() -> None:
    """Enclosed mass, circular speed and the virial ratio are consistent."""
    mass = 1.0e12
    scale = 30.0
    radii = np.array([0.25, 0.5, 1.0, 2.0, 5.0]) * scale
    for radius in radii:
        numerical, _ = quad(
            lambda r: 4.0 * np.pi * r * r
            * float(ics_galaxy.hernquist_density(r, mass, scale)),
            1.0e-6, float(radius), limit=200)
        analytic = float(ics_galaxy.hernquist_mass(radius, mass, scale))
        assert abs(numerical - analytic) / analytic < 0.05

        dr = 1.0e-4 * radius
        phi_p = float(ics_galaxy.hernquist_potential(radius + dr, mass, scale))
        phi_m = float(ics_galaxy.hernquist_potential(radius - dr, mass, scale))
        vc2_numeric = radius * (phi_p - phi_m) / (2.0 * dr)
        vc2_analytic = config.G_KPC * analytic / radius
        assert abs(vc2_numeric - vc2_analytic) / vc2_analytic < 0.05

    rng = np.random.default_rng(5)
    n = 2000
    r_max = 100.0
    pos = ics_galaxy.sample_hernquist(n, mass, scale, r_max, rng)
    radii = np.linalg.norm(pos, axis=1)
    sigma = ics_galaxy.jeans_sigma(radii, mass, scale, r_max)
    vel = rng.normal(0.0, 1.0, (n, 3)) * sigma[:, None]
    m_particle = np.full(n, float(ics_galaxy.hernquist_mass(
        r_max, mass, scale)) / n)
    kinetic, potential = direct_energy(pos, vel, m_particle,
                                       config.G_KPC, 0.0)
    virial = 2.0 * kinetic / abs(potential)
    assert 0.7 < virial < 1.3, virial


# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------
def test_growth_factor_normalisation_and_eds() -> None:
    """D(1) is about 0.78 unnormalised, D(1)=1 normalised, and EdS gives D=a."""
    raw = float(units.growth_factor(1.0, normalized=False))
    assert raw == pytest.approx(0.779, abs=0.02)
    assert float(units.growth_factor(1.0)) == pytest.approx(1.0, rel=1.0e-6)
    for a in (0.25, 0.5, 1.0):
        d_eds = float(units.growth_factor(a, 1.0, 0.0, normalized=False))
        assert d_eds == pytest.approx(a, rel=1.0e-4)


def test_cosmic_age() -> None:
    """The age of the default flat LambdaCDM universe is 13.5-13.9 Gyr."""
    age = float(units.cosmic_age_gyr(1.0))
    assert 13.5 <= age <= 13.9, age


def test_sigma8_normalisation() -> None:
    """The Eisenstein-Hu spectrum is normalised to sigma_8 = 0.8."""
    sigma = ics_cosmo.sigma_r(lambda k: ics_cosmo.linear_power(k), 8.0)
    assert sigma == pytest.approx(config.SIGMA8, rel=0.02)

```

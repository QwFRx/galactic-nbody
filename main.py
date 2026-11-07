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

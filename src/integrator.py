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

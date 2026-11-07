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

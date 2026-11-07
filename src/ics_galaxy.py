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

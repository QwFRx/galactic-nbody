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

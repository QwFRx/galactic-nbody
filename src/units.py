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

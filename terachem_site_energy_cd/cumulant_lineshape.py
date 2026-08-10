"""Second-order cumulant lineshape from the QM/MM gap correlation function.

Replaces the phenomenological T2* Lorentzian. The manuscript currently applies
T2* = 60 fs as homogeneous dephasing AND treats the site-energy spread as static
detuning -- but those are the same physical motion, so counting both double-
counts it. The cumulant separates them by timescale instead:

    fast (tau_c ~ tens of fs, resolved in the 4 fs / 8 ps segments)
        -> homogeneous dephasing through g(t)
    slow (only sampled over the 1 ns ensemble)
        -> static inhomogeneous detuning, averaged frame by frame

Classical high-temperature second-order cumulant:

    g(t) = int_0^t dtau (t - tau) C(tau) / hbar^2      [C in angular freq^2]
    L(t) = exp(-g(t))
    I(w) = (1/pi) Re int_0^inf dt exp(i (w - w0) t) L(t)

`amplitude_scale` rescales C(t) without touching its shape: the archived C(t)
comes from a first-order Delta_q . V electrostatic gap, which underestimates the
QM/MM detuning, while the correlation TIME is far more transferable than the
amplitude.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = ["load_gap_correlation", "cumulant_g", "lineshape_from_g", "fwhm_of"]

TWO_PI_C = 2.0 * np.pi * 2.99792458e-5      # cm^-1 -> rad/fs


def load_gap_correlation(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (lag_fs, C_difference) in fs and cm^-2 from a correlations.csv."""
    rows = np.genfromtxt(path, delimiter=",", names=True)
    return np.asarray(rows["lag_fs"], float), np.asarray(rows["C_difference_cm2"], float)


def cumulant_g(lag_fs: np.ndarray, c_cm2: np.ndarray,
               amplitude_scale: float = 1.0) -> np.ndarray:
    """g(t), dimensionless, from a classical correlation function in cm^-2."""
    c_rad2 = c_cm2 * amplitude_scale * TWO_PI_C**2      # (rad/fs)^2
    t = np.asarray(lag_fs, float)
    g = np.zeros_like(t)
    for i in range(1, len(t)):
        tau = t[: i + 1]
        g[i] = np.trapezoid((t[i] - tau) * c_rad2[: i + 1], tau)
    return g


def lineshape_from_g(lag_fs, g, grid_cm, centre_cm, extra_hwhm_cm=0.0):
    """Absorption lineshape I(w) by half-Fourier transform of exp(-g(t))."""
    t = np.asarray(lag_fs, float)
    decay = np.exp(-np.asarray(g, float))
    if extra_hwhm_cm > 0.0:
        decay = decay * np.exp(-TWO_PI_C * extra_hwhm_cm * t)
    grid = np.asarray(grid_cm, float)
    dw = (grid - centre_cm) * TWO_PI_C                  # rad/fs
    phase = np.exp(1j * dw[:, None] * t[None, :])
    integral = np.trapezoid(phase * decay[None, :], t, axis=1)
    return np.real(integral) / np.pi


def fwhm_of(grid_cm, y) -> float:
    y = np.asarray(y, float)
    peak = y.max()
    above = np.where(y >= 0.5 * peak)[0]
    if len(above) < 2:
        return float("nan")
    return float(grid_cm[above[-1]] - grid_cm[above[0]])

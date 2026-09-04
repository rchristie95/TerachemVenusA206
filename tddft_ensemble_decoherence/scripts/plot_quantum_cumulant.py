#!/usr/bin/env python3
"""Figure: quantum vs classical cumulant band shape against the measured band.

Left  -- the spectral density the 1 ps record actually measures, with kT marked,
         showing that only 8% of the variance lies where the old partition cut.
Right -- the bands. Both curves carry ONE parameter, sigma_static, fixed by the
         measured FWHM and required to lie inside the range the 1 ns ensemble
         permits. The blue/red asymmetry is then free, and it is the only thing
         separating the quantum and classical treatments.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v8dense_quantum_cumulant import (                                # noqa: E402
    dense_series, power_spectrum, cumulant_g, lineshape, solve_static,
    band_metrics, TWO_PI_C, KT_CM, EXP_PEAK_CM, EXP_FWHM_CM, EXP_ASYM)

DIG = REPO / "terachem_site_energy_cd" / "results" / "digitised" / \
    "venus_absorption_digitised.npz"


def main() -> None:
    a, b = dense_series()
    beta = 1.0 / (KT_CM * TWO_PI_C)
    t = np.arange(0.0, 400.0, 0.5)
    grid = np.arange(16500.0, 22500.0, 2.0)

    d = np.load(DIG)
    nu_exp = 1e7 / d["wavelength_nm"]
    o = np.argsort(nu_exp)
    nu_exp, y_exp = nu_exp[o], d["normalised"][o]
    y_exp = y_exp / y_exp.max()

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # ---- left: the measured spectral density ------------------------------
    for tag, x, col in (("site A", a, "#1f4e79"), ("site B", b, "#c0504d")):
        om, s, var = power_spectrum(x)
        nu = om / TWO_PI_C
        w = s.copy(); w[1:-1] *= 2.0
        j = 0.5 * beta * om * w / TWO_PI_C          # J(omega), cm^-1
        k = np.ones(5) / 5.0
        ax[0].plot(nu[1:], np.convolve(j[1:], k, "same"), color=col, lw=1.1,
                   label=tag)
    ax[0].axvline(KT_CM, color="0.35", ls="--", lw=1.0)
    ax[0].text(KT_CM * 1.15, ax[0].get_ylim()[1] * 0.92,
               "$k_BT$\n(the old cut)", fontsize=8, color="0.35", va="top")
    ax[0].axvspan(1500, 2200, color="0.85", zorder=0)
    ax[0].text(1850, ax[0].get_ylim()[1] * 0.58,
               "38% of $\\sigma^2$, 40% of $\\lambda$.\n"
               "$S\\approx0.19$, so it cannot\nmove the half-max points\n"
               "-- it makes the wing", fontsize=7.2, ha="center", color="0.3")
    ax[0].set_xlim(0, 2600)
    ax[0].set_xlabel("wavenumber (cm$^{-1}$)")
    ax[0].set_ylabel("spectral density $J(\\omega)$  (cm$^{-1}$)")
    ax[0].set_title("the spectral density, measured with no partition",
                    fontsize=10)
    ax[0].legend(frameon=False, fontsize=8)

    # ---- right: the bands --------------------------------------------------
    ax[1].plot(nu_exp - EXP_PEAK_CM, y_exp, color="k", lw=2.0,
               label="experiment (Nguyen Fig 2e)")
    styles = {("site A", True): ("#1f4e79", "-"), ("site B", True): ("#4f81bd", "-"),
              ("site A", False): ("#c0504d", "--")}
    for tag, x in (("site A", a), ("site B", b)):
        om, s, var = power_spectrum(x)
        for quantum in (True, False):
            if (tag, quantum) not in styles:
                continue
            g = cumulant_g(t, om, s, beta, quantum=quantum)
            ss = solve_static(t, g, grid, EXP_PEAK_CM, EXP_FWHM_CM)
            y = lineshape(t, g, grid, EXP_PEAK_CM, ss)
            m = band_metrics(grid, y)
            col, ls = styles[(tag, quantum)]
            lab = (f"{tag}, quantum  (blue/red {m['asymmetry']:.2f})"
                   if quantum else
                   f"{tag}, classical  (blue/red {m['asymmetry']:.2f})")
            # peak-aligned: the residual origin shift is reported separately as
            # the Stokes-shift test, and is not a band-SHAPE statement
            ax[1].plot(grid - grid[y.argmax()], y / y.max(), color=col, ls=ls,
                       lw=1.4, label=lab)
    ax[1].axvline(0.0, color="0.7", lw=0.8, zorder=0)
    for off in (800, 1600, 2000):
        ax[1].axvline(off, color="0.88", lw=0.7, zorder=0)
    ax[1].set_xlim(-1800, 2600)
    ax[1].set_ylim(0, 1.12)
    ax[1].set_xlabel("wavenumber relative to the band maximum (cm$^{-1}$)")
    ax[1].set_ylabel("normalised absorption")
    ax[1].set_title("width matched; asymmetry and vibronic wing free",
                    fontsize=10)
    ax[1].annotate("vibronic wing:\nquantum tracks it,\nclassical is 8x low",
                   xy=(2000, 0.15), xytext=(1450, 0.62), fontsize=7.5,
                   color="0.25",
                   arrowprops=dict(arrowstyle="->", color="0.45", lw=0.8))
    ax[1].legend(frameon=False, fontsize=7.5, loc="upper right")

    for a_ in ax:
        a_.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = ROOT / "results" / "quantum_cumulant_band.png"
    fig.savefig(out, dpi=190)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

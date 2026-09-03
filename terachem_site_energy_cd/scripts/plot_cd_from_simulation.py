#!/usr/bin/env python3
"""Draw the forward-predicted CD couplet against the measured one.

Consumes results/cd_from_simulation.npz written by cd_from_simulation.py.
Panel (a) is the absolute molar-ellipticity comparison; panel (b) is the
dissymmetry g = Delta_eps/eps, which removes the transition-dipole magnitude
and is the comparison that survives the |mu| systematic.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
NM = lambda nu: 1e7 / nu          # noqa: E731


def main() -> None:
    d = np.load(RES / "cd_from_simulation.npz")
    meta = json.loads((RES / "cd_from_simulation.json").read_text())
    grid = d["grid_cm"]
    lam = NM(grid)
    win = (lam > 460) & (lam < 580)

    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.8))

    a = ax[0]
    a.axhline(0.0, lw=0.6, color="0.6")
    a.plot(lam[win], d["experiment_theta"][win] / 1e5, color="k", lw=2.0,
           label="measured couplet (Nguyen Table S3)")
    a.plot(lam[win], d["primary_slow_theta"][win] / 1e5, color="C3", lw=1.8,
           label="predicted, slow disorder (primary)")
    a.plot(lam[win], d["raw_detuning_theta"][win] / 1e5, color="C0", lw=1.2,
           ls="--", label="predicted, raw detuning (double-counted)")
    a.plot(lam[win], d["degenerate_theta"][win] / 1e5, color="0.55", lw=1.0,
           ls=":", label=r"predicted, $\Delta=0$")
    a.set_xlabel("wavelength (nm)")
    a.set_ylabel(r"$\Delta[\theta]$  ($10^{5}$ deg cm$^2$ dmol$^{-1}$)")
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title("(a) absolute couplet, no fitted parameters", fontsize=9)

    b = ax[1]
    b.axhline(0.0, lw=0.6, color="0.6")
    g_pred = (d["primary_slow_theta"] / 3298.2) / d["primary_slow_eps"]
    g_exp = (d["experiment_theta"] / 3298.2) / meta["cases"]["experiment"]["eps_max"]
    core = (lam > 495) & (lam < 545)
    b.plot(lam[core], 1e4 * g_exp[core], color="k", lw=2.0,
           label=r"measured (constant $\varepsilon_{\max}$)")
    b.plot(lam[core], 1e4 * g_pred[core], color="C3", lw=1.8, label="predicted")
    b.set_xlabel("wavelength (nm)")
    b.set_ylabel(r"dissymmetry $g\times10^{4}$")
    b.legend(fontsize=8, frameon=False)
    b.set_title(r"(b) $g=\Delta\varepsilon/\varepsilon$ near the band centre",
                fontsize=9)
    b.text(0.02, 0.03,
           "measured trace divides by the peak $\\varepsilon$ only;\n"
           "a measured absorption band is needed in the wings",
           transform=b.transAxes, fontsize=6.5, color="0.35")

    p = meta["cases"]["primary_slow"]
    e = meta["cases"]["experiment"]
    fig.suptitle(
        f"predicted separation {p['separation_cm']:.0f} cm$^{{-1}}$ vs measured "
        f"{e['separation_cm']:.0f}, with $2|J| = 54$ cm$^{{-1}}$",
        fontsize=9, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"cd_from_simulation.{ext}", dpi=200,
                    bbox_inches="tight")
    print("wrote", RES / "cd_from_simulation.pdf")


if __name__ == "__main__":
    main()

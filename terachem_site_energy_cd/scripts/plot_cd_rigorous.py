#!/usr/bin/env python3
"""Plot the rigorous forward CD against the measured couplet.

Consumes results/cd_rigorous.{npz,json}. Every lineshape parameter in the
prediction is measured; nothing is fitted to the CD spectrum being compared.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path(__file__).resolve().parent.parent / "results"


def main() -> None:
    d = np.load(RES / "cd_rigorous.npz")
    meta = json.loads((RES / "cd_rigorous.json").read_text())
    grid = d["grid_cm"]
    lam = 1e7 / grid
    win = (lam > 455) & (lam < 590)

    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.7))

    a = ax[0]
    a.axhline(0.0, lw=0.6, color="0.6")
    a.plot(lam[win], d["experiment_theta"][win] / 1e5, color="k", lw=2.0,
           label="measured (Nguyen Table S3)")
    a.plot(lam[win], d["rigorous_theta"][win] / 1e5, color="C3", lw=1.9,
           label="predicted")
    a.plot(lam[win], d["degenerate_control_theta"][win] / 1e5, color="0.55",
           lw=1.0, ls=":", label=r"control: no static disorder")
    a.set_xlabel("wavelength (nm)")
    a.set_ylabel(r"$\Delta[\theta]$  ($10^{5}$ deg cm$^2$ dmol$^{-1}$)")
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title("(a) absolute couplet", fontsize=9)

    b = ax[1]
    b.plot(lam[win], d["rigorous_eps"][win] / 1e5, color="C3", lw=1.9,
           label="predicted absorption")
    b.set_xlabel("wavelength (nm)")
    b.set_ylabel(r"$\varepsilon$  ($10^{5}$ M$^{-1}$cm$^{-1}$)")
    b.axhline(meta["cases"]["experiment"]["eps_max"] / 1e5, color="k", lw=1.2,
              ls="--", label=r"measured $\varepsilon_{\max}$")
    b.legend(fontsize=7, frameon=False)
    b.set_title("(b) absorption, same model", fontsize=9)

    c = ax[2]
    dd = meta["disorder"]
    labels, vals, colors = [], [], []
    for band in dd["bands"]:
        if band["treatment"] == "Franck-Condon":
            labels.append(f"{band['omega_cm']:.0f}")
            vals.append(band["huang_rhys"])
            colors.append("C0")
    c.bar(labels, vals, color=colors)
    c.set_xlabel(r"mode $\omega_i$ (cm$^{-1}$), from the measured spectral density")
    c.set_ylabel(r"Huang--Rhys $S_i=\lambda_i/\omega_i$")
    c.set_title("(c) vibronic modes, not fitted", fontsize=9)
    c.text(0.97, 0.93,
           f"static disorder\n$\\sigma_{{\\rm site}}$ = "
           f"{dd['sigma_site_static_cm']:.0f} cm$^{{-1}}$\n"
           f"homogeneous width: none",
           transform=c.transAxes, ha="right", va="top", fontsize=7.5)

    r = meta["cases"]["rigorous"]
    e = meta["cases"]["experiment"]
    fig.suptitle(
        f"predicted separation {r['separation_cm']:.0f} cm$^{{-1}}$ vs measured "
        f"{e['separation_cm']:.0f}   |   $2|J|$ = 54 cm$^{{-1}}$   |   "
        f"amplitude {r['peak_to_peak']/e['peak_to_peak']:.2f}$\\times$, "
        f"$g_{{\\rm abs}}$ {r['g_abs_peak']/e['g_abs_peak']:.2f}$\\times$",
        fontsize=9, y=1.03)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"cd_rigorous.{ext}", dpi=200, bbox_inches="tight")
    print("wrote", RES / "cd_rigorous.pdf")


if __name__ == "__main__":
    main()

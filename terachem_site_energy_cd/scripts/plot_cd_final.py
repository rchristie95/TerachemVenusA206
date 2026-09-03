#!/usr/bin/env python3
"""Final CD figure: the couplet, the absorption, and the width-free moment."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path(__file__).resolve().parent.parent / "results"
C = 19322.0
UNIT = 1e5


def main() -> None:
    d = np.load(RES / "cd_final.npz"); g = d["grid_cm"]; lam = 1e7 / g
    L = lambda a, c, w: a / (1 + ((g - c) / w) ** 2)
    obs = sum(L(-a * UNIT, c, w) for a, c, w in
              ((3.59, C - 130.79, 293.51), (-2.27, C + 130.79, 439.25),
               (-0.34, 20774.0, 565.79)))
    mom = lambda y: np.trapezoid(y * (g - C) / g, g)
    top = np.load(RES / "digitised" / "venus_abs_top_digitised.npz")
    win = (lam > 440) & (lam < 600)

    fig, ax = plt.subplots(1, 3, figsize=(13.4, 3.9))

    a = ax[0]
    a.axhline(0, lw=.6, color="0.6", ls=":")
    a.plot(lam[win], -obs[win] / UNIT, color="k", lw=2.0, label="measured (Table S3)")
    a.plot(lam[win], -d["mu_exp_theta"][win] / UNIT, color="C3", lw=1.9, ls="--",
           label=r"this work, $|\mu|$ = 7.74 D (measured)")
    a.plot(lam[win], -d["abinitio_theta"][win] / UNIT, color="C0", lw=1.1, ls=":",
           label=r"this work, $|\mu|$ = 9.40 D (TDDFT)")
    a.set_xlabel("Wavelength (nm)"); a.set_xlim(440, 600)
    a.set_ylabel(r"$\Delta$-[$\Theta$] $\times10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title("(a) difference CD", fontsize=9)

    b = ax[1]
    e = d["mu_exp_eps"]; e0 = d["abinitio_eps"]
    b.plot(top["wavelength_nm"], top["normalised"] * 1.844, color="k", lw=2.0,
           label="measured band top")
    b.plot(lam, e / 1e5, color="C3", lw=1.9, ls="--", label=r"$|\mu|$ = 7.74 D")
    b.plot(lam, e0 / 1e5, color="C0", lw=1.1, ls=":", label=r"$|\mu|$ = 9.40 D")
    b.axhline(1.844, color="0.5", lw=.8, ls="-.", label=r"measured $\varepsilon_{\max}$")
    b.set_xlim(470, 570); b.set_ylim(0, 3.2)
    b.set_xlabel("Wavelength (nm)")
    b.set_ylabel(r"$\varepsilon$ ($10^{5}$ M$^{-1}$cm$^{-1}$)")
    b.legend(fontsize=7, frameon=False); b.set_title("(b) absorption", fontsize=9)

    c = ax[2]
    mo = mom(obs)
    vals = [mom(d["abinitio_theta"]) / mo, mom(d["mu_exp_theta"]) / mo, 1.0]
    names = [r"$|\mu|$ 9.40 D", r"$|\mu|$ 7.74 D", "measured"]
    c.bar(names, vals, color=["C0", "C3", "k"])
    c.axhline(1.0, color="0.4", lw=.8, ls="--")
    c.set_ylabel(r"$M_1$ / $M_1^{\rm meas}$")
    c.set_title(r"(c) first moment: width- and $\Delta$-free", fontsize=9)
    for i, v in enumerate(vals):
        c.text(i, v + 0.05, f"{v:.2f}", ha="center", fontsize=9)
    c.set_ylim(0, 2.7)
    c.text(0.5, 0.62, r"$M_1=-P\pi\langle\bar\nu J T\rangle$" "\n"
           "exact for any detuning\nand any lineshape",
           transform=c.transAxes, ha="center", fontsize=7.5, color="0.3")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"Fig_CD_Final.{ext}", dpi=200, bbox_inches="tight")
    print(f"  M1 ratios: abinitio {vals[0]:.2f}, mu_exp {vals[1]:.2f}")
    print(f"  wrote {RES/'Fig_CD_Final.pdf'}")


if __name__ == "__main__":
    main()

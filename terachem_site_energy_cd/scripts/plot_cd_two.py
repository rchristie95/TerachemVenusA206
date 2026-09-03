#!/usr/bin/env python3
"""The two predicted CD spectra, and the test for a protonated population.

(a) both predicted couplets against the measured one, over the full CD window;
(b) the measured absorption over 350-600 nm. A protonated (neutral) chromophore
    absorbs near 400 nm; that window is flat, which bounds the neutral fraction.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path(__file__).resolve().parent.parent / "results"
C, UNIT = 19322.0, 1e5


def main() -> None:
    d = np.load(RES / "cd_final.npz"); g = d["grid_cm"]; lam = 1e7 / g
    L = lambda A, c, w: A / (1 + ((g - c) / w) ** 2)
    obs = sum(L(-A * UNIT, c, w) for A, c, w in
              ((3.59, C - 130.79, 293.51), (-2.27, C + 130.79, 439.25),
               (-0.34, 20774.0, 565.79)))
    full = np.load(RES / "digitised" / "venus_absorption_digitised.npz")

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.0))
    w = (lam > 400) & (lam < 600)
    a = ax[0]
    a.axhline(0, lw=.6, color="0.6", ls=":")
    a.plot(lam[w], -obs[w] / UNIT, color="k", lw=2.2, label="measured (Table S3)")
    a.plot(lam[w], -d["mu_exp_theta"][w] / UNIT, color="C3", lw=2.0, ls="--",
           label=r"predicted, $|\mu|$ = 7.74 D (measured absorption)")
    a.plot(lam[w], -d["abinitio_theta"][w] / UNIT, color="C0", lw=1.4, ls=":",
           label=r"predicted, $|\mu|$ = 9.40 D (TDDFT)")
    a.set_xlim(400, 600); a.set_xlabel("Wavelength (nm)")
    a.set_ylabel(r"$\Delta$-[$\Theta$] $\times10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    a.legend(fontsize=7.5, frameon=False, loc="lower right")
    a.set_title("(a) the two predicted couplets", fontsize=9)

    b = ax[1]
    lm, yv = full["wavelength_nm"], full["normalised"]
    b.plot(lm, yv, color="k", lw=1.8, label="dVenus-TD (Fig. 2e)")
    b.axvspan(380, 420, color="C1", alpha=.20,
              label="protonated form would absorb here")
    b.set_xlim(350, 600); b.set_ylim(-0.03, 1.05)
    b.set_xlabel("Wavelength (nm)"); b.set_ylabel("normalised absorbance")
    b.legend(fontsize=7.5, frameon=False, loc="upper left")
    b.set_title("(b) no neutral-form band: 380-420 nm is flat", fontsize=9)
    b.text(0.97, 0.55, "380-420 nm mean\n" r"$-0.06\%$ of the peak" "\n"
           r"(max $+0.3\%$)", transform=b.transAxes, ha="right", fontsize=8)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"Fig_CD_TwoSpectra.{ext}", dpi=200, bbox_inches="tight")
    print("wrote", RES / "Fig_CD_TwoSpectra.pdf")


if __name__ == "__main__":
    main()

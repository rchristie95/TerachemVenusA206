#!/usr/bin/env python3
"""Round-4 CD figure: the prediction, and the published analysis run on it.

Panel (b) uses no experimental input at all -- it fits the simulated spectrum
with the published model and reports what that model recovers, on data where the
coupling is known exactly. That is the claim round 4 leads with.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path(__file__).resolve().parent.parent / "results"
G1, UNIT, JTRUE = 19322.0, 1e5, 27.06
L = lambda x, A, c, w: A / (1 + ((x - c) / w) ** 2)


def triple(x, a1, w1, a2, w2, dl, a3, c3, w3):
    return L(x, a1, G1 - dl, w1) + L(x, a2, G1 + dl, w2) + L(x, a3, c3, w3)


def main() -> None:
    d = np.load(RES / "cd_final.npz"); g = d["grid_cm"]; lam = 1e7 / g
    th = -d["mu_exp_theta"]; th0 = -d["abinitio_theta"]        # published sign
    fitw = (g > 17800) & (g < 21800)
    x, y = g[fitw], th[fitw] / np.abs(th[fitw]).max()
    best = None
    for d0 in (60., 130., 232., 350., 500.):
        for sgn in (1., -1.):
            try:
                pr, _ = curve_fit(triple, x, y,
                    p0=[sgn, 293.5, -.63*sgn, 439.3, d0, -.09*sgn, 20774., 565.8],
                    bounds=([-20,40,-20,40,5,-20,19800,40],
                            [20,3000,20,3000,3000,20,22000,3000]), maxfev=400000)
            except Exception:
                continue
            r = 1 - np.sum((y - triple(x, *pr))**2)/np.sum((y - y.mean())**2)
            if best is None or r > best[1]:
                best = (pr, r)
    pr, r2 = best
    sc = np.abs(th[fitw]).max()

    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.0))
    w = (lam > 430) & (lam < 600)
    a = ax[0]
    a.axhline(0, lw=.6, color="0.6", ls=":")
    a.plot(lam[w], th0[w]/UNIT, color="C0", lw=1.9, ls="-",
           label=r"$|\mu| = 9.40$ D (computed here)")
    a.plot(lam[w], th[w]/UNIT, color="C3", lw=1.9, ls="--",
           label=r"$|\mu| = 7.74$ D (experimental moment)")
    a.set_xlim(430, 600); a.set_xlabel("Wavelength (nm)")
    a.set_ylabel(r"$\Delta$-[$\Theta$] $\times10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    a.legend(fontsize=8, frameon=False, loc="lower right")
    a.set_title(r"(a) predicted couplet; amplitude $\propto|\mu|^{4}$", fontsize=9)

    b = ax[1]
    b.axhline(0, lw=.6, color="0.6", ls=":")
    b.plot(lam[w], th[w]/UNIT, color="k", lw=2.0, label="predicted spectrum")
    b.plot(lam[fitw], triple(x, *pr)*sc/UNIT, color="C3", lw=1.8, ls="--",
           label="published triple-Lorentzian fit")
    for i, (A, c, ww) in enumerate(((pr[0], G1-pr[4], pr[1]),
                                    (pr[2], G1+pr[4], pr[3]),
                                    (pr[5], pr[6], pr[7]))):
        b.plot(lam[fitw], L(x, A, c, ww)*sc/UNIT, lw=.9, color="C0", alpha=.7,
               label="fitted components" if i == 0 else None)
    b.set_xlim(430, 600); b.set_ylim(-7.5, 7.5)
    b.set_xlabel("Wavelength (nm)")
    b.set_ylabel(r"$\Delta$-[$\Theta$] $\times10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    b.legend(fontsize=8, frameon=False, loc="upper left")
    b.set_title("(b) the published analysis, applied to it", fontsize=9)
    b.text(.97, .06, f"apparent $U$ = {pr[4]:.0f} cm$^{{-1}}$\n"
           f"true $|J|$ = {JTRUE:.1f} cm$^{{-1}}$\n"
           f"overstated {pr[4]/JTRUE:.1f}$\\times$   ($R^2$ = {r2:.3f})",
           transform=b.transAxes, ha="right", va="bottom", fontsize=8.5)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"Fig_CD_Round4.{ext}", dpi=200, bbox_inches="tight")
    print(f"  apparent U = {pr[4]:.1f} cm^-1 from |J| = {JTRUE} "
          f"({pr[4]/JTRUE:.1f}x), R2 = {r2:.4f}")
    (RES / "cd_round4_fit.json").write_text(json.dumps(
        dict(apparent_U_cm=float(pr[4]), true_J_cm=JTRUE, r2=float(r2)), indent=2))
    print("  wrote", RES / "Fig_CD_Round4.pdf")


if __name__ == "__main__":
    main()

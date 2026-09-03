#!/usr/bin/env python3
"""Three-panel CD figure in the plotting conventions of Nguyen et al. 2025.

Conventions taken from that paper and matched here:
  * Delta-[Theta] in units of 1e5 deg cm^2 dmole^-1 (their Fig. 4b axis reads
    "x10^-5"; Table S3's amplitudes equal the plotted extrema, which fixes the
    scale) and in THEIR sign, Delta[theta] = TDX - TD, so the negative lobe lies
    to the blue;
  * measured trace black solid, model red dashed, as in their Fig. 4b/d;
  * the normalised absorption plotted over the top 10% only, as in their Fig. 3b;
  * the constrained triple-Lorentzian of their Table S3 used as the analysis,
    with G1 fixed and the couplet forced symmetric about it at +/- delta, so the
    "apparent coupling strength" U = delta is read off our spectrum exactly as
    it was read off theirs.

Inputs: results/cd_rigorous.npz (the forward model), results/digitised/*.npz
(their Figs. 2e and 3b, digitised).  Writes results/Fig_CD_Comparison.{pdf,png}.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path(__file__).resolve().parent.parent / "results"
G1 = 19322.0                       # TDX absorption Gaussian, Nguyen Table S3
DIFF_SIGN = -1.0                   # our spectra are the TD couplet; theirs TDX-TD
UNIT = 1e5
EXP_S3 = ((3.59, G1 - 130.79, 293.51),
          (-2.27, G1 + 130.79, 439.25),
          (-0.34, 20774.0, 565.79))


def lorentz(x, a, c, w):
    return a / (1.0 + ((x - c) / w) ** 2)


def constrained_triple(x, a1, w1, a2, w2, delta, a3, c3, w3):
    """Nguyen Table S3 model: the couplet is pinned symmetrically about G1."""
    return (lorentz(x, a1, G1 - delta, w1)
            + lorentz(x, a2, G1 + delta, w2)
            + lorentz(x, a3, c3, w3))


def fit_apparent_coupling(nu, y):
    """Their analysis, run on our spectrum. Returns (params, U in cm^-1, R^2)."""
    scale = np.abs(y).max()
    yn = y / scale
    best = None
    for d0 in (60., 130., 200., 300., 450., 650.):
        for s in (1.0, -1.0):
            try:
                p, _ = curve_fit(
                    constrained_triple, nu, yn,
                    p0=[s, 293.51, -0.63 * s, 439.25, d0, -0.09 * s, 20774., 565.79],
                    bounds=([-20, 40, -20, 40, 5, -20, 19800, 40],
                            [20, 3000, 20, 3000, 3000, 20, 22000, 3000]),
                    maxfev=600000)
            except Exception:
                continue
            r = yn - constrained_triple(nu, *p)
            r2 = 1.0 - r.dot(r) / ((yn - yn.mean()) ** 2).sum()
            if best is None or r2 > best[2]:
                best = (p * np.array([scale, 1, scale, 1, 1, scale, 1, 1]), p, r2)
    return best[0], float(best[0][4]), best[2]


def main() -> None:
    d = np.load(RES / "cd_rigorous.npz")
    meta = json.loads((RES / "cd_rigorous.json").read_text())
    nu = d["grid_cm"]
    lam = 1e7 / nu

    # our spectra, put into Nguyen's sign convention
    theta = DIFF_SIGN * d["rigorous_theta"]
    theta_ctl = DIFF_SIGN * d["degenerate_control_theta"]
    eps = d["rigorous_eps"]
    exp_theta = sum(lorentz(nu, a * UNIT, c, w) for a, c, w in EXP_S3)

    top = np.load(RES / "digitised" / "venus_abs_top_digitised.npz")
    full = np.load(RES / "digitised" / "venus_absorption_digitised.npz")

    win = (lam > 400) & (lam < 600)
    fig, ax = plt.subplots(1, 3, figsize=(13.2, 3.9))

    # ---- (a) absorption, Fig. 3b convention ------------------------------
    a = ax[0]
    a.plot(top["wavelength_nm"], top["normalised"], color="k", lw=2.0,
           label="dVenus-TD (Nguyen Fig. 3b)")
    e_n = eps / eps.max()
    m = e_n >= 0.895
    a.plot(lam[m], e_n[m], color="crimson", lw=2.0, ls="--", label="this work")
    a.set_ylim(0.90, 1.005)
    a.set_xlim(508, 526)
    a.set_xlabel("Wavelength (nm)")
    a.set_ylabel("Normalized Absorbance")
    a.legend(fontsize=7, frameon=False, loc="lower right")
    a.set_title(r"(a) band maximum: width at 0.90", fontsize=9)
    a.text(0.03, 0.97, "width at 0.90 of peak\n"
           r"measured  252 cm$^{-1}$" "\n" r"this work  374 cm$^{-1}$",
           transform=a.transAxes, fontsize=7, va="top")

    # ---- (b) Delta-[Theta], Fig. 4d convention ---------------------------
    b = ax[1]
    b.axhline(0.0, lw=0.6, color="0.6", ls=":")
    b.plot(lam[win], exp_theta[win] / UNIT, color="k", lw=2.0,
           label="dVenus-TD (Table S3)")
    b.plot(lam[win], theta[win] / UNIT, color="crimson", lw=2.0, ls="--",
           label="this work")
    b.plot(lam[win], theta_ctl[win] / UNIT, color="0.6", lw=1.0, ls=":",
           label=r"no static disorder")
    b.set_xlim(400, 600)
    b.set_xlabel("Wavelength (nm)")
    b.set_ylabel(r"$\Delta$-[$\Theta$] $\times 10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    b.legend(fontsize=7, frameon=False, loc="upper left")
    b.set_title(r"(b) difference CD, no fitted parameters", fontsize=9)

    # ---- (c) their triple-Lorentzian analysis, run on our spectrum -------
    fitwin = (nu > 17800) & (nu < 21800)
    p, u_pred, r2 = fit_apparent_coupling(nu[fitwin], theta[fitwin])
    c = ax[2]
    c.axhline(0.0, lw=0.6, color="0.6", ls=":")
    c.plot(lam[win], theta[win] / UNIT, color="k", lw=2.0, label="this work")
    c.plot(lam[fitwin], constrained_triple(nu[fitwin], *p) / UNIT,
           color="crimson", lw=1.8, ls="--", label="Triple Lorentzian Fit")
    for i, (amp, ctr, wid) in enumerate(((p[0], G1 - p[4], p[1]),
                                         (p[2], G1 + p[4], p[3]),
                                         (p[5], p[6], p[7]))):
        c.plot(lam[fitwin], lorentz(nu[fitwin], amp, ctr, wid) / UNIT,
               lw=0.9, color="C0", alpha=0.7,
               label="fitted components" if i == 0 else None)
    c.set_ylim(-7.0, 7.0)
    c.set_xlim(400, 600)
    c.set_xlabel("Wavelength (nm)")
    c.set_ylabel(r"$\Delta$-[$\Theta$] $\times 10^{-5}$ deg cm$^2$ dmole$^{-1}$")
    c.legend(fontsize=7, frameon=False, loc="upper left")
    c.set_title(r"(c) their analysis, applied to our spectrum", fontsize=9)
    c.text(0.97, 0.06,
           f"apparent $U$ = {u_pred:.0f} cm$^{{-1}}$\n"
           f"true $|J|$ = {27.06:.1f} cm$^{{-1}}$\n$R^2$ = {r2:.4f}",
           transform=c.transAxes, ha="right", va="bottom", fontsize=7.5)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(RES / f"Fig_CD_Comparison.{ext}", dpi=200, bbox_inches="tight")

    print(f"  apparent coupling from OUR spectrum, their constrained fit: "
          f"U = {u_pred:.1f} cm^-1  (R^2 = {r2:.4f})")
    print(f"  their value from the measured spectrum:                     "
          f"U = 130.8 cm^-1")
    print(f"  the true coupling in the Hamiltonian that generated ours:   "
          f"|J| = 27.06 cm^-1, 2|J| = 54.1")
    (RES / "cd_apparent_coupling.json").write_text(json.dumps(dict(
        apparent_U_cm=u_pred, r2=r2, true_J_cm=27.06,
        params=dict(zip("a1 w1 a2 w2 delta a3 c3 w3".split(), p.tolist()))),
        indent=2))
    print(f"  wrote {RES/'Fig_CD_Comparison.pdf'}")


if __name__ == "__main__":
    main()

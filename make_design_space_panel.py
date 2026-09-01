#!/usr/bin/env python3
"""Figure 5(d): the coupling/detuning design space for a beta-barrel FP dimer.

Answers "could any fluorescent-protein dimer be engineered to delocalise?" by
putting the measured 500-frame ensemble in the (|Delta|, 2|J|) plane together
with the ceiling each engineering lever can reach.

Nothing here is a fit. The three levers are read off the ensemble itself:

  orientation  kappa is recovered from the stored point-dipole coupling,
               kappa = J_PDA R^3 / (mu_A mu_B), and its ceiling is the
               head-to-tail |kappa| = 2;
  separation   J ~ R^-3 (justified because <J_PDA/J_TDC> = 0.92 here, i.e. the
               far field already holds), and the floor on R is barrel contact,
               so R_min is the barrel diameter quoted in the text;
  dipole       |mu|^2, whose ceiling means abandoning the FP chromophore, so it
               is drawn as an out-of-class extrapolation rather than a lever.

The delocalisation criteria are the two used elsewhere in the manuscript:
|Delta| < 2|J| (a true 50/50 window) and the eigenstate mixing
sin 2theta = 2|J|/Omega > 0.5.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MANIFESTS = [REPO / "tddft_ensemble_decoherence/manifests" / f"v5camp_{w}_frames.csv"
             for w in "abc"]
EV = 8065.544005
HARTREE = 219474.6313632
BOHR = 0.529177210903
DEBYE_PER_AU = 2.541746473


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--barrel-diameter-A", type=float, default=22.0,
                   help="centroid floor at barrel contact; the manuscript quotes >=22 A")
    p.add_argument("--cyanine-dipole-D", type=float, default=14.0,
                   help="out-of-class transition dipole for the dashed extrapolation")
    p.add_argument("--out-dir", type=Path, default=HERE)
    return p.parse_args()


def load() -> dict[str, np.ndarray]:
    rows = [r for m in MANIFESTS for r in csv.DictReader(m.open())
            if r["convergence_status"] == "complete"]
    g = lambda k: np.array([float(r[k]) for r in rows])  # noqa: E731
    ea, eb = g("site_A_energy_eV"), g("site_B_energy_eV")
    return {
        "delta": (ea - eb) * EV,
        "J": g("J_screened_cm-1"),
        "J_pda_vac": g("J_pda_vacuum_cm-1"),
        "J_pda_scr": g("J_pda_screened_cm-1"),
        "R": g("separation_A"),
        "muA": np.sqrt(1.5 * g("site_A_oscillator_strength") / (ea / 27.211386)),
        "muB": np.sqrt(1.5 * g("site_B_oscillator_strength") / (eb / 27.211386)),
    }


def main() -> None:
    args = parse_args()
    d = load()
    absD, twoJ = np.abs(d["delta"]), 2.0 * d["J"]
    mixing = lambda t: t / np.hypot(d["delta"], t)          # noqa: E731  sin 2theta
    frac_mix = lambda t: float((mixing(t) > 0.5).mean())    # noqa: E731
    frac_win = lambda t: float((absD < t).mean())           # noqa: E731

    kappa = np.abs(d["J_pda_vac"] / HARTREE * (d["R"] / BOHR) ** 3 / (d["muA"] * d["muB"]))
    far_field = float(np.mean(d["J_pda_scr"] / d["J"]))
    R_mean, k_mean = float(d["R"].mean()), float(kappa.mean())
    mu_D = float(d["muA"].mean() * DEBYE_PER_AU)

    lev_kappa = 2.0 / k_mean
    lev_R = (R_mean / args.barrel_diameter_A) ** 3
    lev_mu = (args.cyanine_dipole_D / mu_D) ** 2
    base = float(twoJ.mean())
    barrel = base * lev_kappa * lev_R
    out_of_class = barrel * lev_mu
    need = 0.577 * float(np.median(absD))                   # mixing 0.5 at the median frame

    print(f"  far-field check <J_PDA/J_TDC> = {far_field:.3f}")
    print(f"  kappa = {k_mean:.3f} (ceiling 2, headroom {lev_kappa:.2f}x)")
    print(f"  R = {R_mean:.2f} A -> {args.barrel_diameter_A:.0f} A, headroom {lev_R:.2f}x")
    print(f"  |mu| = {mu_D:.2f} D -> {args.cyanine_dipole_D:.0f} D, headroom {lev_mu:.2f}x")
    for lab, t in (("as built", base), ("optimal barrel", barrel), ("out of class", out_of_class)):
        print(f"  {lab:15s} 2|J| = {t:6.0f}  |D|<2|J| {100*frac_win(t):5.1f}%  "
              f"mixing>0.5 {100*frac_mix(t):5.1f}%")
    print(f"  median frame needs 2|J| > {need:.0f}; optimal barrel short by {need/barrel:.1f}x")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.formatter.use_mathtext": True,
        "axes.labelsize": 12, "legend.fontsize": 8.5,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
    })
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    xlo, xhi, ylo, yhi = 1.0, 3.0e3, 42.0, 4.2e2
    x = np.logspace(np.log10(xlo), np.log10(xhi), 400)
    med = float(np.median(absD))

    ax.fill_between(x, np.minimum(0.577 * x, yhi), yhi,
                    color="#1f6aa5", alpha=0.06, lw=0)
    ax.plot(x, x, color="#b02020", lw=1.2, ls="--", zorder=3,
            label=r"$|\Delta| = 2|J|$")
    ax.plot(x, 0.577 * x, color="#1f6aa5", lw=1.2, ls="-.", zorder=3,
            label=r"mixing $\sin 2\theta = 0.5$")
    ax.scatter(np.clip(absD, xlo, None), twoJ, s=6, color="0.4", alpha=0.4,
               lw=0, zorder=2, label=f"{absD.size} QM/MM frames")

    # the ladder at the median detuning
    ax.plot([med, med], [ylo, need], color="0.35", lw=0.9, ls=":", zorder=4)
    for t, lab, col in ((base, "as built", "#111111"),
                        (barrel, "optimal barrel", "#1f6aa5"),
                        (out_of_class, "cyanine", "0.5")):
        ax.axhline(t, color=col, lw=1.1, zorder=4,
                   ls="-" if col != "0.5" else ":")
        ax.text(xhi * 0.96, t * 1.05, f"{lab}: {t:.0f}", ha="right", va="bottom",
                fontsize=8.5, color=col, zorder=6)
    ax.annotate("", xy=(med, barrel), xytext=(med, base), zorder=5,
                arrowprops=dict(arrowstyle="-|>", color="#1f6aa5", lw=1.4))
    ax.text(med * 1.12, np.sqrt(base * barrel), r"$1.9\times$", fontsize=9.5,
            color="#1f6aa5", va="center", zorder=6)
    ax.annotate("", xy=(med, need), xytext=(med, barrel), zorder=5,
                arrowprops=dict(arrowstyle="-|>", color="#b02020", lw=1.4))
    ax.text(med * 1.12, np.sqrt(barrel * need), r"$2.8\times$ short", fontsize=9.5,
            color="#b02020", va="center", zorder=6)

    ax.plot([med], [base], marker="o", ms=6.5, mfc="#e8a33d", mec="black",
            mew=0.9, zorder=7)
    ax.plot([med], [need], marker="*", ms=12, mfc="#b02020", mec="black",
            mew=0.7, zorder=7)
    ax.text(med * 0.50, need, r"required, $%.0f$" % need, fontsize=8.5,
            color="#b02020", ha="right", va="center", zorder=6)
    ax.text(8.0, 145.0, "delocalised", fontsize=10, color="#1f6aa5",
            ha="center", va="center", zorder=6)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
    ax.set_xlabel(r"$|\Delta|$ (cm$^{-1}$)")
    ax.set_ylabel(r"$2|J|$ (cm$^{-1}$)")
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([50, 100, 200, 400]))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["50", "100", "200", "400"]))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.grid(alpha=0.2, lw=0.6, which="major")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    for s_ in ax.spines.values():
        s_.set_linewidth(0.8)
    ax.legend(loc="upper left", frameon=True, framealpha=0.96, edgecolor="0.75",
              fancybox=False, borderpad=0.4, handlelength=1.8)
    fig.tight_layout()
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(args.out_dir / f"Fig_Design_Space{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print(f"  wrote {args.out_dir/'Fig_Design_Space.pdf'}")


if __name__ == "__main__":
    main()

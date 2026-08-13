#!/usr/bin/env python3
"""Figure 7 panels and the Figure 6(b) site-energy histogram.

Figure 7 must show exactly the same plots as Figure 1, so it calls the same
generator, open_quantum_dynamics.regenerate_base_figures, changing only the two
site energies. The two frames are the extreme members of the QM/MM site-energy
ensemble:

  frame 465   Delta =   +12.9 cm^-1   near-degenerate, maximal delocalisation
  frame 1005  Delta = -1825.0 cm^-1   maximal detuning, one chromophore favoured

Figure 6(b) is replaced by a single histogram of the TDDFT site-energy
differences over the QM/MM ensemble.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

REPO = Path("/home/robson/PetaChem")
sys.path.insert(0, str(REPO))
import open_quantum_dynamics as oqd  # noqa: E402

E0 = oqd.E1                      # common centre; only the difference matters
FRAMES = {"deloc": ("frame 465", 12.9), "loc": ("frame 1005", -1825.0)}
OUT = REPO / "manuscript"


def panels():
    for tag, (label, delta) in FRAMES.items():
        out = OUT / f"frame_{tag}"
        out.mkdir(exist_ok=True)
        p = oqd.make_params(e1=E0 + delta / 2.0, e2=E0 - delta / 2.0)
        oqd.regenerate_base_figures(p, out, oqd.TF, oqd.DT, 20260618)
        om = np.hypot(delta, 2 * oqd.J_OPT)
        th = 0.5 * np.arctan2(2 * oqd.J_OPT, delta)
        print(f"  {label}: Delta={delta:8.1f}  Omega={om:7.1f}  "
              f"2J/Omega={2*oqd.J_OPT/om:.4f}  "
              f"minor={100*min(np.sin(th)**2, np.cos(th)**2):.2f}%  -> {out.name}/")


def histogram():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.load(REPO / "terachem_site_energy_cd/results/ensembles/ens_v2_all.npz")
    delta = d["e_a_cm"] - d["e_b_cm"]
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.hist(delta, bins=18, color="0.55", edgecolor="black", linewidth=0.7)
    ax.axvline(0.0, color="black", lw=1.0, ls=":")
    for x, c, lab in ((2 * oqd.J_OPT, "tab:red", r"$+2|J|$"),
                      (-2 * oqd.J_OPT, "tab:red", r"$-2|J|$")):
        ax.axvline(x, color=c, lw=1.4, ls="--")
    ax.set_xlabel(r"$\Delta = E_A - E_B$ (cm$^{-1}$)")
    ax.set_ylabel("frames")
    ax.text(0.03, 0.95,
            rf"$\langle|\Delta|\rangle$ = {np.abs(delta).mean():.0f} cm$^{{-1}}$"
            "\n"
            rf"median = {np.median(np.abs(delta)):.0f} cm$^{{-1}}$"
            "\n"
            rf"$2|J|$ = {2*oqd.J_OPT:.1f} cm$^{{-1}}$",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(OUT / "Fig_Detuning_Histogram.pdf")
    plt.close(fig)
    n = len(delta)
    print(f"  histogram: n={n}, <|D|>={np.abs(delta).mean():.0f}, "
          f"median={np.median(np.abs(delta)):.0f}, "
          f"frames with |D|<2|J|: {(np.abs(delta) < 2*oqd.J_OPT).sum()} "
          f"({100*(np.abs(delta) < 2*oqd.J_OPT).mean():.1f}%)")


if __name__ == "__main__":
    print("Figure 7 panels (same generator as Figure 1):")
    panels()
    print("Figure 6(b) replacement:")
    histogram()

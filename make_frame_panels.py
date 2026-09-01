#!/usr/bin/env python3
"""Figure 7 panels and the Figure 6(b) site-energy histogram.

Figure 7 must show exactly the same plots as Figure 1, so it calls the same
generator, open_quantum_dynamics.regenerate_base_figures, changing only the two
site energies. The two frames are the extreme members of the QM/MM site-energy
ensemble:

  466 ps   Delta =    -0.17 cm^-1   near-degenerate, maximal delocalisation
  714 ps   Delta = -2495.52 cm^-1   maximal detuning, one chromophore favoured

Figure 6(b) is a histogram of the TDDFT site-energy differences over the QM/MM
ensemble.

Both now come from the v5 campaign: 500 frames over a fresh 1 ns trajectory,
def2-SVP QM/MM with full-system 75,035-charge embedding, and transition
densities recomputed per frame and per site rather than one static density
rigidly fitted onto every frame.

Superseded on 2026-08-31. This script previously used frames 465 and 1005 of the
1000-frame rigid-fit ensemble (Delta = +12.9 and -1825.0), and the histogram
read terachem_site_energy_cd/results/ensembles/ens_v2_all.npz, which carries the
per-site 12 A embedding cutoff that Section "boundary" shows is worth
573-1205 cm^-1 rms. Each frame also carries its OWN recomputed J now, instead of
the ensemble mean.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

REPO = Path("/home/robson/PetaChem")
sys.path.insert(0, str(REPO))
import open_quantum_dynamics as oqd  # noqa: E402

E0 = oqd.E1                      # common centre; only the difference matters
# tag -> (label, Delta / cm^-1, J / cm^-1). J is the frame's own recomputed
# coupling, which the per-frame densities make available for the first time.
FRAMES = {"deloc": ("466 ps", -0.17, 24.55), "loc": ("714 ps", -2495.52, 24.65)}
OUT = REPO / "manuscript"
# The v5 manifests, merged. This is the same 500-frame ensemble behind Figure 5.
V5_MANIFESTS = [REPO / "tddft_ensemble_decoherence/manifests" / f"v5camp_{w}_frames.csv"
                for w in "abc"]


def panels():
    for tag, (label, delta, j) in FRAMES.items():
        out = OUT / f"frame_{tag}"
        out.mkdir(exist_ok=True)
        p = oqd.make_params(e1=E0 + delta / 2.0, e2=E0 - delta / 2.0, j_opt=j)
        oqd.regenerate_base_figures(p, out, oqd.TF, oqd.DT, 20260618)
        om = np.hypot(delta, 2 * j)
        th = 0.5 * np.arctan2(2 * j, delta)
        print(f"  {label}: Delta={delta:9.2f}  J={j:.2f}  Omega={om:7.1f}  "
              f"2J/Omega={2*j/om:.4f}  "
              f"minor={100*min(np.sin(th)**2, np.cos(th)**2):.2f}%  -> {out.name}/")


def histogram():
    """Figure 5(b). Styled to match Figure 5(c) by importing that panel's helpers,
    so the two histograms of the same 500-frame ensemble read as one figure."""
    import make_tandem_panels as mtp
    import csv
    plt = mtp.mpl()
    rows = [r for m in V5_MANIFESTS for r in csv.DictReader(m.open())
            if r["convergence_status"] == "complete"]
    delta = np.array([float(r["site_A_energy_cm-1"]) - float(r["site_B_energy_cm-1"])
                      for r in rows])
    two_j = 2 * oqd.J_OPT

    fig, ax = plt.subplots(figsize=mtp.FIGSIZE)
    counts, _, _ = ax.hist(delta, bins=mtp.PANEL_BINS, color=mtp.C_ABS, alpha=0.80,
                           edgecolor="white", lw=0.4)
    # No 2|J| marker here: this axis spans +-2500 cm^-1 while panel (c) spans
    # 15-32.5, so the scale separation is already visible from the two axes.
    # mean line plus +-sigma span, in the same order, colours and label format
    # that panel (c) uses for J
    mean_d, std_d = delta.mean(), delta.std(ddof=1)
    ax.axvline(mean_d, color=mtp.C_MEAN, lw=1.8,
               label=fr"$\bar\Delta = {mean_d:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.axvspan(mean_d - std_d, mean_d + std_d, color=mtp.C_MEAN, alpha=0.10,
               label=fr"$\pm\sigma = {std_d:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.set_xlabel(r"$\Delta = E_1 - E_2$ (cm$^{-1}$)")
    ax.set_ylabel("MD snapshots")
    mtp.panel_ylim(ax, counts, "Fig_Detuning_Histogram")
    mtp.style(ax)
    mtp.leg(ax, loc="upper left")
    fig.tight_layout()
    mtp.save_panel(fig, OUT / "Fig_Detuning_Histogram")
    plt.close(fig)

    n = len(delta)
    inside = int((np.abs(delta) < two_j).sum())
    print(f"  histogram: n={n}, <|D|>={np.abs(delta).mean():.0f}, "
          f"median={np.median(np.abs(delta)):.0f}, "
          f"frames with |D|<2|J|: {inside} ({100*inside/n:.1f}%)")


if __name__ == "__main__":
    print("Figure 7 panels (same generator as Figure 1):")
    panels()
    print("Figure 6(b) replacement:")
    histogram()

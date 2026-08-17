#!/usr/bin/env python3
"""Screening-convention reconciliation: ours vs Cusick, like for like.

THE PROBLEM THIS EXISTS TO FIX. Our headline J = 32.82 cm^-1 is computed with
epsilon = 1.77 and is therefore SCREENED. Cusick explicitly do not screen --
they disregard dielectric attenuation by the protein matrix on the grounds that
the system is only nanometres across. So the apparent bullseye (32.82 against
their 32-40) compares a screened number to an unscreened one. Unscreen ours and
it is ~58 cm^-1; screen theirs and it is ~18-23. Either way the agreement is not
what it looks like on the page, and quoting the coincidence would not survive a
referee who checks the conventions.

The honest -- and more interesting -- statement is that at this geometry the
methods agree to within the screening convention, and the choice of epsilon is
now the DOMINANT systematic, larger than the near-field correction originally
reported. That is a real methodological result and it is defensible.

A SECOND CORRECTION, independent of screening. The point-dipole partner of
J = 32.82 is NOT the 13.31 cm^-1 quoted in manuscript/Submit-JPCL-21April2026.tex.
The genuine partner is on disk: J_pda_cm over the same 1000 frames, mean
27.64 cm^-1, a ratio of 1.19 -- not 5.6. The 13.31 value belongs to a different
structure (crystal venus_dimer.pdb), a different density (a TDDFT cube) and a
single static frame; it pairs with the superseded 74.38, which reproduce_paper.py:115
already records as corrected to 20.83. Pairing 32.82 with 13.31 would restate a
2.5x-inflated near-field claim in new clothes.

Consequence for the near-field story: at R ~ 25 A the point-dipole approximation
is accurate to about 19%, so cheap point-dipole estimates are adequate for
fluorescent-protein dimers at biological separations. Unglamorous, but true, and
worth writing down.

Read-only. Writes results/screening_reconciliation.{md,json}.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PRODUCTION = ROOT / "coupling_nvt_production_cr2_1000_20260721"
OUT_DIR = Path(__file__).resolve().parent / "results"

EPSILON = 1.77

# Cusick et al. 2026. Their couplings are UNSCREENED by construction.
CUSICK = {
    "vdw": {
        "label": "dVenus-vdW (crystal contact)",
        "J_cm": (28.0, 38.0),
        "method": "point-dipole from structure",
        "R_A": 25.4,
        "delta_deg": (31.0, 31.0),
    },
    "tandem": {
        "label": "dVenus-TD (tandem, linkered)",
        "J_cm": (32.0, 40.0),
        "method": "from the Omega two-photon analysis",
        "R_A": 26.0,
        "delta_deg": (14.0, 20.0),
    },
}
# Table S2 of the same paper quotes a slightly wider span for the 1MYW dimer.
CUSICK_TABLE_S2_CM = (27.0, 43.0)

# The superseded pair, kept here only so the table can name and disown it.
LEGACY = {
    "J_tdc_cm": 74.38,
    "J_pda_cm": 13.31,
    "J_tdc_corrected_cm": 20.8285213577,   # reproduce_paper.py:115
    "structure": "venus_dimer.pdb (crystal)",
    "density": "TDDFT S2 cube",
    "frames": 1,
}


def load_production():
    rows = list(csv.DictReader(open(PRODUCTION / "coupling_samples.csv")))
    j = np.array([float(r["J_cm"]) for r in rows])
    pda = np.array([float(r["J_pda_cm"]) for r in rows])
    sep = np.array([float(r["separation_A"]) for r in rows])
    return j, pda, sep


def fmt_range(lo, hi, nd=1):
    return f"{lo:.{nd}f}" if lo == hi else f"{lo:.{nd}f}-{hi:.{nd}f}"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    j, pda, sep = load_production()
    n = len(j)

    ours = {
        "n_frames": n,
        "R_A_mean": float(sep.mean()),
        "R_A_std": float(sep.std()),
        "tdc_screened_cm": float(j.mean()),
        "tdc_screened_std_cm": float(j.std()),
        "tdc_unscreened_cm": float(j.mean() * EPSILON),
        "pda_screened_cm": float(pda.mean()),
        "pda_screened_std_cm": float(pda.std()),
        "pda_unscreened_cm": float(pda.mean() * EPSILON),
        "tdc_over_pda": float(j.mean() / pda.mean()),
        "sigma_over_mean": float(j.std() / j.mean()),
    }

    theirs = {}
    for key, c in CUSICK.items():
        lo, hi = c["J_cm"]
        theirs[key] = {
            "label": c["label"],
            "method": c["method"],
            "R_A": c["R_A"],
            "delta_deg": c["delta_deg"],
            "unscreened_cm": [lo, hi],
            "screened_cm": [lo / EPSILON, hi / EPSILON],
        }

    payload = {
        "epsilon": EPSILON,
        "ours_production_tandem": ours,
        "cusick": theirs,
        "cusick_table_s2_unscreened_cm": CUSICK_TABLE_S2_CM,
        "cusick_table_s2_screened_cm": [
            CUSICK_TABLE_S2_CM[0] / EPSILON, CUSICK_TABLE_S2_CM[1] / EPSILON
        ],
        "superseded_pair": LEGACY,
        "verdict": (
            "The 32.82-vs-32-40 agreement compares a screened number to an "
            "unscreened one. In a single convention the methods differ by "
            "roughly the screening factor itself. epsilon is now the dominant "
            "systematic."
        ),
    }
    (OUT_DIR / "screening_reconciliation.json").write_text(json.dumps(payload, indent=2))

    L = []
    A = L.append
    A("# Screening reconciliation: this work vs Cusick et al. 2026\n")
    A("All couplings in cm^-1. `epsilon = 1.77` is the optical dielectric used")
    A("throughout this repo; Cusick apply no screening at all.\n")
    A("## The comparison in both conventions\n")
    A("| Quantity | Structure | R (A) | Screened (eps=1.77) | Unscreened |")
    A("|---|---|---|---|---|")
    A(f"| This work, TDC | tandem MD, n={n} | {ours['R_A_mean']:.2f} +/- {ours['R_A_std']:.2f} "
      f"| **{ours['tdc_screened_cm']:.2f} +/- {ours['tdc_screened_std_cm']:.2f}** "
      f"| {ours['tdc_unscreened_cm']:.2f} |")
    A(f"| This work, PDA | tandem MD, n={n} | {ours['R_A_mean']:.2f} +/- {ours['R_A_std']:.2f} "
      f"| {ours['pda_screened_cm']:.2f} +/- {ours['pda_screened_std_cm']:.2f} "
      f"| {ours['pda_unscreened_cm']:.2f} |")
    for key in ("vdw", "tandem"):
        t = theirs[key]
        A(f"| Cusick, {t['label']} | {t['method']} | {t['R_A']:.1f} "
          f"| {fmt_range(*t['screened_cm'])} | **{fmt_range(*t['unscreened_cm'])}** |")
    lo, hi = CUSICK_TABLE_S2_CM
    A(f"| Cusick Table S2 (1MYW) | -- | -- | {fmt_range(lo/EPSILON, hi/EPSILON)} "
      f"| {fmt_range(lo, hi)} |")
    A("")
    A("Bold marks the number each group actually reports. Reading across the bold")
    A("cells compares a screened value against an unscreened one -- that is the")
    A("entire content of the apparent agreement.\n")
    A("## Like for like\n")
    A(f"- Screened throughout: ours {ours['tdc_screened_cm']:.1f} vs theirs "
      f"{fmt_range(*theirs['tandem']['screened_cm'])} -- we are HIGH by ~"
      f"{ours['tdc_screened_cm']/np.mean(theirs['tandem']['screened_cm']):.2f}x.")
    A(f"- Unscreened throughout: ours {ours['tdc_unscreened_cm']:.1f} vs theirs "
      f"{fmt_range(*theirs['tandem']['unscreened_cm'])} -- we are HIGH by ~"
      f"{ours['tdc_unscreened_cm']/np.mean(theirs['tandem']['unscreened_cm']):.2f}x.")
    A("- The discrepancy is the same in either convention, as it must be: the")
    A("  screening factor cancels from the ratio. What does not survive is the")
    A("  claim of a ~1% coincidence, which existed only because the two numbers")
    A("  were quoted in different conventions.\n")
    A("## The point-dipole partner\n")
    A(f"- Genuine partner of {ours['tdc_screened_cm']:.2f}: **{ours['pda_screened_cm']:.2f} cm^-1** "
      f"(same {n} frames, same density, same convention), ratio "
      f"**{ours['tdc_over_pda']:.2f}**.")
    A(f"- NOT {LEGACY['J_pda_cm']}, which pairs with the superseded "
      f"{LEGACY['J_tdc_cm']} on {LEGACY['structure']}, a {LEGACY['density']}, "
      f"{LEGACY['frames']} frame. reproduce_paper.py:115 records that TDC value as "
      f"corrected to {LEGACY['J_tdc_corrected_cm']:.2f}.")
    A(f"- So the near-field enhancement is {ours['tdc_over_pda']:.2f}x, not 5.6x. At "
      f"R ~ {ours['R_A_mean']:.1f} A the point-dipole approximation is accurate to "
      f"~{abs(1-ours['tdc_over_pda'])*100:.0f}%, which makes cheap point-dipole")
    A("  estimates adequate for fluorescent-protein dimers at biological separations.\n")
    A("## What survives untouched\n")
    A(f"- sigma/<J> = {ours['sigma_over_mean']*100:.1f}% over {n} frames at 293 K. This is a")
    A("  ratio, so it is invariant under both the units correction and the choice")
    A("  of epsilon. The coupling is remarkably rigid against thermal motion.\n")
    (OUT_DIR / "screening_reconciliation.md").write_text("\n".join(L))

    print("\n".join(L))
    print(f"\nwrote {OUT_DIR/'screening_reconciliation.md'}")
    print(f"wrote {OUT_DIR/'screening_reconciliation.json'}")


if __name__ == "__main__":
    main()

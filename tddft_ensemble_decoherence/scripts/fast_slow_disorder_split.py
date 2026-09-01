#!/usr/bin/env python3
"""Same-trajectory fast/slow split of the site-energy disorder.

This replaces the approximate split that `cd_width_with_computed_disorder.py`
carries at :23-30, which differences a 1 ns run against a 592 fs run and warns
in its own caveats that the two come from DIFFERENT trajectories, so its
`sigma_site = 159` is approximate.

Here both halves share one lineage: the v6t2 windows are 1 fs continuations
seeded from frames of the very v5 trajectory that supplies the 1 ns statistics,
built from the same solvated system (verified by potential-energy identity to
~0.1 kJ/mol out of 1.2e6). They are continuations rather than literal
sub-segments -- v5 saved only every 0.5 ps, so it holds no 1 fs data -- but the
system, force field, thermostat and box are identical, which the old comparison
could not claim.

Two things are measured rather than assumed:

  * sigma_gap directly. The CD forward model currently builds it as
    sqrt(2) * sigma_site, which presumes the two sites fluctuate independently.
    The gap is what the exciton actually sees, so it is reported here directly
    AND as the sqrt(2) construction, with the site-site correlation printed so
    the assumption can be checked instead of inherited.

  * the fast/slow division at the window length, from
    sigma_slow^2 = sigma_total^2 - sigma_fast^2.

NOTE ON WHAT THIS IS FOR. It is deliberately NOT a dephasing time. Modes with
hbar*w < kT have period > 160 fs and are quasi-static over any plausible
coherence decay; modes fast enough to dephase within ~20 fs lie above 834 cm^-1
where hbar*w/kT > 4 and they are quantum-frozen, producing vibronic structure
rather than homogeneous broadening. The classical trajectory therefore supports
a disorder split, which feeds an experimentally validated lineshape model, but
not a homogeneous T2.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EV = 8065.544005
WINDOWS = {"v6t2a": 900, "v6t2b": 1400, "v6t2c": 1900}


def load(window: str) -> dict:
    path = ROOT / "manifests" / f"{window}_frames.csv"
    with path.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle)
                if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))
    return {
        "ea": np.array([float(r["site_A_energy_eV"]) for r in rows]) * EV,
        "eb": np.array([float(r["site_B_energy_eV"]) for r in rows]) * EV,
        "n": len(rows),
    }


def stats(ea: np.ndarray, eb: np.ndarray) -> dict:
    gap = ea - eb
    r = float(np.corrcoef(ea, eb)[0, 1]) if len(ea) > 2 else float("nan")
    sa, sb = float(ea.std(ddof=1)), float(eb.std(ddof=1))
    return {
        "sigma_site_A_cm": sa,
        "sigma_site_B_cm": sb,
        "sigma_site_rms_cm": float(np.sqrt((sa**2 + sb**2) / 2.0)),
        "site_site_r": r,
        "sigma_gap_measured_cm": float(gap.std(ddof=1)),
        "sigma_gap_sqrt2_construction_cm": float(np.sqrt(2) * np.sqrt((sa**2 + sb**2) / 2.0)),
        "n": int(len(ea)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/fast_slow_disorder_split.json")
    args = parser.parse_args()

    # ---- FAST: within the 1 fs windows -------------------------------------
    per_window = {}
    ea_all, eb_all = [], []
    for w in WINDOWS:
        d = load(w)
        if d["n"] < 100:
            raise RuntimeError(f"{w}: only {d['n']}/100 complete; run after the campaign")
        per_window[w] = stats(d["ea"], d["eb"])
        # centre each window on its own mean before pooling: the windows sit at
        # different detunings by construction, and that offset is slow disorder,
        # not fast.
        ea_all.append(d["ea"] - d["ea"].mean())
        eb_all.append(d["eb"] - d["eb"].mean())
    fast = stats(np.concatenate(ea_all), np.concatenate(eb_all))

    # ---- TOTAL: the 1 ns parent ensemble -----------------------------------
    with (ROOT / "manifests/v5camp_frames.csv").open(newline="") as handle:
        v5 = [r for r in csv.DictReader(handle) if r["convergence_status"] == "complete"]
    total = stats(np.array([float(r["site_A_energy_eV"]) for r in v5]) * EV,
                  np.array([float(r["site_B_energy_eV"]) for r in v5]) * EV)

    def split(tot: float, fst: float) -> float:
        return float(np.sqrt(max(tot**2 - fst**2, 0.0)))

    slow_gap = split(total["sigma_gap_measured_cm"], fast["sigma_gap_measured_cm"])
    slow_site = split(total["sigma_site_rms_cm"], fast["sigma_site_rms_cm"])

    print("=== same-trajectory fast/slow disorder split ===\n")
    print(f"{'':34s} {'site A':>9s} {'site B':>9s} {'site rms':>9s} "
          f"{'gap meas':>9s} {'sqrt2*site':>11s} {'r(A,B)':>8s}")
    for lab, s in (("FAST  (within 100 fs, 3 windows)", fast),
                   ("TOTAL (1 ns v5 ensemble)", total)):
        print(f"{lab:34s} {s['sigma_site_A_cm']:9.1f} {s['sigma_site_B_cm']:9.1f} "
              f"{s['sigma_site_rms_cm']:9.1f} {s['sigma_gap_measured_cm']:9.1f} "
              f"{s['sigma_gap_sqrt2_construction_cm']:11.1f} {s['site_site_r']:8.3f}")
    print(f"{'SLOW  (quadrature difference)':34s} {'':9s} {'':9s} "
          f"{slow_site:9.1f} {slow_gap:9.1f}")

    print("\n  per window (fast only):")
    for w, s in per_window.items():
        print(f"    {w} (seed {WINDOWS[w]:4d}): sigma_gap {s['sigma_gap_measured_cm']:6.1f}  "
              f"r(A,B) {s['site_site_r']:+.3f}")

    print("\n=== what this replaces in the CD forward model ===")
    print(f"  cd_width_with_computed_disorder.py:60  SIGMA_SITE_SLOW = 159.0  (approximate,")
    print(f"     from two different trajectories)")
    print(f"  measured here, same trajectory:        sigma_site_slow = {slow_site:.1f}")
    print(f"  sigma_gap fed to the model:")
    print(f"     current construction sqrt(2)*159.0 = {np.sqrt(2)*159.0:.1f} cm-1 "
          f"(assumes uncorrelated sites)")
    print(f"     measured slow gap directly         = {slow_gap:.1f} cm-1")
    rt = total["site_site_r"]
    print(f"     site-site r over 1 ns = {rt:+.3f}; the sqrt(2) construction is exact only at r = 0")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({
        "fast_within_100fs": fast,
        "total_1ns": total,
        "slow_gap_cm": slow_gap,
        "slow_site_cm": slow_site,
        "per_window_fast": per_window,
        "supersedes": {
            "file": "terachem_site_energy_cd/scripts/cd_width_with_computed_disorder.py",
            "constant": "SIGMA_SITE_SLOW = 159.0",
            "why": "that value differences a 1 ns run against a 592 fs run from a "
                   "different trajectory; both halves here share one lineage",
        },
        "not_a_dephasing_time": "hbar*w < kT implies period > 160 fs (quasi-static over "
                                "any plausible coherence); modes fast enough to dephase "
                                "lie above 834 cm^-1 where they are quantum-frozen and "
                                "produce vibronic structure, not homogeneous broadening",
    }, indent=2))
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()

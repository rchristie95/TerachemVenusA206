#!/usr/bin/env python3
"""Ensemble absolute rotational strength of the CD couplet, from the v5 frames.

The manuscript quotes an integrated ensemble rotational strength for the exciton
couplet, computed from the superseded rigid-fit ensemble. This recomputes it
from the 500-frame per-frame-density campaign using the same absolute_cd.py
machinery, so the convention is unchanged.

THE PHASE PROBLEM. Per frame this needs the triple product
R_AB . (mu_A x mu_B), and a TDDFT transition dipole carries an arbitrary global
sign: TeraChem's eigenvector phase is not fixed between frames. Taking the
dipoles as printed makes the triple product's sign random -- it came out
negative in 46% of frames, i.e. the ensemble mean was averaging a coin flip and
was meaningless.

The coupling does not have this problem, because it is computed from densities
that the state tracking aligns against a common reference, which is why every
stored J is positive. To combine the two, each dipole's sign is fixed here in
the CHROMOPHORE's own frame, by requiring a positive projection on the vector
between two fixed heteroatoms (indices 6 and 22, 8.55 A apart). Atom ordering is
identical in all 1000 geometry files, so this is the same physical convention at
both sites and in every frame.

Validation: with the phases fixed this way, the point-dipole coupling rebuilt
from the phased dipoles has a single consistent sign in 100% of frames and its
ratio to the stored J_pda is constant to 2% (the residual is the QM-atom
centroid standing in for the transition-density centroid). The relative phase is
therefore the one the density calculation used.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
sys.path.insert(0, str(REPO / "terachem_site_energy_cd"))
from absolute_cd import rotational_strengths, triple_product_cgs  # noqa: E402

PHASE_ATOMS = (6, 22)     # fixed heteroatom pair defining the molecular axis
N_PHYSICAL = 41           # the last three atoms are link caps


def geometry(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    return np.array([list(map(float, l.split()[1:4])) for l in lines[2:2 + n]])


def main() -> None:
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/v5camp_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))

    recs = []
    for r in rows:
        d = ROOT / "frames/v5camp" / f"frame_{int(r['frame_index']):04d}"
        s = json.loads((d / "frame_summary.json").read_text())
        ga = geometry(d / "site_A/geometry.xyz")
        gb = geometry(d / "site_B/geometry.xyz") + np.asarray(
            s["tdc"]["site_B_lattice_shift_A"], float)
        k1, k2 = PHASE_ATOMS
        mu_a = np.asarray(s["site_A"]["transition_dipole_au"], float)
        mu_b = np.asarray(s["site_B"]["transition_dipole_au"], float)
        mu_a = mu_a * np.sign(np.dot(mu_a, ga[k2] - ga[k1]))
        mu_b = mu_b * np.sign(np.dot(mu_b, gb[k2] - gb[k1]))
        r_a, r_b = ga[:N_PHYSICAL].mean(0), gb[:N_PHYSICAL].mean(0)
        recs.append((float(r["site_A_energy_cm-1"]), float(r["site_B_energy_cm-1"]),
                     float(r["J_screened_cm-1"]),
                     triple_product_cgs(r_a, r_b, mu_a, mu_b)))

    T = np.array([t for *_, t in recs])
    J = np.array([j for _, _, j, _ in recs])
    D = np.array([a - b for a, b, _, _ in recs])
    Om = np.hypot(D, 2 * J)
    Rp = np.array([rotational_strengths(a, b, j, t)[2] for a, b, j, t in recs])
    sem = Rp.std(ddof=1) / np.sqrt(Rp.size)

    print(f"  n = {Rp.size} frames")
    print(f"  triple product  <T>       = {T.mean():+.4e} cgs, "
          f"same sign in {100*max(np.mean(T > 0), np.mean(T < 0)):.0f}% of frames")
    print(f"  suppression     <2|J|/Om> = {np.mean(2*J/Om):.4f}  "
          f"(median {np.median(2*J/Om):.4f})")
    print(f"  rotational strength <R_+> = {Rp.mean():+.4e} +/- {sem:.2e} cgs (sem)")
    print(f"                       sd   = {Rp.std(ddof=1):.4e} cgs")
    print(f"                       |mean|/sem = {abs(Rp.mean())/sem:.1f}")
    print(f"  manuscript quoted from the rigid-fit ensemble: -2.2e-38 cgs")


if __name__ == "__main__":
    main()

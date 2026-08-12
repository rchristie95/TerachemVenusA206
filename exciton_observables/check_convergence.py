#!/usr/bin/env python3
"""Is the long unbiased MD converging to the crystal geometry and to 57 ps?

Runs the full chain on any protein-only DCD from run_overnight.py:
  DCD -> per-frame rigid CR2 fits -> rigid STEOM density placement
      -> cos(alpha), J -> thermally corrected superradiant lifetime shift.

The angle logged live in <arm>_metrics.csv is the cheap route-3 axis
(OH minus imidazolinone-ring centroid) and carries a geometry-dependent offset
from the STEOM-density route used for the observable -- +3.9 deg at the crystal,
+0.7 deg in the MD. Monitoring is fine on route 3; any number that gets quoted
has to come through this script.

Usage:
  python check_convergence.py                       # the running continuation
  python check_convergence.py --dcd path/to.dcd
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path("/home/robson/PetaChem")
HERE = REPO / "exciton_observables"
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.63
KT_CM = 208.509
PHI = 0.57
TAU_PS = 3026.0
EPSILON = 1.77
TDC_OVER_PDA = 1.1872
DELTA_CM = 570.0          # computed QM/MM detuning (549 / 576 / 581 across ensembles)
MEASURED_PS = 57.0
MEASURED_SD_PS = 4.0
CRYSTAL = {"alpha": 110.44, "cos": -0.349, "J": 30.82, "dtau": 57.0}

DENSITY = REPO / "neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz"
MONOMER = REPO / "tc_simple_old/classical_relaxed.pdb"
TOPOLOGY = REPO / "tc_tandem_nvt_v3_ff19sb_opc/protein_topology.pdb"


def predicted_dtau(cos_alpha, j_cm, delta_cm=DELTA_CM):
    omega = np.sqrt(delta_cm ** 2 + 4.0 * j_cm ** 2)
    return PHI * np.tanh(omega / (2.0 * KT_CM)) * (2.0 * j_cm / omega) * (-cos_alpha) * TAU_PS


def geometry_from_transforms(npz_path, box_a=None):
    d = np.load(DENSITY)
    pts, q = d["pts_ang"], d["q"]
    mu_local = ((pts - pts.mean(0)) * q[:, None]).sum(0) * ANGSTROM_TO_BOHR
    origin = pts.mean(0)

    t = np.load(npz_path)
    rot, trans = t["rotation"], t["translation"]
    mu_a = np.einsum("fij,j->fi", rot[:, 0], mu_local)
    mu_b = np.einsum("fij,j->fi", rot[:, 1], mu_local)
    o_a = np.einsum("fij,j->fi", rot[:, 0], origin) + trans[:, 0]
    o_b = np.einsum("fij,j->fi", rot[:, 1], origin) + trans[:, 1]

    cos_alpha = np.sum(mu_a * mu_b, axis=1) / (
        np.linalg.norm(mu_a, axis=1) * np.linalg.norm(mu_b, axis=1))
    sep_vec = o_b - o_a
    if box_a:  # protein-only DCDs are written wrapped; a two-chain dimer can
        # land in different periodic images, which inflates the separation by
        # ~box and collapses J. The angle is unaffected (rotations only).
        sep_vec = sep_vec - box_a * np.round(sep_vec / box_a)
    r_vec = sep_vec * ANGSTROM_TO_BOHR
    r = np.linalg.norm(r_vec, axis=1)
    r_hat = r_vec / r[:, None]
    jdd = (np.sum(mu_a * mu_b, axis=1)
           - 3.0 * np.sum(mu_a * r_hat, axis=1) * np.sum(mu_b * r_hat, axis=1))
    j_cm = jdd / (r ** 3 * EPSILON) * HARTREE_TO_CM * TDC_OVER_PDA
    return cos_alpha, j_cm, np.linalg.norm(sep_vec, axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dcd", type=Path,
                    default=REPO / "overnight_linker_release/control_long/control_long.dcd")
    ap.add_argument("--blocks", type=int, default=6, help="Report this many time blocks.")
    ap.add_argument("--topology", type=Path, default=TOPOLOGY,
                    help="Protein-only topology matching the DCD atom order.")
    ap.add_argument("--frame-ps", type=float, default=25.0, help="Spacing between DCD frames.")
    ap.add_argument("--box-a", type=float, default=None,
                    help="Cubic box edge for minimum-image correction. Required for\n"
                         "multi-chain dimers, whose chains can wrap into different images.")
    args = ap.parse_args()

    if not args.dcd.exists():
        sys.exit(f"No trajectory at {args.dcd}")

    transforms = HERE / f"transforms_{args.dcd.stem}.npz"
    print(f"[*] extracting CR2 transforms from {args.dcd.name} ...")
    subprocess.run([sys.executable, str(REPO / "extract_cr2_transforms.py"),
                    "--traj", str(args.dcd), "--topology", str(args.topology),
                    "--monomer", str(MONOMER), "--out", str(transforms)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    cos_alpha, j_cm, sep = geometry_from_transforms(transforms, args.box_a)
    dtau = predicted_dtau(cos_alpha, j_cm)
    angle = np.degrees(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))
    n = len(cos_alpha)

    print(f"\n{n} frames (25 ps apart = {n*args.frame_ps/1000:.1f} ns)\n")
    print(f"{'block':>14}{'alpha':>9}{'cos a':>9}{'J':>8}{'sep':>8}{'pred dtau':>11}")
    for i, sl in enumerate(np.array_split(np.arange(n), args.blocks)):
        lo, hi = sl[0]*args.frame_ps/1000, (sl[-1]+1)*args.frame_ps/1000
        print(f"{f'{lo:.0f}-{hi:.0f} ns':>14}{angle[sl].mean():9.2f}"
              f"{cos_alpha[sl].mean():+9.3f}{j_cm[sl].mean():8.2f}"
              f"{sep[sl].mean():8.2f}{dtau[sl].mean():11.1f}")

    tail = slice(n // 2, None)
    blocks = np.array_split(dtau[tail], 5)
    sem = np.std([b.mean() for b in blocks], ddof=1) / np.sqrt(5)
    print(f"\n  last half:   alpha {angle[tail].mean():.2f} deg, "
          f"predicted Dtau {dtau[tail].mean():.1f} +/- {sem:.1f} ps")
    print(f"  crystal:     alpha {CRYSTAL['alpha']:.2f} deg, "
          f"predicted Dtau {CRYSTAL['dtau']:.1f} ps")
    print(f"  measured:                    Dtau {MEASURED_PS:.1f} +/- {MEASURED_SD_PS:.1f} ps")

    gap = CRYSTAL["alpha"] - angle[tail].mean()
    print(f"\n  still {gap:+.2f} deg from the crystal interface; "
          f"{'converging' if gap > 0 else 'at or past'} it.")


if __name__ == "__main__":
    main()

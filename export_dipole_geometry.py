#!/usr/bin/env python3
r"""
export_dipole_geometry.py  --  Export the physical two-dipole geometry of the
excitonic dimer for the CD/absorption lineshape (absorption_cd_spectra.py).

The CD rotational strength of the exciton couplet is
    R_+- = -+ (pi nu0 / 2) * R_AB . (mu_A x mu_B),
so it depends on the two monomer transition dipoles and centroids placed at the
real chromophore sites. This script produces those four vectors from the STEOM
transition density using the SAME Kabsch/`super` placement as the coupling
pipeline (coupling_ensemble.py rigid mode / align_steom_density.py), so the
lineshape geometry is consistent with the reported J.

Placement (identical to coupling_ensemble.coupling_for_frame):
  1. Load the spectroscopically-normalised STEOM density (npz: pts_ang, q, mu_au),
     already expressed in the "old" monomer frame by align_steom_density.py.
  2. PyMOL `super` the monomer onto chain A and chain B of the reference dimer
     to get the two site transforms, and map the fixed density onto both sites.
  3. r_A, r_B  = density centroids (Angstrom);
     mu_A, mu_B = point-charge transition dipoles about each centroid (a.u.).

Output JSON (read by absorption_cd_spectra.py --geometry-json):
    { "mu_A":[...], "mu_B":[...], "r_A":[...], "r_B":[...],   # a.u. / Angstrom
      "separation_A": ..., "angle_deg": ..., "mu_debye": ... }

Example:
    python export_dipole_geometry.py \
        --density neo_model/orca_steom/steom_transdens_capmasked_oldframe.npz \
        --monomer tc_simple_old/classical_relaxed.pdb \
        --dimer   venus_dimer.pdb \
        --out     coupling_paper_steom_thermal/dipole_geometry.json
"""
import argparse
import json
from pathlib import Path

import numpy as np

from coupling_core import (
    get_super_matrices_with_pymol,
    apply_pymol_matrix,
    transition_dipole_au,
)

# 1 a.u. of dipole moment = 2.541746 Debye.
AU_TO_DEBYE = 2.541746230


def export_geometry(density_npz, monomer_pdb, dimer_pdb):
    d = np.load(density_npz)
    pts = np.ascontiguousarray(d["pts_ang"], dtype=float)
    q = np.ascontiguousarray(d["q"], dtype=float)

    matrix_A, matrix_B, aln_A, aln_B, err = get_super_matrices_with_pymol(
        str(monomer_pdb), str(dimer_pdb))
    if err:
        raise RuntimeError(f"PyMOL super alignment failed: {err}")

    pts_A = apply_pymol_matrix(pts, matrix_A)
    pts_B = apply_pymol_matrix(pts, matrix_B)
    r_A = pts_A.mean(axis=0)
    r_B = pts_B.mean(axis=0)

    mu_A = transition_dipole_au(pts_A, q, origin_angstrom=r_A)
    mu_B = transition_dipole_au(pts_B, q, origin_angstrom=r_B)

    sep = float(np.linalg.norm(r_A - r_B))
    cosang = float(np.dot(mu_A, mu_B) /
                   (np.linalg.norm(mu_A) * np.linalg.norm(mu_B) + 1e-30))
    angle = float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))
    mu_debye = float(0.5 * (np.linalg.norm(mu_A) + np.linalg.norm(mu_B)) * AU_TO_DEBYE)

    return {
        "mu_A": mu_A.tolist(),
        "mu_B": mu_B.tolist(),
        "r_A": r_A.tolist(),
        "r_B": r_B.tolist(),
        "separation_A": sep,
        "angle_deg": angle,
        "mu_debye": mu_debye,
        "aln_A_rms": float(aln_A[0]) if aln_A else float("nan"),
        "aln_B_rms": float(aln_B[0]) if aln_B else float("nan"),
        "source_density": str(density_npz),
        "monomer": str(monomer_pdb),
        "dimer": str(dimer_pdb),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--density",
                   default="neo_model/orca_steom/steom_transdens_capmasked_oldframe.npz",
                   help="STEOM density npz in the old-monomer frame (align_steom_density.py).")
    p.add_argument("--monomer", default="tc_simple_old/classical_relaxed.pdb",
                   help="Monomer reference PDB (chain A) for the super alignment.")
    p.add_argument("--dimer", default="venus_dimer.pdb",
                   help="Reference dimer PDB (chains A and B).")
    p.add_argument("--out", type=Path,
                   default=Path("coupling_paper_steom_thermal/dipole_geometry.json"))
    args = p.parse_args(argv)

    info = export_geometry(args.density, args.monomer, args.dimer)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(info, f, indent=2)
        f.write("\n")

    print(f"[dipole-geometry] separation |r_A - r_B| = {info['separation_A']:.2f} A")
    print(f"[dipole-geometry] inter-dipole angle      = {info['angle_deg']:.2f} deg")
    print(f"[dipole-geometry] |mu| (monomer)          = {info['mu_debye']:.2f} D")
    print(f"[dipole-geometry] super RMSD A/B           = "
          f"{info['aln_A_rms']:.3f}/{info['aln_B_rms']:.3f} A")
    print(f"[dipole-geometry] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

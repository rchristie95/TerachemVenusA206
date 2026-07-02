#!/usr/bin/env python3
"""
align_steom_density.py

Place the STEOM-CCSD spec-normalised transition density into the SAME coordinate
frame the TeraChem/TDDFT coupling pipeline uses, so the Davydov coupling J can be
computed for STEOM and TDDFT at identical dimer geometries (only the density
differs). Without this, aligning the STEOM density via the anionic monomer
mis-places it ~4 A and inflates J in the near field.

Method: rigid Kabsch fit of the chromophore (CR2) heavy atoms shared between the
anionic monomer (the frame the STEOM density was built in) and the "old" monomer
(the frame the NVT dimer chains were built from), applied to the density points.
A pure rotation+translation preserves |mu|.

Usage:
    python align_steom_density.py \
        --density   neo_model/orca_steom/steom_transdens_specnorm.npz \
        --anion-pdb tc_simple_anionic/monomer_relaxed.pdb \
        --old-pdb   tc_simple_old/classical_relaxed.pdb \
        --out       neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz
"""
import argparse
from pathlib import Path

import numpy as np


def cr2_atoms(pdb_path):
    """Return {atom_name: xyz(Angstrom)} for CR2 residue atoms in a PDB."""
    d = {}
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() in ("ATOM", "HETATM") and line[17:20].strip() == "CR2":
                d[line[12:16].strip()] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    return d


def kabsch(P, Q):
    """Rigid R,t mapping P -> Q (least-squares). Returns (R, t, rmsd)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    rmsd = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def match_density_to_frame(density_npz, anion_pdb, old_pdb, out_npz):
    anion, old = cr2_atoms(anion_pdb), cr2_atoms(old_pdb)
    common = sorted(set(anion) & set(old))
    if len(common) < 4:
        raise RuntimeError(f"Too few shared CR2 atoms ({len(common)}) for a rigid fit.")
    P = np.array([anion[n] for n in common])
    Q = np.array([old[n] for n in common])
    R, t, rmsd = kabsch(P, Q)

    d = np.load(density_npz)
    pts, q, mu = d["pts_ang"], d["q"], d["mu_au"]
    pts2 = (R @ pts.T).T + t
    mu2 = R @ mu
    Path(out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_npz, pts_ang=pts2, q=q, mu_au=mu2)
    return {
        "n_common_cr2": len(common),
        "fit_rmsd_A": rmsd,
        "mu_before": float(np.linalg.norm(mu)),
        "mu_after": float(np.linalg.norm(mu2)),
        "n_points": int(pts2.shape[0]),
        "out": str(out_npz),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--density", default="neo_model/orca_steom/steom_transdens_specnorm.npz")
    p.add_argument("--anion-pdb", default="tc_simple_anionic/monomer_relaxed.pdb")
    p.add_argument("--old-pdb", default="tc_simple_old/classical_relaxed.pdb")
    p.add_argument("--out", default="neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz")
    args = p.parse_args(argv)
    info = match_density_to_frame(args.density, args.anion_pdb, args.old_pdb, args.out)
    print(f"[matched-density] {info['n_common_cr2']} shared CR2 atoms, "
          f"fit RMSD {info['fit_rmsd_A']:.3f} A, "
          f"|mu| {info['mu_before']:.4f} -> {info['mu_after']:.4f} au (preserved)")
    print(f"[matched-density] wrote {info['out']} ({info['n_points']} points)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

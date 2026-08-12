#!/usr/bin/env python3
"""Inter-chromophore transition-dipole angle for any Venus dimer geometry.

Reproduces the geometry half of the production coupling run (rigid STEOM density
placed by a Kabsch fit on the shared CR2 heavy atoms) with numpy only, so the
angle can be evaluated on structures that were never put through the full TDC
pipeline -- in particular the crystal A206 dimer that Kim measured.

Motivation: Nguyen's limiting-anisotropy drop (r0 0.52 -> 0.30, Fig. 7) implies
|cos alpha| = 0.660, i.e. alpha = 48.7 or 131.3 deg, whereas the production
tandem ensemble sits at 100.78 +/- 2.72 deg (|cos alpha| = 0.187).  The angle
enters the detuning-free absorption red shift as J*cos(alpha), so a factor ~3.5
on |cos alpha| matters.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

ANGSTROM_TO_BOHR = 1.8897259886


def cr2_heavy_atoms(pdb_path, chain=None, resid=None):
    """Map atom name -> xyz for one CR2 residue (heavy atoms only)."""
    atoms = {}
    with open(pdb_path) as handle:
        for line in handle:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip() != "CR2":
                continue
            if chain is not None and line[21] != chain:
                continue
            if resid is not None and line[22:26].strip() != str(resid):
                continue
            element = line[76:78].strip().upper()
            name = line[12:16].strip()
            if element == "H" or (not element and name.upper().startswith("H")):
                continue
            atoms[name] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])], float
            )
    return atoms


def cr2_sites(pdb_path):
    """Every CR2 residue in the file, as (chain, resid, {name: xyz})."""
    sites = []
    seen = []
    with open(pdb_path) as handle:
        for line in handle:
            if line[:6].strip() in ("ATOM", "HETATM") and line[17:20].strip() == "CR2":
                key = (line[21], line[22:26].strip())
                if key not in seen:
                    seen.append(key)
    for chain, resid in seen:
        sites.append((chain, resid, cr2_heavy_atoms(pdb_path, chain, resid)))
    return sites


def kabsch(source, target):
    """Rotation+translation carrying source onto target, plus the fit RMSD."""
    sc = source - source.mean(0)
    tc = target - target.mean(0)
    u, _, vt = np.linalg.svd(sc.T @ tc)
    handedness = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, handedness]) @ u.T
    trans = target.mean(0) - rot @ source.mean(0)
    fitted = source @ rot.T + trans
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, 1))))
    return rot, trans, rmsd


def reference_dipole(density_npz, monomer_pdb):
    """Transition dipole (a.u.) in the reference monomer's own frame."""
    dens = np.load(density_npz)
    pts, q = dens["pts_ang"], dens["q"]
    mu = ((pts - pts.mean(0)) * q[:, None]).sum(0) * ANGSTROM_TO_BOHR
    ref = cr2_heavy_atoms(monomer_pdb)
    return mu, ref


def place(ref_atoms, ref_mu, site_atoms):
    """Rotate the reference dipole onto one CR2 site; return (mu, centroid, rmsd)."""
    names = sorted(set(ref_atoms) & set(site_atoms))
    if len(names) < 12:
        raise RuntimeError(f"only {len(names)} shared CR2 heavy atoms: {names}")
    src = np.array([ref_atoms[n] for n in names])
    dst = np.array([site_atoms[n] for n in names])
    rot, _, rmsd = kabsch(src, dst)
    return rot @ ref_mu, dst.mean(0), rmsd, len(names)


def report(label, pdb_path, ref_atoms, ref_mu):
    sites = cr2_sites(pdb_path)
    if len(sites) != 2:
        print(f"{label:<34} SKIP: {len(sites)} CR2 residues")
        return None
    placed = [place(ref_atoms, ref_mu, s[2]) for s in sites]
    (muA, rA, rmsA, nA), (muB, rB, rmsB, _) = placed
    cos = float(np.dot(muA, muB) / (np.linalg.norm(muA) * np.linalg.norm(muB)))
    ang = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
    sep = float(np.linalg.norm(rA - rB))
    print(
        f"{label:<34} angle {ang:7.2f} deg   cos {cos:+.3f}   "
        f"sep {sep:6.2f} A   fit RMSD {rmsA:.2f}/{rmsB:.2f} A ({nA} atoms)"
    )
    return ang, cos, sep


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--density", default="neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz")
    ap.add_argument("--monomer", default="tc_simple_old/classical_relaxed.pdb")
    ap.add_argument("structures", nargs="+")
    args = ap.parse_args()

    ref_mu, ref_atoms = reference_dipole(args.density, args.monomer)
    print(f"reference monomer : {args.monomer}")
    print(f"reference density : {args.density}")
    print(f"|mu_ref| = {np.linalg.norm(ref_mu):.4f} a.u. "
          f"({np.linalg.norm(ref_mu) * 2.541746:.3f} D), "
          f"{len(ref_atoms)} CR2 heavy atoms\n")
    for path in args.structures:
        report(Path(path).name, path, ref_atoms, ref_mu)


if __name__ == "__main__":
    main()

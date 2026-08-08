#!/usr/bin/env python3
"""Check 1: is J(STEOM)=154 robust, or dominated by low-amplitude diffuse tails?
Sweep the density threshold (trim low-|q| voxels = the diffuse tails) and a grid stride,
recomputing J each time. Stable J across thresholds => physical; sharp drop => tail artifact."""
import numpy as np, sys
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc
HARTREE_CM = 219474.6314; EPS = 1.77
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"

d = np.load("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.npz")
pts, q = d["pts_ang"], d["q"]
mA, mB, aA, aB, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)
ptsA = cc.apply_pymol_matrix(pts, mA); ptsB = cc.apply_pymol_matrix(pts, mB)
qmax = np.abs(q).max()
print(f"{'rel-thresh':>10} {'npts':>8} {'|mu| au':>8} {'J cm^-1':>9}  (full set: {len(q)} pts)")
for relthr in [1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3]:
    keep = np.abs(q) > relthr * qmax
    if keep.sum() < 10: continue
    J = cc.calculate_coupling(ptsA[keep], q[keep], ptsB[keep], q[keep]) / EPS * HARTREE_CM
    mu = np.linalg.norm(cc.transition_dipole_au(pts[keep], q[keep]))
    print(f"{relthr:>10.0e} {int(keep.sum()):>8d} {mu:>8.3f} {J:>9.2f}")

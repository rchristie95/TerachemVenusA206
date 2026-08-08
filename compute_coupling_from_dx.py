#!/usr/bin/env python3
"""Coupling J from any TeraChem transition-density .dx (same alignment as the STEOM calc).
Usage: compute_coupling_from_dx.py <transdens.dx>"""
import numpy as np, sys
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc
HARTREE_CM = 219474.6314; EPS = 1.77
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"

dxfile = sys.argv[1]
pts, q = cc.read_dx(dxfile, threshold=1e-6, stride=1)
mu = np.linalg.norm(cc.transition_dipole_au(pts, q))
mA, mB, aA, aB, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)
ptsA = cc.apply_pymol_matrix(pts, mA); ptsB = cc.apply_pymol_matrix(pts, mB)
J = cc.calculate_coupling(ptsA, q, ptsB, q) / EPS * HARTREE_CM
print(f"RESULT {dxfile}: |mu|={mu:.3f} au | {len(q)} pts | J_TDC={J:.2f} cm^-1 (eps={EPS})")

#!/usr/bin/env python3
"""Excitonic coupling J from the STEOM-CCSD bright-state transition density.
Reuses the unit-corrected coupling_core pipeline (44-atom TDDFT reference:
J_TDC=22.1 and J_PDA=18.0 cm^-1):
align the monomer (in the STEOM density frame) onto dimer chains A/B, transform the density,
Coulomb double-sum. Optical-limit screening eps=1.77 (t=0), as in the numerical model."""
import numpy as np, sys
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc

HARTREE_CM = 219474.6314
EPS = 1.77
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"  # co-frame w/ STEOM density
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"

d = np.load("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.npz")
pts, q = d["pts_ang"], d["q"]

# --- units/consistency check: the loaded density must reproduce the STEOM dipole (~3.75 au) ---
mu = cc.transition_dipole_au(pts, q)
print(f"[check] STEOM density |mu| = {np.linalg.norm(mu):.3f} au (expect ~3.77); "
      f"sum q = {q.sum():+.2e} (expect ~0); {len(q)} grid points")

# --- align monomer -> dimer sites A/B, transform the density ---
matrix_A, matrix_B, aln_A, aln_B, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)
if err:
    print("ALIGNMENT ERROR:", err); sys.exit(1)
print(f"[align] super RMSD A={aln_A}, B={aln_B}")
pts_A = cc.apply_pymol_matrix(pts, matrix_A)
pts_B = cc.apply_pymol_matrix(pts, matrix_B)

# --- coupling ---
J_ha = cc.calculate_coupling(pts_A, q, pts_B, q) / EPS
J_cm = J_ha * HARTREE_CM
print("")
print(f"=== EXCITONIC COUPLING from STEOM-CCSD transition density ===")
print(f"  J_TDC(STEOM) = {J_cm:.2f} cm^-1   (eps_opt={EPS})")
print(f"  2|J|         = {2*abs(J_cm):.1f} cm^-1  (Davydov splitting)")
print(f"  -- revised 44-atom TDDFT reference: J_TDC=22.1 cm^-1, J_PDA=18.0 cm^-1 --")

# point-dipole estimate from the transformed site dipoles
muA = cc.transition_dipole_au(pts_A, q); muB = cc.transition_dipole_au(pts_B, q)
cenA = (pts_A * np.abs(q)[:,None]).sum(0)/np.abs(q).sum()
cenB = (pts_B * np.abs(q)[:,None]).sum(0)/np.abs(q).sum()
R = (cenB - cenA) * cc.ANGSTROM_TO_BOHR; Rn = np.linalg.norm(R); Rhat = R/Rn
Jdd = (np.dot(muA,muB) - 3*np.dot(muA,Rhat)*np.dot(muB,Rhat)) / Rn**3 / EPS
print(f"  J_PDA(STEOM) = {Jdd*HARTREE_CM:.2f} cm^-1   (centroid sep {Rn*cc.BOHR_TO_ANGSTROM:.1f} A)"
      if hasattr(cc,'BOHR_TO_ANGSTROM') else f"  J_PDA(STEOM) = {Jdd*HARTREE_CM:.2f} cm^-1")
print(f"  finite-density correction TDC/PDA = {abs(J_cm/(Jdd*HARTREE_CM)):.2f}x")

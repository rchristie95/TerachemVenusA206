#!/usr/bin/env python3
"""STEOM-CCSD excitonic coupling J — numba-parallel CPU kernel (GPU driver is down).
Kernel is identical to coupling_core.coupling_cuda_kernel (r-floor 0.1 Ang),
final result * ANGSTROM_TO_BOHR — directly comparable to published J_TDDFT=74.38 cm^-1."""
import numpy as np, sys, time
from numba import njit, prange, set_num_threads
set_num_threads(16)  # leave cores for the running ORCA job (8 ranks)
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc

HARTREE_CM = 219474.6314
EPS = 1.78
A2B = cc.ANGSTROM_TO_BOHR
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"

@njit(parallel=True, fastmath=True, cache=True)
def coupling_numba(p1, q1, p2, q2):
    n1 = p1.shape[0]; n2 = p2.shape[0]; total = 0.0
    for i in prange(n1):
        xi = p1[i,0]; yi = p1[i,1]; zi = p1[i,2]; qi = q1[i]
        acc = 0.0
        for j in range(n2):
            dx = xi-p2[j,0]; dy = yi-p2[j,1]; dz = zi-p2[j,2]
            r = (dx*dx+dy*dy+dz*dz)**0.5
            if r < 0.1: r = 0.1
            acc += qi*q2[j]/r
        total += acc
    return total

d = np.load("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.npz")
pts, q = np.ascontiguousarray(d["pts_ang"],np.float64), np.ascontiguousarray(d["q"],np.float64)
mu = cc.transition_dipole_au(pts, q)
print(f"[check] STEOM |mu|={np.linalg.norm(mu):.3f} au (expect ~3.77); sum q={q.sum():+.2e}; {len(q)} pts", flush=True)

matrix_A, matrix_B, aln_A, aln_B, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)
if err: print("ALIGNMENT ERROR:", err); sys.exit(1)
print(f"[align] RMSD A={aln_A[0]:.3f} A, B={aln_B[0]:.3f} A", flush=True)
pts_A = np.ascontiguousarray(cc.apply_pymol_matrix(pts, matrix_A), np.float64)
pts_B = np.ascontiguousarray(cc.apply_pymol_matrix(pts, matrix_B), np.float64)

print("[run] JIT compile + double-sum...", flush=True)
t0 = time.time()
jsum = coupling_numba(pts_A, q, pts_B, q)
J_ha = jsum * A2B / EPS
J_cm = J_ha * HARTREE_CM
print(f"[run] done in {time.time()-t0:.1f}s", flush=True)

# point-dipole cross-check
muA = cc.transition_dipole_au(pts_A, q); muB = cc.transition_dipole_au(pts_B, q)
cenA = (pts_A*np.abs(q)[:,None]).sum(0)/np.abs(q).sum()
cenB = (pts_B*np.abs(q)[:,None]).sum(0)/np.abs(q).sum()
R = (cenB-cenA)*A2B; Rn = np.linalg.norm(R); Rhat = R/Rn
Jdd = (np.dot(muA,muB)-3*np.dot(muA,Rhat)*np.dot(muB,Rhat))/Rn**3/EPS

print("", flush=True)
print("=== EXCITONIC COUPLING from STEOM-CCSD transition density ===", flush=True)
print(f"  J_TDC(STEOM) = {J_cm:.2f} cm^-1   (eps_opt={EPS})", flush=True)
print(f"  2|J|         = {2*abs(J_cm):.1f} cm^-1  (Davydov splitting)", flush=True)
print(f"  J_PDA(STEOM) = {Jdd*HARTREE_CM:.2f} cm^-1   (centroid sep {Rn/A2B:.1f} A)", flush=True)
print(f"  near-field enhancement TDC/PDA = {abs(J_cm/(Jdd*HARTREE_CM)):.1f}x", flush=True)
print(f"  -- reference: J_TDC(TDDFT, published)=74.38 cm^-1, 2|J|=149, exp 131-186 --", flush=True)

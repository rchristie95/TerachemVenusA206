#!/usr/bin/env python3
"""Validate the unit-corrected numba kernel on the TDDFT .dx densities.

The Coulomb sum uses inverse-Angstrom distances and therefore multiplies by
BOHR_TO_ANGSTROM. The revised 44-atom bright-state reference is 22.1 cm^-1.
"""
import numpy as np, sys, time, glob
from numba import njit, prange, set_num_threads
set_num_threads(16)
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc

HARTREE_CM = 219474.6314; EPS = 1.77; INV_A_TO_INV_BOHR = cc.BOHR_TO_ANGSTROM
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"

@njit(parallel=True, fastmath=True, cache=True)
def coupling_numba(p1, q1, p2, q2):
    n1=p1.shape[0]; n2=p2.shape[0]; total=0.0
    for i in prange(n1):
        xi=p1[i,0]; yi=p1[i,1]; zi=p1[i,2]; qi=q1[i]; acc=0.0
        for j in range(n2):
            dx=xi-p2[j,0]; dy=yi-p2[j,1]; dz=zi-p2[j,2]
            r=(dx*dx+dy*dy+dz*dz)**0.5
            if r<0.1: r=0.1
            acc += qi*q2[j]/r
        total += acc
    return total

# alignment is the same for every density (monomer frame -> dimer A/B)
mA, mB, aA, aB, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)
if err: print("ALIGN ERR", err); sys.exit(1)
print(f"[align] RMSD A={aA[0]:.3f} B={aB[0]:.3f}", flush=True)

def J_of(pts, q):
    pA=np.ascontiguousarray(cc.apply_pymol_matrix(pts,mA),np.float64)
    pB=np.ascontiguousarray(cc.apply_pymol_matrix(pts,mB),np.float64)
    return coupling_numba(pA,q,pB,q)*INV_A_TO_INV_BOHR/EPS*HARTREE_CM

print("\n=== TDDFT states (tc_tddft_44/transdens_N.dx), revised reference 22.1 cm^-1 ===", flush=True)
for f in sorted(glob.glob("/home/robson/PetaChem/tc_tddft_44/transdens_*.dx")):
    pts,q = cc.read_dx(f, threshold=1e-6, stride=1)
    pts=np.ascontiguousarray(pts,np.float64); q=np.ascontiguousarray(q,np.float64)
    mu=np.linalg.norm(cc.transition_dipole_au(pts,q))
    t0=time.time(); J=J_of(pts,q)
    print(f"  {f.split('/')[-1]:18s} |mu|={mu:5.3f} au  npts={len(q):7d}  J_TDC={J:8.2f} cm^-1  ({time.time()-t0:.1f}s)", flush=True)

# STEOM for side-by-side (already-built npz)
d=np.load("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.npz")
pts,q=np.ascontiguousarray(d["pts_ang"],np.float64),np.ascontiguousarray(d["q"],np.float64)
mu=np.linalg.norm(cc.transition_dipole_au(pts,q))
J=J_of(pts,q)
print(f"\n  STEOM bright state   |mu|={mu:5.3f} au  npts={len(q):7d}  J_TDC={J:8.2f} cm^-1", flush=True)

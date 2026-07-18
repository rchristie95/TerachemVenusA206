#!/usr/bin/env python3
"""Check 2 finisher: for each 44-atom TDDFT run, find the bright state (max f), read its
transition-density .dx, and compute J with the same alignment used for STEOM."""
import numpy as np, sys, glob
sys.path.insert(0, "/home/robson/PetaChem")
import coupling_core as cc
HARTREE_CM = 219474.6314; EPS = 1.78
MONOMER = "/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb"
DIMER   = "/home/robson/PetaChem/venus_dimer.pdb"
mA, mB, aA, aB, err = cc.get_super_matrices_with_pymol(MONOMER, DIMER)

def bright(outfile):
    txt = open(outfile, encoding="latin-1").read().splitlines()
    st = []
    for i, l in enumerate(txt):
        if "Final Excited State Results" in l:
            for j in range(i+3, i+40):
                p = txt[j].split()
                if len(p) >= 4 and p[0].isdigit(): st.append((int(p[0]), float(p[2]), float(p[3])))
                elif p and not p[0].lstrip("-").isdigit(): break
            break
    return max(st, key=lambda s: s[2]) if st else (None, None, None)

print(f"{'model':<32}{'bright':>10}{'f':>6}{'|mu|':>7}{'J cm^-1':>9}")
for d, lbl in [("tc_tddft_44", "TDDFT 6-311G** (no diffuse)"),
               ("tc_tddft_44_diff", "TDDFT 6-311++G** (diffuse)")]:
    r, ev, f = bright(f"/home/robson/PetaChem/{d}/td.out")
    cand = glob.glob(f"/home/robson/PetaChem/{d}/transdens_{r}.dx")
    if not cand:
        print(f"{lbl:<32} root {r} ({1239.84/ev:.0f}nm f={f:.2f}) -- transdens .dx missing"); continue
    pts, q = cc.read_dx(cand[0], threshold=1e-6, stride=1)
    ptsA = cc.apply_pymol_matrix(pts, mA); ptsB = cc.apply_pymol_matrix(pts, mB)
    J = cc.calculate_coupling(ptsA, q, ptsB, q) / EPS * HARTREE_CM
    mu = np.linalg.norm(cc.transition_dipole_au(pts, q))
    print(f"{lbl:<32}{1239.84/ev:>8.0f}nm{f:>6.2f}{mu:>7.2f}{J:>9.1f}")
print("--- reference: STEOM-CCSD/def2-SVPD (44-atom) J=154.8 ; published TDDFT/6-311G** (238-atom) J=74.4 ---")
print("CHECK2 DONE")

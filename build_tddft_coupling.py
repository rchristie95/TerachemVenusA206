#!/usr/bin/env python3
"""Build the 44-atom ORCA TDDFT bright-state (S2) transition density from the CIS
transition-density cube, normalise to the ORCA transition dipole (the same
spectroscopic-normalisation step as the STEOM build), reframe into the
dimer-chain frame, and compute J_TDC / J_PDA
on the identical dimer geometry used for the STEOM coupling. This is the
single-excitation cross-check of the near-field correction (revised values:
J_TDC=22.1 and J_PDA=18.0 cm^-1, 1.23x; the definitive STEOM density couples ~1.38x more
strongly). The revised TDC value uses the reciprocal-length conversion
1 Angstrom^-1 = 0.529177 bohr^-1.

Regenerate the S2 transition-density cube from the ORCA TDDFT run (neo_model/orca_dft,
wB97X-D3/6-311G** TDA, DoNTO) with:
    orca_plot tddft_wb97xd3.gbw -i   # 4->100 grid, 7 (transition density), y, 2 (S2), 11, 12
The resulting *.cistp02.cube is the input below (committed as tddft_S2_transdens.cube)."""
import os, sys, numpy as np
REPO = "/home/robson/PetaChem"
sys.path.insert(0, REPO)
import coupling_core as cc
from align_steom_density import match_density_to_frame

BOHR = 0.529177210903
H = 219474.6314; EPS = 1.77
SD = os.path.join(REPO, "neo_model/orca_dft")            # cube + output npz live here
CUBE = os.path.join(SD, "tddft_S2_transdens.cube")
MU_TDDFT = np.array([1.06125, -1.88217, 2.75746])   # au, ORCA S2 transition dipole (|mu|=3.503)

def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()
        toks = f.readline().split(); natom = int(toks[0]); origin = np.array(list(map(float, toks[1:4])))
        n=[]; vec=[]
        for _ in range(3):
            t=f.readline().split(); n.append(int(t[0])); vec.append(list(map(float,t[1:4])))
        n=np.array(n); vec=np.array(vec)
        for _ in range(abs(natom)): f.readline()
        if natom < 0: f.readline()
        data=np.fromstring(" ".join(f.read().split()), sep=" ")
    return origin, n, vec, data.reshape(n)

o,n,vec,rho = read_cube(CUBE)
dvol = abs(np.linalg.det(vec))
ii,jj,kk = np.meshgrid(np.arange(n[0]),np.arange(n[1]),np.arange(n[2]),indexing="ij")
pts_bohr = o[None,None,None,:] + ii[...,None]*vec[0] + jj[...,None]*vec[1] + kk[...,None]*vec[2]
q = (rho*dvol).reshape(-1); pts_bohr=pts_bohr.reshape(-1,3)
mu_raw = -(pts_bohr*q[:,None]).sum(0)
scale = np.dot(MU_TDDFT, mu_raw)/np.dot(mu_raw, mu_raw)
mu_scaled = mu_raw*scale
print(f"grid {tuple(n)}  voxel {dvol:.4f} bohr^3   sum q = {q.sum():+.2e}")
print(f"raw cistp dipole |mu| = {np.linalg.norm(mu_raw):.3f} au  (ORCA target {np.linalg.norm(MU_TDDFT):.3f})")
print(f"Spectroscopic dipole scale = {scale:.4f}   cos(theta)={np.dot(mu_raw,MU_TDDFT)/(np.linalg.norm(mu_raw)*np.linalg.norm(MU_TDDFT)):+.4f}")
q_norm = q*scale; pts_ang = pts_bohr*BOHR
thr = 1e-6*np.abs(q_norm).max(); keep=np.abs(q_norm)>thr
npz_anion = f"{SD}/tddft_transdens_specnorm.npz"
np.savez(npz_anion, pts_ang=pts_ang[keep], q=q_norm[keep], mu_au=mu_scaled)
print(f"saved {keep.sum()} pts (anion/cube frame) -> tddft_transdens_specnorm.npz")

# reframe cube(anion) -> old (dimer-chain) frame, exactly like STEOM
npz_old = f"{SD}/tddft_transdens_specnorm_oldframe.npz"
info = match_density_to_frame(npz_anion, "tc_simple_anionic/monomer_relaxed.pdb",
                              "tc_simple_old/classical_relaxed.pdb", npz_old)
print(f"reframe: {info['n_common_cr2']} CR2 atoms, RMSD {info['fit_rmsd_A']:.3f} A, "
      f"|mu| {info['mu_before']:.3f}->{info['mu_after']:.3f} (preserved)")

# couple: OLD_MONOMER -> dimer placement + Coulomb double sum + PDA (same as STEOM)
d = np.load(npz_old); pts,qq = d["pts_ang"].astype(float), d["q"].astype(float)
mA,mB,aA,aB,e = cc.get_super_matrices_with_pymol("tc_simple_old/classical_relaxed.pdb","venus_dimer.pdb"); assert not e,e
pA=cc.apply_pymol_matrix(pts,mA); pB=cc.apply_pymol_matrix(pts,mB)
J = cc.calculate_coupling(pA,qq,pB,qq,backend="opencl")*H/EPS
muA=cc.transition_dipole_au(pA,qq); muB=cc.transition_dipole_au(pB,qq)
cA=pA.mean(0)*cc.ANGSTROM_TO_BOHR; cB=pB.mean(0)*cc.ANGSTROM_TO_BOHR; R=cB-cA; Rn=np.linalg.norm(R); Rh=R/Rn
pda=(np.dot(muA,muB)-3*np.dot(muA,Rh)*np.dot(muB,Rh))/Rn**3*H/EPS
print("")
print(f"=== 44-ATOM ORCA TDDFT (S2) TRANSITION-DENSITY COUPLING (eps={EPS}) ===")
print(f"  J_TDC(TDDFT) = {J:.2f} cm^-1     J_PDA(TDDFT) = {pda:.2f} cm^-1     TDC/PDA = {abs(J/pda):.2f}x")
print(f"  (definitive STEOM on same geometry: J_TDC=30.5, J_PDA=24.8; STEOM/TDDFT = {30.453973388943076/abs(J):.2f}x)")

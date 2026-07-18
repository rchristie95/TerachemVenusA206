#!/usr/bin/env python3
"""Rigorous STEOM transition density: rho_0n = sum_k sigma_k phi_hole_k phi_part_k
over the 4 significant NTO pairs (98.6% of S1), NO dipole-forcing. Check it reproduces
the ORCA dipole on its own, then compute J on the identical coupling pipeline."""
import numpy as np, sys, time
from numba import njit, prange, set_num_threads
set_num_threads(16); sys.path.insert(0,"/home/robson/PetaChem")
import coupling_core as cc
BOHR=0.52917721067; HC=219474.6314; EPS=1.78; A2B=cc.ANGSTROM_TO_BOHR
MU_STEOM=np.array([1.08027,-1.98184,3.01365])  # au, ORCA SVPD state-1
D="/home/robson/PetaChem/neo_model/orca_steom/"

def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()
        t=f.readline().split(); na=abs(int(t[0])); o=np.array(list(map(float,t[1:4])))
        n=[];v=[]
        for _ in range(3):
            r=f.readline().split(); n.append(int(r[0])); v.append(list(map(float,r[1:4])))
        n=np.array(n); v=np.array(v)
        for _ in range(na): f.readline()
        if int(t[0])<0: f.readline()
        data=np.fromstring(" ".join(f.read().split()),sep=" ")
    return o,n,v,data.reshape(n)

# NTO pairs: (hole, particle, occupation n_k)
pairs=[(92,93,0.94764284),(91,94,0.02484332),(90,95,0.00866944),(89,96,0.00518175)]
o,n,v,_=read_cube(D+"steom_phenol_svpd.s1.mo92a.cube")
dvol=abs(np.linalg.det(v))
ii,jj,kk=np.meshgrid(np.arange(n[0]),np.arange(n[1]),np.arange(n[2]),indexing="ij")
pts_bohr=(o[None,None,None,:]+ii[...,None]*v[0]+jj[...,None]*v[1]+kk[...,None]*v[2]).reshape(-1,3)
muhat=MU_STEOM/np.linalg.norm(MU_STEOM)

rho_tot=np.zeros(n[0]*n[1]*n[2])
print("per-pair:  sigma   |d_k|(au)  sign  proj(d_k.muhat)")
for h,p,nk in pairs:
    _,_,_,ph=read_cube(D+f"steom_phenol_svpd.s1.mo{h}a.cube")
    _,_,_,pp=read_cube(D+f"steom_phenol_svpd.s1.mo{p}a.cube")
    rho_k=(ph*pp).reshape(-1)
    q_k=rho_k*dvol
    d_k=-(pts_bohr*q_k[:,None]).sum(0)          # raw dipole of unit product
    sig=np.sqrt(nk); s=np.sign(np.dot(d_k,MU_STEOM))
    print(f"  {h}->{p}: {sig:6.4f}  {np.linalg.norm(d_k):7.3f}   {s:+.0f}   {np.dot(d_k,muhat):+7.3f}")
    rho_tot += s*sig*rho_k

q_tot=rho_tot*dvol
mu_tot=-(pts_bohr*q_tot[:,None]).sum(0)
print(f"\nSUMMED multipair dipole  mu={mu_tot.round(3)} au |mu|={np.linalg.norm(mu_tot):.3f}")
print(f"ORCA target              mu={MU_STEOM.round(3)} au |mu|={np.linalg.norm(MU_STEOM):.3f}")
implied_scale=np.dot(MU_STEOM,mu_tot)/np.dot(mu_tot,mu_tot)
print(f"implied residual scale = {implied_scale:.4f} (want ~1.0; old single-pair build needed 1.268)")

# save (no forcing) + threshold like the original
pts_ang=pts_bohr*BOHR
thr=1e-6*np.abs(q_tot).max(); keep=np.abs(q_tot)>thr
np.savez(D+"steom_transdens_multipair.npz",pts_ang=pts_ang[keep],q=q_tot[keep],mu_au=mu_tot)
print(f"saved {keep.sum()} pts; sum q={q_tot[keep].sum():+.2e}")

# ---- coupling on the identical pipeline ----
@njit(parallel=True,fastmath=True,cache=True)
def K(p1,q1,p2,q2):
    n1=p1.shape[0];n2=p2.shape[0];t=0.0
    for i in prange(n1):
        xi=p1[i,0];yi=p1[i,1];zi=p1[i,2];qi=q1[i];a=0.0
        for j in range(n2):
            dx=xi-p2[j,0];dy=yi-p2[j,1];dz=zi-p2[j,2];r=(dx*dx+dy*dy+dz*dz)**0.5
            if r<0.1:r=0.1
            a+=qi*q2[j]/r
        t+=a
    return t
mA,mB,aA,aB,err=cc.get_super_matrices_with_pymol("/home/robson/PetaChem/tc_simple_anionic/monomer_relaxed.pdb","/home/robson/PetaChem/venus_dimer.pdb")
pk=np.ascontiguousarray(pts_ang[keep],np.float64); qk=np.ascontiguousarray(q_tot[keep],np.float64)
pA=np.ascontiguousarray(cc.apply_pymol_matrix(pk,mA),np.float64)
pB=np.ascontiguousarray(cc.apply_pymol_matrix(pk,mB),np.float64)
J=K(pA,qk,pB,qk)*A2B/EPS*HC
print(f"\n=== STEOM coupling, RIGOROUS multipair density ===")
print(f"  J_TDC(STEOM,multipair) = {J:.2f} cm^-1   2|J| = {2*abs(J):.1f}")
print(f"  (old single-pair forced = 154.82 ; matched-geom TDDFT full = 117.8 ; exp 2|J|~130)")

import numpy as np
D="/home/robson/PetaChem/neo_model/orca_steom/"
MU=np.array([1.08027,-1.98184,3.01365])  # ORCA spectroscopic, |mu|=3.765
d=np.load(D+"steom_transdens_multipair.npz")
pts,q,mu_tot=d["pts_ang"],d["q"],d["mu_au"]
scale=np.dot(MU,mu_tot)/np.dot(mu_tot,mu_tot)        # TrEsp projection to spectroscopic dipole
q2=q*scale
mu2=-(0)  # recompute below
# new dipole (charges already include dvol; mu = -sum r*q in bohr)
BOHR=0.52917721067
mu_new=-(pts/BOHR*q2[:,None]).sum(0)
np.savez(D+"steom_transdens_specnorm.npz",pts_ang=pts,q=q2,mu_au=mu_new)
Junf=102.12  # J of the unforced multipair density (numba, eps=1.78)
print(f"TrEsp scale to spectroscopic dipole = {scale:.4f}")
print(f"new |mu| = {np.linalg.norm(mu_new):.3f} au (target 3.765)")
print(f"J(specnorm) = {Junf*scale**2:.1f} cm^-1  (= J_unforced 102.12 * scale^2)")
print(f"saved -> steom_transdens_specnorm.npz  ({len(q2)} pts, sum q={q2.sum():+.2e})")

#!/usr/bin/env python3
"""Build the STEOM-CCSD bright-state transition density from the NTO hole/particle cubes,
verify it against the STEOM transition dipole, and save it as points+charges (au) for the
coupling pipeline. rho_0n(r) ~ sqrt(lambda) * phi_hole(r) * phi_particle(r) (dominant NTO pair).
We fix the absolute normalization by matching the transition dipole to the ORCA value."""
import numpy as np

BOHR = 0.529177210903  # Angstrom per bohr
MU_STEOM = np.array([1.08027, -1.98184, 3.01365])  # au, from SVPD spectrum (state 1)

def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()                      # 2 comment lines
        toks = f.readline().split(); natom = int(toks[0])
        origin = np.array(list(map(float, toks[1:4])))  # bohr
        n = []; vec = []
        for _ in range(3):
            t = f.readline().split(); n.append(int(t[0])); vec.append(list(map(float, t[1:4])))
        n = np.array(n); vec = np.array(vec)            # bohr
        na = abs(natom)
        for _ in range(na): f.readline()                # atom lines
        if natom < 0: f.readline()                      # MO-cube: orbital-index line
        data = np.fromstring(" ".join(f.read().split()), sep=" ")
    return origin, n, vec, data.reshape(n)

o1, n1, v1, hole = read_cube("neo_model/orca_steom/steom_phenol_svpd.s1.mo92a.cube")
o2, n2, v2, part = read_cube("neo_model/orca_steom/steom_phenol_svpd.s1.mo93a.cube")
assert np.allclose(o1, o2) and (n1 == n2).all() and np.allclose(v1, v2), "grid mismatch"

rho = hole * part                                       # transition density (unnormalized), shape n
dvol = abs(np.linalg.det(v1))                           # voxel volume, bohr^3

# grid point coordinates (bohr)
ii, jj, kk = np.meshgrid(np.arange(n1[0]), np.arange(n1[1]), np.arange(n1[2]), indexing="ij")
pts_bohr = (o1[None,None,None,:] + ii[...,None]*v1[0] + jj[...,None]*v1[1] + kk[...,None]*v1[2])

# transition dipole of the raw NTO product:  mu = -integral r * rho dr  (electron charge -1)
q = rho * dvol                                          # "charge" per voxel
mu_raw = -(pts_bohr.reshape(-1,3) * q.reshape(-1,1)).sum(axis=0)   # au
# normalization factor to match the ORCA transition dipole (sign+magnitude via projection)
scale = np.dot(MU_STEOM, mu_raw) / np.dot(mu_raw, mu_raw)
mu_scaled = mu_raw * scale

print(f"grid {tuple(n1)}  voxel {dvol:.4f} bohr^3  origin {o1.round(2)}")
print(f"raw NTO-product dipole   mu = {mu_raw.round(3)} au, |mu|={np.linalg.norm(mu_raw):.3f}")
print(f"STEOM (ORCA) dipole      mu = {MU_STEOM.round(3)} au, |mu|={np.linalg.norm(MU_STEOM):.3f}")
print(f"scaled NTO-product       mu = {mu_scaled.round(3)} au, |mu|={np.linalg.norm(mu_scaled):.3f}")
print(f"normalization scale = {scale:.4f}  | direction cos(theta)="
      f"{np.dot(mu_raw,MU_STEOM)/(np.linalg.norm(mu_raw)*np.linalg.norm(MU_STEOM)):.4f}")

# save normalized transition density as points(Angstrom)+charges(au) above a threshold
q_norm = (q * scale).reshape(-1)
pts_ang = pts_bohr.reshape(-1,3) * BOHR
thr = 1e-6 * np.abs(q_norm).max()
keep = np.abs(q_norm) > thr
np.savez("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.npz",
         pts_ang=pts_ang[keep], q=q_norm[keep], mu_au=mu_scaled)
print(f"saved {keep.sum()} grid points (>|{thr:.2e}|) to steom_transdens.npz; "
      f"sum q = {q_norm[keep].sum():+.3e} (should be ~0)")

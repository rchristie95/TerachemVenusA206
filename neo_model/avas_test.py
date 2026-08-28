#!/usr/bin/env python3
"""Validate the pi active space for CASSCF: reorient the (planar) chromophore so
the molecular-plane normal = z, then AVAS on C/N/O 2pz picks the pi/pi* space."""
import numpy as np
from pyscf import gto, scf
from pyscf.mcscf import avas

L = open('cr2_chromophore.xyz').readlines(); n = int(L[0].split()[0])
syms = [L[2+i].split()[0] for i in range(n)]
xyz = np.array([[float(c) for c in L[2+i].split()[1:4]] for i in range(n)])
heavy = [i for i, s in enumerate(syms) if s != 'H']
P = xyz[heavy]; c = P.mean(0)
_, _, Vt = np.linalg.svd(P - c)          # Vt[2] = smallest-variance axis = plane normal
xyz_rot = (xyz - c) @ Vt.T               # now the normal lies along z  -> pi = 2pz
rms_oop = np.sqrt((xyz_rot[heavy, 2]**2).mean())
print(f"heavy-atom out-of-plane RMS after reorient: {rms_oop:.3f} A (planar if small)")

mol = gto.M(atom=[(syms[i], *xyz_rot[i]) for i in range(n)],
            basis='def2-svp', charge=-1, spin=0, verbose=2)
mf = scf.RHF(mol).density_fit(); mf.kernel()
print(f"RHF converged={mf.converged}  E={mf.e_tot:.6f}")

for thr in (0.20, 0.15, 0.10, 0.30):
    ncas, nelec, _ = avas.avas(mf, ['C 2pz', 'N 2pz', 'O 2pz'],
                               threshold=thr, canonicalize=False, verbose=0)
    print(f"AVAS threshold={thr:>4}:  ncas={ncas:2d}  nelecas={nelec}")

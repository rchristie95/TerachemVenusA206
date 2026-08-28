#!/usr/bin/env python3
"""Pick a (12,12) frontier-pi active space by ranking MOs on their Mulliken
population in the heavy-atom 2pz AOs (reoriented so pi = 2pz)."""
import numpy as np
from pyscf import gto, scf

L = open('cr2_chromophore.xyz').readlines(); n = int(L[0].split()[0])
syms = [L[2+i].split()[0] for i in range(n)]
xyz = np.array([[float(c) for c in L[2+i].split()[1:4]] for i in range(n)])
heavy = [i for i, s in enumerate(syms) if s != 'H']
P = xyz[heavy]; c = P.mean(0)
_, _, Vt = np.linalg.svd(P - c)
xyz_rot = (xyz - c) @ Vt.T
mol = gto.M(atom=[(syms[i], *xyz_rot[i]) for i in range(n)],
            basis='def2-svp', charge=-1, spin=0, verbose=2)
mf = scf.RHF(mol).density_fit(); mf.kernel()

pz = mol.search_ao_label(['C 2pz', 'N 2pz', 'O 2pz'])
S = mf.get_ovlp(); C = mf.mo_coeff; SC = S @ C
pi = np.einsum('mi,mi->i', C[pz, :], SC[pz, :])    # pi-population per MO (0..1)
nocc = mol.nelectron // 2
occ = np.arange(nocc); vir = np.arange(nocc, mol.nao)
act_occ = sorted(occ[np.argsort(-pi[occ])][:6])
act_vir = sorted(vir[np.argsort(-pi[vir])][:6])
active = act_occ + act_vir
print(f"HOMO = MO {nocc-1};  selected (12,12) frontier-pi active space:")
for i in active:
    print(f"  MO{i:3d} {'occ' if i < nocc else 'vir'}  e={mf.mo_energy[i]*27.2114:7.2f} eV  pi-pop={pi[i]:.2f}")
print("contiguous around gap?  occ:", act_occ, " vir:", act_vir)
print("min pi-pop in active:", round(min(pi[active]), 2),
      "| max pi-pop OUTSIDE active (gap-adjacent leak):",
      round(max(pi[[k for k in range(mol.nao) if k not in active]]), 2))
np.save('pi_active_idx.npy', np.array(active))

#!/usr/bin/env python3
"""ADC(2): chromophore + Tyr203-phenol (44 atoms) in the LIGHTWEIGHT protein field
(~2.4k pocket charges, not the 30k solvent box). MM-min geometry (deterministic).
Frozen 1s cores + top virtuals for memory; nroots=4 (bright state splits)."""
import time, numpy as np
from pyscf import gto, scf, adc, qmmm
t0 = time.time()

mol = gto.M(atom='qm_cthrp.xyz', basis='def2-svp', charge=-1, spin=0,
            max_memory=40000, verbose=4)
D = open('mm_light_cthrp.dat').readlines(); nmm = int(D[0].split()[0])
mm = np.array([[float(x) for x in D[2+i].split()] for i in range(nmm)])
print(f"[{time.time()-t0:.0f}s] {mol.nao} AOs, {nmm} lightweight MM charges (net {mm[:,0].sum():+.3f})", flush=True)

mf = qmmm.mm_charge(scf.RHF(mol).density_fit(), mm[:, 1:4], mm[:, 0], unit='Angstrom')
mf.conv_tol = 1e-9; mf.kernel()
if not mf.converged:
    mf = mf.newton(); mf.kernel()
nocc, nao = mol.nelectron // 2, mol.nao
print(f"[{time.time()-t0:.0f}s] RHF+MM E={mf.e_tot:.6f} conv={mf.converged} nao={nao} nocc={nocc}", flush=True)

ncore, nfvir = 26, 140
frozen = list(range(ncore)) + list(range(nao - nfvir, nao))
print(f"active occ={nocc-ncore} active vir={nao-nocc-nfvir}", flush=True)
myadc = adc.ADC(mf, frozen=frozen); myadc.method = "adc(2)"; myadc.method_type = "ee"
e = np.atleast_1d(myadc.kernel(nroots=4)[0])
print(f"[{time.time()-t0:.0f}s] ADC(2) done", flush=True)
myadc.analyze()

ev = 27.211386; nm = 1239.84193
print("\n===== ADC(2) chromophore+Tyr203-phenol + lightweight protein field =====")
for i, en in enumerate(e):
    print(f"  root {i}: {en*ev:.3f} eV  {nm/(en*ev):.1f} nm")
print("  bright = max-osc root in analyze() above (check 2p2h vs CT)")
print("  ref: chromophore-alone ADC(2) 476 nm ; TDDFT composite ~505 ; exp ~515 nm")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

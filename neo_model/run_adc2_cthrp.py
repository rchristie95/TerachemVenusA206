#!/usr/bin/env python3
"""ADC(2) confirmation of the Tyr203 pi-stacking red-shift. QM = chromophore +
Tyr203 phenol ring (44 atoms), gas-phase (MM effect ~0), MM-min geometry.
Frozen 1s cores + top virtuals to fit 62 GB. nroots=5 (the bright state splits).
Baselines: chromophore-alone ADC(2) 476 nm ; TDDFT composite ~505 nm ; exp 515 nm."""
import time, numpy as np
from pyscf import gto, scf, adc
t0 = time.time()

mol = gto.M(atom='qm_cthrp.xyz', basis='def2-svp', charge=-1, spin=0,
            max_memory=40000, verbose=4)
mf = scf.RHF(mol).density_fit(); mf.conv_tol = 1e-9
mf.kernel()
if not mf.converged:
    mf = mf.newton(); mf.kernel()
nocc = mol.nelectron // 2; nao = mol.nao
print(f"[{time.time()-t0:.0f}s] RHF E={mf.e_tot:.6f} conv={mf.converged} nao={nao} nocc={nocc}", flush=True)

ncore, nfvir = 26, 100                      # freeze C/N/O 1s (19C? actually 26) + top 100 virtuals
frozen = list(range(ncore)) + list(range(nao - nfvir, nao))
print(f"frozen core={ncore} frozen vir={nfvir} -> active occ={nocc-ncore} active vir={nao-nocc-nfvir}", flush=True)

myadc = adc.ADC(mf, frozen=frozen); myadc.method = "adc(2)"; myadc.method_type = "ee"
e = np.atleast_1d(myadc.kernel(nroots=5)[0])
print(f"[{time.time()-t0:.0f}s] ADC(2) done", flush=True)
myadc.analyze()

ev = 27.211386; nm = 1239.84193
print("\n========= ADC(2) chromophore + Tyr203-phenol (gas, MM-min geom) =========")
for i, en in enumerate(e):
    print(f"  root {i}: {en*ev:.3f} eV  {nm/(en*ev):.1f} nm")
print("  bright = max-osc root in analyze() above (check 2p2h vs CT character)")
print("  ref: chromophore-alone ADC(2) 476 nm ; exp in-protein ~515 nm")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

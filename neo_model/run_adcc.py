#!/usr/bin/env python3
"""ADC(2) via adcc (optimized, fully-threaded ADC engine) on the pyscf RHF(+MM).
adcc fixes the parallel/memory inefficiency of pyscf's native ADC. Full virtuals
on the 44-atom system is still ~o^2 v^2 large, so frozen_virtual is tunable; adcc's
efficiency lets us keep MORE virtuals (less blue bias) in the same RAM.
Usage: run_adcc.py <qm.xyz> <mm.dat|gas> <n_frozen_core> <n_frozen_virtual> <nstates>
Env: adcc_env."""
import sys, time, numpy as np
import adcc
from pyscf import gto, scf, qmmm

t0 = time.time()
geom = sys.argv[1]
mmf = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != 'gas' else None
nfc = int(sys.argv[3]) if len(sys.argv) > 3 else 26
nfv = int(sys.argv[4]) if len(sys.argv) > 4 else 0
nst = int(sys.argv[5]) if len(sys.argv) > 5 else 4

mol = gto.M(atom=geom, basis='def2-svp', charge=-1, spin=0, max_memory=45000, verbose=4)
mf = scf.RHF(mol).density_fit()
if mmf:
    D = open(mmf).readlines(); nmm = int(D[0].split()[0])
    mm = np.array([[float(x) for x in D[2+i].split()] for i in range(nmm)])
    mf = qmmm.mm_charge(mf, mm[:, 1:4], mm[:, 0], unit='Angstrom')
mf.conv_tol = 1e-10
mf.kernel()
print(f"[{time.time()-t0:.0f}s] RHF E={mf.e_tot:.6f} conv={mf.converged} nao={mol.nao}", flush=True)

adcc.set_n_threads(32)
kw = dict(n_singlets=nst, conv_tol=1e-5)
if nfc: kw['frozen_core'] = nfc
if nfv: kw['frozen_virtual'] = nfv
print(f"adcc ADC(2): frozen_core={nfc} frozen_virtual={nfv} n_singlets={nst} threads={adcc.get_n_threads()}", flush=True)
state = adcc.adc2(mf, **kw)
print(f"[{time.time()-t0:.0f}s] ADC(2) done", flush=True)
print(state.describe())

ev = 27.211386; nm = 1239.84193
print("\n===== adcc ADC(2) =====")
for i in range(len(state.excitation_energy)):
    e = float(state.excitation_energy[i]); f = float(state.oscillator_strength[i])
    print(f"  S{i+1}: {e*ev:.3f} eV  {nm/(e*ev):.1f} nm  f={f:.3f}")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

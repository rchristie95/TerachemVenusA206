#!/usr/bin/env python3
"""ADC(2) on the anionic chromophore in the protein MM point-charge field (QM/MM).
Same QM geometry/basis as the gas-phase run (cr2_chromophore.xyz, def2-SVP) so the
shift vs gas-phase (bright 2.60 eV / 476 nm) is the clean protein electrostatic effect.
His148 is in the MM field. Env: pyscfneo + LD_LIBRARY_PATH anaconda3/lib."""
import time, numpy as np
from pyscf import gto, scf, adc, qmmm
t0 = time.time()

mol = gto.M(atom='cr2_chromophore.xyz', basis='def2-svp', charge=-1, spin=0,
            max_memory=38000, verbose=4)
D = open('mm_field.dat').readlines(); nmm = int(D[0].split()[0])
mm = np.array([[float(x) for x in D[2+i].split()] for i in range(nmm)])
print(f"[{time.time()-t0:.0f}s] {mol.nao} AOs, {nmm} MM charges (net {mm[:,0].sum():+.3f})", flush=True)

mf = qmmm.mm_charge(scf.RHF(mol).density_fit(), mm[:, 1:4], mm[:, 0], unit='Angstrom')
mf.conv_tol = 1e-9
mf.kernel()
if not mf.converged:
    mf = mf.newton(); mf.kernel()
print(f"[{time.time()-t0:.0f}s] RHF+MM E={mf.e_tot:.6f} conv={mf.converged}", flush=True)

myadc = adc.ADC(mf, frozen=19); myadc.method = "adc(2)"; myadc.method_type = "ee"
e = np.atleast_1d(myadc.kernel(nroots=3)[0])
print(f"[{time.time()-t0:.0f}s] ADC(2)+MM done", flush=True)
myadc.analyze()                       # prints per-root oscillator strengths -> bright state

ev = 27.211386; nm = 1239.84193
print("\n================ ADC(2)+MM (protein field) ================")
for i, en in enumerate(e):
    print(f"  root {i}: {en*ev:.3f} eV  {nm/(en*ev):.1f} nm")
print("  bright = max-osc root in analyze() above")
print("  GAS-PHASE reference: bright 2.60 eV / 476 nm ; exp in-protein ~2.41 eV / 515 nm")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

#!/usr/bin/env python3
"""RI-ADC(2) vs TDDFT on the bare anionic CR2 chromophore (His stripped).
The ADC(2) - TDDFT difference at the SAME geometry/basis isolates the cyanine
(GFP-chromophore) error that makes LR-TDDFT overestimate S1 by ~0.7-0.9 eV.
Env: pyscfneo ; LD_LIBRARY_PATH must include anaconda3/lib (MKL)."""
import time, numpy as np
from pyscf import gto, scf, adc, dft, tddft

t0 = time.time()
mol = gto.M(atom='cr2_chromophore.xyz', basis='def2-svp', charge=-1, spin=0,
            max_memory=38000, verbose=4)   # cap below 62 GB (chrome etc. share RAM)
print(f"[{time.time()-t0:.0f}s] mol: {mol.natm} atoms, {mol.nao} basis fns, {mol.nelectron} electrons", flush=True)

# --- density-fitted RHF reference (=> RI-ADC) ---
mf = scf.RHF(mol).density_fit()
mf.conv_tol = 1e-9
mf.kernel()
if not mf.converged:                       # anion: fall back to second-order SCF
    print("[!] DIIS not converged; retrying with Newton (SOSCF)", flush=True)
    mf = mf.newton(); mf.kernel()
print(f"[{time.time()-t0:.0f}s] RHF E = {mf.e_tot:.8f}  converged={mf.converged}", flush=True)

# --- RI-ADC(2) excitation energies (frozen core via constructor!) ---
myadc = adc.ADC(mf, frozen=19)             # freeze C/N/O 1s cores (13+3+3); o: 68->49
myadc.method = "adc(2)"
myadc.method_type = "ee"
e_adc = np.atleast_1d(myadc.kernel(nroots=3)[0])   # 3 roots keeps the doubles-vector RAM bounded
print(f"[{time.time()-t0:.0f}s] ADC(2) done", flush=True)
try:
    myadc.analyze()                        # prints spectroscopic factors if available
except Exception as exc:
    print(f"  (adc analyze skipped: {exc})", flush=True)

# --- matched TDDFT (CAM-B3LYP) at the same geometry/basis ---
mfd = dft.RKS(mol, xc='camb3lyp').density_fit()
mfd.conv_tol = 1e-9
mfd.kernel()
td = tddft.TDDFT(mfd); td.nstates = 5
td.kernel()
f_td = td.oscillator_strength()
print(f"[{time.time()-t0:.0f}s] TDDFT done", flush=True)

ev = 27.211386; nm = 1239.84193
print("\n================ RESULTS (def2-SVP, bare CR2 anion) ================")
print(f"{'state':>5} {'ADC(2) eV':>10} {'ADC(2) nm':>10}   {'TDDFT eV':>9} {'TDDFT nm':>9} {'f_TDDFT':>8}")
for i in range(len(e_adc)):
    a_ev = e_adc[i]*ev; a_nm = nm/a_ev
    t_ev = td.e[i]*ev;  t_nm = nm/t_ev
    print(f"{i+1:>5} {a_ev:>10.3f} {a_nm:>10.1f}   {t_ev:>9.3f} {t_nm:>9.1f} {f_td[i]:>8.3f}")
# bright TDDFT state = max f; its cyanine shift vs the nearest ADC(2) state
ib = int(np.argmax(f_td))
print(f"\nBright TDDFT state: S{ib+1} = {td.e[ib]*ev:.3f} eV ({nm/(td.e[ib]*ev):.1f} nm), f={f_td[ib]:.3f}")
print(f"ADC(2) S1        : {e_adc[0]*ev:.3f} eV ({nm/(e_adc[0]*ev):.1f} nm)")
print(f"Cyanine shift (TDDFT_bright - ADC2_S1) = {(td.e[ib]-e_adc[0])*ev:+.3f} eV")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

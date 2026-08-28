#!/usr/bin/env python3
"""Fast TDDFT (CAM-B3LYP/def2-SVP) screen on the MM-min geometry to (1) measure the
Tyr203 pi-stacking shift and the MM-field shift, and (2) flag chromophore<->Tyr203
charge-transfer contamination (a new low-lying DARK state in chr+Tyr absent in chr).
TDDFT absolute is cyanine-biased (bluer); read the SHIFTS, then add to ADC(2) 476 nm."""
import os, numpy as np
from pyscf import gto, scf, dft, tddft, qmmm
os.chdir('/home/robson/PetaChem/neo_model')
ev = 27.211386; nm = 1239.84193

def run(tag, qmxyz, mmfile):
    mol = gto.M(atom=qmxyz, basis='def2-svp', charge=-1, spin=0, max_memory=30000, verbose=0)
    rks = dft.RKS(mol, xc='camb3lyp').density_fit()
    if mmfile:
        D = open(mmfile).readlines(); nmm = int(D[0].split()[0])
        mm = np.array([[float(x) for x in D[2+i].split()] for i in range(nmm)])
        rks = qmmm.mm_charge(rks, mm[:, 1:4], mm[:, 0], unit='Angstrom')
    rks.kernel()
    td = tddft.TDDFT(rks); td.nstates = 6; td.kernel()
    f = td.oscillator_strength(); ib = int(np.argmax(f))
    print(f"\n[{tag}]  bright S{ib+1}: {td.e[ib]*ev:.3f} eV  {nm/(td.e[ib]*ev):.1f} nm  f={f[ib]:.3f}", flush=True)
    for i in range(6):
        mk = " <-BRIGHT" if i == ib else (" dark" if f[i] < 0.05 else "")
        print(f"    S{i+1}: {td.e[i]*ev:.3f} eV {nm/(td.e[i]*ev):7.1f} nm  f={f[i]:.3f}{mk}", flush=True)
    return td.e[ib]*ev, f[ib]

print("=== TDDFT CAM-B3LYP/def2-SVP screen (MM-min geometry) ===", flush=True)
r = {}
r['chr_gas']  = run('chr      gas', 'qm_chr.xyz',  None)
r['chr_mm']   = run('chr      +MM', 'qm_chr.xyz',  'mm_chr.dat')
r['cthr_gas'] = run('chr+Tyr  gas', 'qm_cthr.xyz', None)
r['cthr_mm']  = run('chr+Tyr  +MM', 'qm_cthr.xyz', 'mm_cthr.dat')
print("\n================ SUMMARY (bright state) ================")
for k, (e, f) in r.items():
    print(f"  {k:9s}: {e:.3f} eV  {nm/e:6.1f} nm  f={f:.2f}")
print(f"\n  Tyr203 pi-stack shift (gas):   {nm/r['chr_gas'][0]:.0f} -> {nm/r['cthr_gas'][0]:.0f} nm  "
      f"({(r['cthr_gas'][0]-r['chr_gas'][0]):+.3f} eV)")
print(f"  MM-field shift (chromophore):  {nm/r['chr_gas'][0]:.0f} -> {nm/r['chr_mm'][0]:.0f} nm  "
      f"({(r['chr_mm'][0]-r['chr_gas'][0]):+.3f} eV)")
print(f"  FULL chr+Tyr+MM: {nm/r['cthr_mm'][0]:.0f} nm  (TDDFT cyanine-biased; exp 515)")
print("ALL DONE")

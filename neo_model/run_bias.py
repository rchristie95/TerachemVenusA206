#!/usr/bin/env python3
"""Bias-correction for the truncation blue-shift. Chromophore-alone (31 atoms) fits
FULL virtuals, so measure bias = E_trunc - E_full at the SAME 62% virtual fraction
used for chr+Tyr (357->221 virt), then subtract it from the chr+Tyr 451 nm number.
Truncation bias is a correlation effect -> ~field-independent, so run gas-phase."""
import adcc, numpy as np
from pyscf import gto, scf
ev = 27.211386; nm = 1239.84193

def run(geom, nfc, nfv, nst=2):
    mol = gto.M(atom=geom, basis='def2-svp', charge=-1, spin=0, max_memory=30000, verbose=0)
    mf = scf.RHF(mol).density_fit(); mf.conv_tol = 1e-10; mf.kernel()
    adcc.set_n_threads(32)
    # cap the Jacobi-Davidson subspace so it restarts instead of hoarding doubles vectors
    kw = dict(n_singlets=nst, conv_tol=1e-5, frozen_core=nfc, max_subspace=2*nst+4)
    if nfv: kw['frozen_virtual'] = nfv
    st = adcc.adc2(mf, **kw)
    e = np.array([float(x) for x in st.excitation_energy])
    f = np.array([float(x) for x in st.oscillator_strength])
    ib = int(np.argmax(f))
    return e[ib], f[ib], mol.nao - mol.nelectron // 2

print("=== chr-alone (C-O 1.263) FULL virtuals ===", flush=True)
e_full, f_full, nvir = run("qm_chr_relaxed.xyz", 19, 0, 2)
print(f"  bright: {e_full*ev:.3f} eV  {nm/(e_full*ev):.1f} nm  f={f_full:.2f}  (nvir={nvir})", flush=True)

nfv = int(round(nvir * (1 - 221/357)))      # match chr+Tyr's 62%-kept fraction
print(f"=== chr-alone 62%-virt (freeze {nfv} of {nvir}) ===", flush=True)
e_tr, f_tr, _ = run("qm_chr_relaxed.xyz", 19, nfv, 2)
print(f"  bright: {e_tr*ev:.3f} eV  {nm/(e_tr*ev):.1f} nm  f={f_tr:.2f}", flush=True)

bias = e_tr - e_full
e_cthr = 2.749 / ev   # chr+Tyr relaxed, pyscf 221-virt = 451 nm = 2.749 eV
corr = e_cthr - bias
print("\n================ BIAS-CORRECTED IN-PROTEIN ESTIMATE ================")
print(f"  truncation bias (62%): {bias*ev*1000:+.0f} meV (blue)")
print(f"  chr+Tyr raw (221-virt): 2.749 eV / 451.0 nm")
print(f"  chr+Tyr BIAS-CORRECTED: {corr*ev:.3f} eV / {nm/(corr*ev):.1f} nm")
print(f"  experiment (in-protein): 2.41 eV / 515 nm")
print("ALL DONE")

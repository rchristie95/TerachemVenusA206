#!/usr/bin/env python3
"""SA-CASSCF(12,10)+NEVPT2 cross-check of the cyanine S1 (bare anionic chromophore).
Active = 6 valence pi + 4 clean valence pi*. nroots=5 (bright singlet sits above S0/T1/dark-S1).
NEVPT2 needs a SEPARATE multi-root CASCI on the SA-CASSCF orbitals (pyscf requirement).
CAVEAT: pi-poor virtual space on a non-planar geometry -> a rough cross-check.
Compare: ADC(2) 2.60 eV/476 nm ; CAM-B3LYP TDDFT 3.18 eV/390 nm ; exp ~2.5 eV/~480 nm."""
import time, numpy as np
from pyscf import gto, scf, mcscf, mrpt
from pyscf.fci import direct_spin1
t0 = time.time()

L = open('cr2_chromophore.xyz').readlines(); n = int(L[0].split()[0])
syms = [L[2+i].split()[0] for i in range(n)]
xyz = np.array([[float(c) for c in L[2+i].split()[1:4]] for i in range(n)])
heavy = [i for i, s in enumerate(syms) if s != 'H']
P = xyz[heavy]; c = P.mean(0); _, _, Vt = np.linalg.svd(P - c)
mol = gto.M(atom=[(syms[i], *((xyz[i]-c) @ Vt.T)) for i in range(n)],
            basis='def2-svp', charge=-1, spin=0, max_memory=34000, verbose=4)
mf = scf.RHF(mol).density_fit(); mf.kernel()
print(f"[{time.time()-t0:.0f}s] RHF E={mf.e_tot:.6f} conv={mf.converged}", flush=True)

active = [56, 60, 62, 65, 66, 67, 68, 72, 74, 77]      # 0-based: 6 occ pi + 4 pi*
ncas, nelecas, nst = 10, 12, 5
mc = mcscf.CASSCF(mf, ncas, nelecas)
mc.state_average_([1.0/nst]*nst)
mc.fcisolver.nroots = nst
mc.fix_spin_(ss=0, shift=0.5)
mc.max_cycle_macro = 100
mc.chkfile = 'casscf_chk.h5'                            # persist converged orbitals
mc.kernel(mc.sort_mo(active, base=0))
print(f"[{time.time()-t0:.0f}s] SA-CASSCF done; state energies (Ha): {np.array(mc.e_states)}", flush=True)

# separate multi-root CASCI on the SA-CASSCF orbitals (required for NEVPT2)
ci = mcscf.CASCI(mf, ncas, nelecas)
ci.fcisolver.nroots = nst
ci.kernel(mc.mo_coeff)
ev = 27.21138; e = np.array(ci.e_tot); e0 = e[0]

# oscillator strengths -> identify the bright singlet (all 5 states are S^2=0)
mo_cas = mc.mo_coeff[:, ci.ncore:ci.ncore+ncas]
dip_ao = mol.intor('int1e_r')
osc = {}
print("  state  dE_CASSCF(eV)    f")
for i in range(1, nst):
    tdm = direct_spin1.trans_rdm1(ci.ci[0], ci.ci[i], ncas, nelecas)
    tdip = -np.einsum('xij,ji->x', dip_ao, mo_cas @ tdm @ mo_cas.T)
    osc[i] = 2./3. * (e[i]-e0) * float(tdip @ tdip)
    print(f"   {i}      {(e[i]-e0)*ev:7.3f}     {osc[i]:6.3f}", flush=True)
ib = max(osc, key=osc.get)
print(f"Bright singlet = state {ib}: CASSCF {(e[ib]-e0)*ev:.3f} eV "
      f"({1239.84193/((e[ib]-e0)*ev):.1f} nm)", flush=True)

# NEVPT2 dynamic correlation on ground + bright
c_gs = mrpt.NEVPT(ci, root=0).kernel()
c_br = mrpt.NEVPT(ci, root=ib).kernel()
dE = ((e[ib]+c_br) - (e0+c_gs)) * ev
print("\n================ CASSCF / NEVPT2 RESULT (bright singlet) ================")
print(f"  CASSCF : {(e[ib]-e0)*ev:.3f} eV ({1239.84193/((e[ib]-e0)*ev):.1f} nm)")
print(f"  NEVPT2 : {dE:.3f} eV ({1239.84193/dE:.1f} nm)")
print(f"  cross-check vs ADC(2) 2.60 eV/476 nm, TDDFT 3.18 eV/390 nm, exp ~2.5 eV/~480 nm")
print(f"[{time.time()-t0:.0f}s] ALL DONE")

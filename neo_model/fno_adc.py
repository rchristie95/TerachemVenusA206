#!/usr/bin/env python3
"""FNO-ADC(2): build frozen natural orbitals from the MP2 virtual density, keep the
correlation-important virtuals (high NO occupation), semicanonicalize, run EE-ADC(2).
Principled virtual truncation -> removes the blue bias of arbitrary energy-freezing,
fits memory, ~faster. Returns excitation energies + the adc object (for analyze)."""
import numpy as np
from pyscf import mp, adc


def make_fno(mf, ncore, n_keep):
    """Return (mo_coeff, mo_energy, frozen_list) with virtuals as semicanonical FNOs;
    frozen = core + the (nvir - n_keep) lowest-occupation FNOs."""
    nmo = mf.mo_coeff.shape[1]
    nocc = int((mf.mo_occ > 0).sum())
    pt = mp.MP2(mf)
    pt.frozen = ncore
    pt.kernel()
    dm = pt.make_rdm1()                       # full-MO MP2 1-RDM
    dvv = dm[nocc:, nocc:]                     # virtual-virtual block
    no_occ, V = np.linalg.eigh(dvv)
    order = np.argsort(no_occ)[::-1]          # descending occupation
    V = V[:, order]; no_occ = no_occ[order]
    vir = mf.mo_coeff[:, nocc:] @ V           # virtuals in NO basis (occ-ordered)
    # semicanonicalize kept and dropped blocks separately (preserve the truncation)
    fock_ao = mf.get_fock()
    def semicanon(C):
        if C.shape[1] == 0:
            return C, np.zeros(0)
        f = C.T @ fock_ao @ C
        e, W = np.linalg.eigh(f)
        return C @ W, e
    kept, e_kept = semicanon(vir[:, :n_keep])
    drop, e_drop = semicanon(vir[:, n_keep:])
    mo_coeff = np.hstack([mf.mo_coeff[:, :nocc], kept, drop])
    mo_energy = np.concatenate([mf.mo_energy[:nocc], e_kept, e_drop])
    frozen = list(range(ncore)) + list(range(nocc + n_keep, nmo))
    first_drop = no_occ[n_keep] if n_keep < len(no_occ) else 0.0
    print(f"  FNO: nvir={nmo-nocc} -> keep {n_keep} (NO occ range {no_occ[0]:.4f}..{no_occ[n_keep-1]:.2e}; "
          f"first dropped {first_drop:.2e})", flush=True)
    return mo_coeff, mo_energy, frozen


def fno_adc2(mf, ncore, n_keep, nroots):
    mo_coeff, mo_energy, frozen = make_fno(mf, ncore, n_keep)
    mf.mo_coeff = mo_coeff; mf.mo_energy = mo_energy   # install FNO orbitals
    a = adc.ADC(mf, frozen=frozen)
    a.method = "adc(2)"; a.method_type = "ee"
    e = np.atleast_1d(a.kernel(nroots=nroots)[0])
    return e, a


if __name__ == "__main__":
    # validation: FNO-ADC(2) vs full-virtual ADC(2) on formaldehyde
    from pyscf import gto, scf
    mol = gto.M(atom='C 0 0 0; O 0 0 1.21; H 0 0.94 -0.54; H 0 -0.94 -0.54',
                basis='def2-svp', verbose=0)
    mf = scf.RHF(mol).run()
    nvir = mol.nao - mol.nelectron // 2
    a_full = adc.ADC(mf, frozen=2); a_full.method = "adc(2)"; a_full.method_type = "ee"
    e_full = np.atleast_1d(a_full.kernel(nroots=3)[0])
    for nk in (nvir, int(nvir*0.7), int(nvir*0.5)):
        e_fno, _ = fno_adc2(scf.RHF(mol).run(), 2, nk, 3)
        ev = 27.2114
        print(f"  n_keep={nk:2d}: S1={e_fno[0]*ev:.3f} eV  (full S1={e_full[0]*ev:.3f} eV, "
              f"err {abs(e_fno[0]-e_full[0])*1000*ev:.0f} meV)")

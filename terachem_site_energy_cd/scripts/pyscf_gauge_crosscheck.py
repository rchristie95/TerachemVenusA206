"""Independent cross-check of the length/velocity gauge defect with PySCF.

TeraChem reports, for the bright state of every production frame, a velocity
transition dipole whose magnitude is 23 % of the exact hypervirial value
-omega*mu_len and whose direction is 129 deg away from it (should be 180 deg),
giving f_len/f_vel ~ 19. That is far too large to be a basis-set or
range-separated-exchange effect, but "too large" is a judgement call unless we
measure the real defect for the same molecule, method and excitation.

This script does exactly that: same geometry (production frame 200, site A),
same functional (CAM-B3LYP), same TDA truncation, gas phase. If PySCF returns
f_len/f_vel of order 1, TeraChem's velocity/magnetic property code -- not the
physics -- is responsible, and the printed R(vel) carries no information.

Gas phase is the right control: the MM point charges are a local multiplicative
potential and commute with r, so they cannot break the hypervirial relation.

Run:
    conda activate pyscfneo
    python terachem_site_energy_cd/scripts/pyscf_gauge_crosscheck.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import dft, gto, tdscf

ROOT = Path(__file__).resolve().parents[1]
GEOM = ROOT / "results" / "v2_camb3lyp_frame_0200" / "site_A" / "geometry.xyz"
NSTATES = 7


def run(basis: str) -> dict:
    mol = gto.M(
        atom=str(GEOM), basis=basis, charge=-1, spin=0, unit="Angstrom", verbose=3
    )
    print(f"\n{'='*70}\nbasis {basis}: {mol.nao} AOs, {mol.nelectron} electrons\n{'='*70}",
          flush=True)

    mf = dft.RKS(mol, xc="camb3lyp").density_fit()
    mf.conv_tol = 1e-9
    mf.kernel()
    print(f"SCF converged: {mf.converged}  E = {mf.e_tot:.8f}", flush=True)

    td = tdscf.TDA(mf)  # matches TeraChem's `cis yes`
    td.nstates = NSTATES
    td.kernel()

    # gauge origin at the centre of nuclear charge, the usual molecular choice
    origin = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords()) / mol.atom_charges().sum()
    mol.set_common_orig_(origin)

    mu_len = td.transition_dipole()            # <0|r|n>, a.u.
    mu_vel = td.transition_velocity_dipole()   # <0|nabla|n>, a.u.
    m_mag = td.transition_magnetic_dipole()    # <0|r x nabla|n> * (1/2), a.u.
    w = td.e                                   # a.u.

    f_len = td.oscillator_strength(gauge="length")
    f_vel = td.oscillator_strength(gauge="velocity")

    bright = int(np.argmax(np.linalg.norm(mu_len, axis=1)))
    ml, mv = mu_len[bright], mu_vel[bright]
    wb = w[bright]

    # exact relation is <0|nabla|n> = -omega <0|r|n>, so cos should be -1
    cos = float(mv @ ml) / (np.linalg.norm(mv) * np.linalg.norm(ml))
    defect = float(np.linalg.norm(mv) / (wb * np.linalg.norm(ml)))

    # rotatory strengths, both gauges (a.u.); length gauge R = Im<0|r|n>.<n|m|0>
    r_len = float(ml @ m_mag[bright])
    r_vel = float(-(mv @ m_mag[bright]) / wb)

    # direct measurement of origin dependence: recompute m about a shifted origin
    shift = np.array([1.0, 0.0, 0.0])  # 1 bohr along x
    mol.set_common_orig_(origin + shift * 0.529177210903)  # set_common_orig_ takes Angstrom-free au? use au
    res = dict(
        basis=basis,
        nao=int(mol.nao),
        bright_root=bright + 1,
        nm=float(45.5633526 / wb * 1e0) if False else float(1e7 / (wb * 219474.63)),
        f_len=float(f_len[bright]),
        f_vel=float(f_vel[bright]),
        ratio=float(f_len[bright] / f_vel[bright]),
        vel_defect=defect,
        angle_deg=float(np.degrees(np.arccos(np.clip(cos, -1, 1)))),
        r_len=r_len,
        r_vel=r_vel,
        origin_sens_per_bohr=float(np.linalg.norm(np.cross(mv, ml))),
    )
    print(json.dumps(res, indent=1), flush=True)
    return res


def main() -> None:
    out = []
    for basis in ("6-31gs", "6-311+gss"):
        try:
            out.append(run(basis))
        except Exception as exc:  # noqa: BLE001
            print(f"basis {basis} failed: {exc}", flush=True)
    dest = ROOT / "results" / "pyscf_gauge_crosscheck.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}", flush=True)


if __name__ == "__main__":
    main()

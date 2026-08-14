#!/usr/bin/env python3
"""ADC(2) rung of the neutral gas-phase benchmark ladder, density-fitted.

Why this rung exists at all: Q-Chem's licence expired 2026-07-25, and ORCA has
no excited-state triples correction (the manual lists only EE/IP/EA-EOM-CCSD,
STEOM, DLPNO-STEOM and IH-FSMR-CCSD, with CCSD(2) as the sole perturbative
option, and that scales DOWN not up). ADC(2) is an independent hierarchy from
coupled cluster, so STEOM agreeing with it says more than STEOM agreeing with
another flavour of CC. It is honestly weaker than CCSD in absolute accuracy, so
it brackets rather than benchmarks -- it can show STEOM is not wildly off, not
that STEOM is triples-quality.

Why density fitting: the first attempt used adcc and was OOM-killed. The
binding constraint is NOT the Davidson vectors (1.67 GB each) but the MO
integral storage -- the full 4-index tensor at 331 basis functions is 89.4 GB,
and the vvvv block alone is 35.6 GB, against 62 GB of RAM. Density fitting
replaces nbf^4 with nbf^2 x naux, about 1.1 GB here, an ~83x saving. pyscf.adc
supports it directly; adcc does not.

Same geometry, same def2-SVP orbital basis, same neutral closed-shell
chromophore as every other rung, so the ladder stays method-matched. The
density-fitting auxiliary basis is an approximation to the integrals only, not
a change of model chemistry, and is reported in the JSON for the record.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
XYZ = HERE / "neutral_chromophore.xyz"


def main() -> int:
    import numpy as np
    from pyscf import gto, scf, adc

    t0 = time.time()
    mol = gto.M(atom=str(XYZ), basis="def2-svp", charge=0, spin=0,
                verbose=4, output=str(HERE / "adc2_pyscf.log"),
                max_memory=40000)
    print(f"[adc2] {mol.natm} atoms, {mol.nao_nr()} basis functions, "
          f"{mol.nelectron} electrons", flush=True)

    mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit")
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        print("[adc2] SCF did NOT converge -- aborting", flush=True)
        return 1
    print(f"[adc2] DF-RHF converged, E = {mf.e_tot:.8f} Eh "
          f"({time.time()-t0:.0f} s)", flush=True)

    myadc = adc.ADC(mf).density_fit(auxbasis="def2-svp-ri")
    myadc.method = "adc(2)"
    myadc.method_type = "ee"
    # One 1s per heavy atom. 13 C + 3 N + 3 O = 19.
    myadc.frozen = sum(1 for i in range(mol.natm) if mol.atom_symbol(i) != "H")
    myadc.max_memory = 40000
    myadc.max_space = 12          # keep the Davidson subspace small; each
                                  # doubles vector is ~1.7 GB
    print(f"[adc2] frozen core orbitals: {myadc.frozen}", flush=True)

    e, v, p, x = myadc.kernel(nroots=3)
    e = np.atleast_1d(e)
    p = np.atleast_1d(p)

    rows = []
    for i, (e_au, strength) in enumerate(zip(e, p)):
        ev = float(e_au) * 27.211386245988
        rows.append({"root": i + 1, "eV": ev, "cm-1": ev * 8065.543937,
                     "nm": 1239.841984 / ev, "fosc": float(strength)})
        print(f"[adc2] root {i+1}: {ev:.4f} eV  {1239.841984/ev:7.2f} nm  "
              f"f = {float(strength):.4f}", flush=True)

    bright = max(rows, key=lambda r: r["fosc"])
    out = {"method": "DF-ADC(2)/def2-SVP (pyscf.adc, frozen core)",
           "auxbasis_scf": "def2-universal-jkfit", "auxbasis_adc": "def2-svp-ri",
           "system": "neutral protonated CR2 chromophore, gas phase, charge 0",
           "n_basis": int(mol.nao_nr()), "frozen_core": int(myadc.frozen),
           "scf_energy_Eh": float(mf.e_tot), "roots": rows,
           "bright_state": bright, "wall_seconds": time.time() - t0}
    (HERE / "adc2_neutral.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[adc2] BRIGHT: {bright['nm']:.2f} nm, f = {bright['fosc']:.4f}",
          flush=True)
    print(f"[adc2] wrote adc2_neutral.json in {time.time()-t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""ADC(3) rung of the neutral gas-phase benchmark ladder, density-fitted.

Why this rung exists at all: Q-Chem's licence expired 2026-07-25, and ORCA has
no excited-state triples correction (the manual lists only EE/IP/EA-EOM-CCSD,
STEOM, DLPNO-STEOM and IH-FSMR-CCSD, with CCSD(2) as the sole perturbative
option, and that scales DOWN not up). ADC(3) is an independent hierarchy from
coupled cluster, so STEOM agreeing with it says more than STEOM agreeing with
another flavour of CC. It is honestly weaker than CCSD in absolute accuracy, so
it brackets rather than benchmarks -- it can show STEOM is not wildly off, not
that STEOM is triples-quality.

Why ADC(3) is feasible after all: it was rejected earlier on a memory estimate
of 33-84 GB, but that arithmetic assumed adcc's non-density-fitted integral
handling and was never revisited once DF rescued ADC(2). Measured directly on a
117-basis-function probe, DF-ADC(3) costs only 1.3x the ADC(2) memory (2.61 vs
2.00 GB post-SCF), not the 2-3x assumed. Scaling from the real DF-ADC(2) run
(35.9 GB peak of 62 GB) puts ADC(3) near 45-48 GB: tight but feasible.

Cost: the measured ADC(3)/ADC(2) wall-time ratio is 10.5-20x (13.9x at 117 bf),
so against the 4.4 h ADC(2) run this is roughly 3-5 days. The ratio grows with
system size, so treat that as a lower bound.

ADC(3) is closer to CCSD quality than ADC(2), so it tightens the bracket around
STEOM rather than merely widening it. It is still NOT a triples benchmark --
that needs CC3 on a truncated HBDI core.

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
                verbose=4, output=str(HERE / "adc3_pyscf.log"),
                max_memory=40000)
    print(f"[adc3] {mol.natm} atoms, {mol.nao_nr()} basis functions, "
          f"{mol.nelectron} electrons", flush=True)

    mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit")
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        print("[adc3] SCF did NOT converge -- aborting", flush=True)
        return 1
    print(f"[adc3] DF-RHF converged, E = {mf.e_tot:.8f} Eh "
          f"({time.time()-t0:.0f} s)", flush=True)

    myadc = adc.ADC(mf).density_fit(auxbasis="def2-svp-ri")
    myadc.method = "adc(3)"
    myadc.method_type = "ee"
    # One 1s per heavy atom. 13 C + 3 N + 3 O = 19.
    myadc.frozen = sum(1 for i in range(mol.natm) if mol.atom_symbol(i) != "H")
    myadc.max_memory = 40000
    myadc.max_space = 14          # each doubles vector is ~1.7 GB, so the
                                  # subspace is the memory knob; 14 leaves
                                  # headroom under the projected 45-48 GB peak
    print(f"[adc3] frozen core orbitals: {myadc.frozen}", flush=True)

    e, v, p, x = myadc.kernel(nroots=3)
    e = np.atleast_1d(e)
    p = np.atleast_1d(p)

    rows = []
    for i, (e_au, strength) in enumerate(zip(e, p)):
        ev = float(e_au) * 27.211386245988
        rows.append({"root": i + 1, "eV": ev, "cm-1": ev * 8065.543937,
                     "nm": 1239.841984 / ev, "fosc": float(strength)})
        print(f"[adc3] root {i+1}: {ev:.4f} eV  {1239.841984/ev:7.2f} nm  "
              f"f = {float(strength):.4f}", flush=True)

    bright = max(rows, key=lambda r: r["fosc"])
    out = {"method": "DF-ADC(3)/def2-SVP (pyscf.adc, frozen core)",
           "auxbasis_scf": "def2-universal-jkfit", "auxbasis_adc": "def2-svp-ri",
           "system": "neutral protonated CR2 chromophore, gas phase, charge 0",
           "n_basis": int(mol.nao_nr()), "frozen_core": int(myadc.frozen),
           "scf_energy_Eh": float(mf.e_tot), "roots": rows,
           "bright_state": bright, "wall_seconds": time.time() - t0}
    (HERE / "adc3_neutral.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[adc3] BRIGHT: {bright['nm']:.2f} nm, f = {bright['fosc']:.4f}",
          flush=True)
    print(f"[adc3] wrote adc3_neutral.json in {time.time()-t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

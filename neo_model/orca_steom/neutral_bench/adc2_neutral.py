#!/usr/bin/env python3
"""ADC(2) rung of the neutral gas-phase benchmark ladder.

Q-Chem's licence expired on 25 July 2026, and ORCA has no excited-state triples
correction of any kind (the manual lists EE/IP/EA-EOM-CCSD, STEOM-CCSD,
DLPNO-STEOM-CCSD and IH-FSMR-CCSD, and the only perturbative option is the
*lower*-scaling CCSD(2)). CC3 and CCSDT are formally available in free codes but
scale as N^7 and N^8; at 331 basis functions that is ~330x and ~10^5x the cost
of the EOM-CCSD rung, so neither is reachable here.

ADC(3) was the first choice but cannot run here: with 57 active occupied and
263 virtual orbitals a single doubles vector is 1.67 GB, so a Davidson subspace
of 20-50 vectors needs 33-84 GB against 62 GB of system RAM, most of which the
concurrent DLPNO-STEOM job is already holding. ADC(2) uses the same vector size
but far fewer of them and scales N^5 rather than N^6, so it fits.

ADC(2) is an independent hierarchy from coupled cluster, which is the point:
STEOM agreeing with it says more than STEOM agreeing with another flavour of
CC. It is honestly weaker than CCSD in absolute accuracy, so it brackets rather
than benchmarks -- it can show STEOM is not wildly off, not that STEOM is
triples-quality. A genuine triples benchmark needs a truncated HBDI-core model,
which is a separate rung. ADC(2) also reports oscillator strengths, which the
Q-Chem runs never produced because CC_EOM_PROP was never set.

This method has already been validated on this machine: ADC(2) gave 476 nm for
the cyanine anion against an experimental ~480, where TDDFT gave 390.

Same geometry, same def2-SVP basis, same neutral (protonated, closed-shell)
chromophore as every other rung, so the ladder is method-matched throughout.
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
    from pyscf import gto, scf
    import adcc

    t0 = time.time()
    mol = gto.M(atom=str(XYZ), basis="def2-svp", charge=0, spin=0,
                verbose=4, output=str(HERE / "adc2_pyscf.log"),
                max_memory=48000)
    print(f"[adc2] {mol.natm} atoms, {mol.nao_nr()} basis functions, "
          f"{mol.nelectron} electrons", flush=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10          # adcc wants a tight reference
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        print("[adc2] SCF did NOT converge -- aborting", flush=True)
        return 1
    print(f"[adc2] RHF converged, E = {mf.e_tot:.8f} Eh "
          f"({time.time()-t0:.0f} s)", flush=True)

    adcc.set_n_threads(8)
    # Freeze the 1s cores; they carry no valence excitation and dominate cost.
    # adcc 0.17 has no automatic frozen-core detection (it raises
    # NotImplementedError on frozen_core=True), so the count is explicit: one
    # 1s per heavy atom, 13 C + 3 N + 3 O = 19.
    n_core = sum(1 for i in range(mol.natm) if mol.atom_symbol(i) != "H")
    print(f"[adc2] freezing {n_core} core orbitals", flush=True)
    state = adcc.adc2(mf, n_singlets=3, conv_tol=1e-5, max_iter=200,
                      frozen_core=n_core)
    print(state.describe(), flush=True)

    rows = []
    for i, (e_au, f_len) in enumerate(zip(state.excitation_energy,
                                          state.oscillator_strength)):
        ev = float(e_au) * 27.211386245988
        rows.append({"root": i + 1, "eV": ev,
                     "cm-1": ev * 8065.543937,
                     "nm": 1239.841984 / ev,
                     "fosc": float(f_len)})
        print(f"[adc2] root {i+1}: {ev:.4f} eV  {1239.841984/ev:7.2f} nm  "
              f"f = {float(f_len):.4f}", flush=True)

    bright = max(rows, key=lambda r: r["fosc"])
    out = {"method": "ADC(2)/def2-SVP (adcc, frozen core)",
           "system": "neutral protonated CR2 chromophore, gas phase, charge 0",
           "n_basis": int(mol.nao_nr()), "scf_energy_Eh": float(mf.e_tot),
           "roots": rows, "bright_state": bright,
           "wall_seconds": time.time() - t0}
    (HERE / "adc2_neutral.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[adc2] BRIGHT: {bright['nm']:.2f} nm, f = {bright['fosc']:.4f}",
          flush=True)
    print(f"[adc2] wrote adc2_neutral.json in {time.time()-t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

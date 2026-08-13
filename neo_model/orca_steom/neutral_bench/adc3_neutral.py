#!/usr/bin/env python3
"""ADC(3) rung of the neutral gas-phase benchmark ladder.

Q-Chem's licence expired on 25 July 2026, and ORCA has no excited-state triples
correction of any kind (the manual lists EE/IP/EA-EOM-CCSD, STEOM-CCSD,
DLPNO-STEOM-CCSD and IH-FSMR-CCSD, and the only perturbative option is the
*lower*-scaling CCSD(2)). CC3 and CCSDT are formally available in free codes but
scale as N^7 and N^8; at 331 basis functions that is ~330x and ~10^5x the cost
of the EOM-CCSD rung, so neither is reachable here.

ADC(3) is the practical substitute and is arguably the better validation
anyway: it is an independent hierarchy from coupled cluster, so agreement
between STEOM and ADC(3) is a stronger statement than agreement between two
flavours of CC. It also reports oscillator strengths, which the Q-Chem runs
never produced because CC_EOM_PROP was never set.

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
                verbose=4, output=str(HERE / "adc3_pyscf.log"),
                max_memory=48000)
    print(f"[adc3] {mol.natm} atoms, {mol.nao_nr()} basis functions, "
          f"{mol.nelectron} electrons", flush=True)

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10          # adcc wants a tight reference
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        print("[adc3] SCF did NOT converge -- aborting", flush=True)
        return 1
    print(f"[adc3] RHF converged, E = {mf.e_tot:.8f} Eh "
          f"({time.time()-t0:.0f} s)", flush=True)

    adcc.set_n_threads(8)
    # Freeze the 1s cores; they carry no valence excitation and dominate cost.
    state = adcc.adc3(mf, n_singlets=3, conv_tol=1e-5, max_iter=200,
                      frozen_core=True)
    print(state.describe(), flush=True)

    rows = []
    for i, (e_au, f_len) in enumerate(zip(state.excitation_energy,
                                          state.oscillator_strength)):
        ev = float(e_au) * 27.211386245988
        rows.append({"root": i + 1, "eV": ev,
                     "cm-1": ev * 8065.543937,
                     "nm": 1239.841984 / ev,
                     "fosc": float(f_len)})
        print(f"[adc3] root {i+1}: {ev:.4f} eV  {1239.841984/ev:7.2f} nm  "
              f"f = {float(f_len):.4f}", flush=True)

    bright = max(rows, key=lambda r: r["fosc"])
    out = {"method": "ADC(3)/def2-SVP (adcc, frozen core)",
           "system": "neutral protonated CR2 chromophore, gas phase, charge 0",
           "n_basis": int(mol.nao_nr()), "scf_energy_Eh": float(mf.e_tot),
           "roots": rows, "bright_state": bright,
           "wall_seconds": time.time() - t0}
    (HERE / "adc3_neutral.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"[adc3] BRIGHT: {bright['nm']:.2f} nm, f = {bright['fosc']:.4f}",
          flush=True)
    print(f"[adc3] wrote adc3_neutral.json in {time.time()-t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

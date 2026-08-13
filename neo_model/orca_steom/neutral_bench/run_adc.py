#!/usr/bin/env python3
"""ADC(2) and ADC(3) on the neutral gas-phase chromophore, for the benchmark ladder.

The EOM-CCSD(fT) row of the ladder cannot be reproduced: ORCA has no
excited-state triples correction of any kind (the ORCA 6.1 manual lists only
EE/IP/EA-EOM-CCSD, STEOM-CCSD, DLPNO-STEOM-CCSD and IH-FSMR-CCSD, with CCSD(2)
as a *lower*-scaling truncation), and the Q-Chem licence expired on 2026-07-25.

ADC(3) is the practical substitute and is arguably a better validation anyway,
because it is an independent method family rather than another flavour of
coupled cluster: agreement between STEOM and ADC(3) tests the same claim
without sharing the CC ansatz. Proper CC triples are out of reach at this size
regardless of code -- 331 basis functions puts CC3 at ~330x the cost of the
EOM-CCSD reference and CCSDT at ~10^5x.

Same neutral (protonated) geometry and def2-SVP basis as every other row, so
the ladder stays method-matched. Reports excitation energies AND oscillator
strengths, the latter being what the Q-Chem rows never produced.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
from pyscf import gto, scf
import adcc

HERE = Path(__file__).resolve().parent
XYZ = HERE / "neutral_chromophore.xyz"
NSTATES = 3

mol = gto.M(atom=str(XYZ), basis="def2-svp", charge=0, spin=0,
            verbose=0, max_memory=48000)
print(f"[scf] {mol.natm} atoms, {mol.nao_nr()} basis functions, "
      f"{mol.nelectron} electrons", flush=True)

t0 = time.time()
mf = scf.RHF(mol)
mf.conv_tol = 1e-10
mf.kernel()
print(f"[scf] E = {mf.e_tot:.8f} Ha, converged = {mf.converged}, "
      f"{time.time()-t0:.0f} s", flush=True)
assert mf.converged, "SCF did not converge"

results = {"system": {"atoms": int(mol.natm), "nbas": int(mol.nao_nr()),
                      "nelec": int(mol.nelectron), "basis": "def2-SVP",
                      "charge": 0, "scf_energy_Ha": float(mf.e_tot)}}

for level, fn in (("adc2", adcc.adc2), ("adc3", adcc.adc3)):
    t0 = time.time()
    print(f"\n[{level}] starting ({NSTATES} singlet states) ...", flush=True)
    try:
        state = fn(mf, n_singlets=NSTATES, conv_tol=1e-5)
        rows = []
        for i in range(len(state.excitation_energy)):
            ev = float(state.excitation_energy[i] * 27.211386245988)
            f = float(state.oscillator_strength[i])
            rows.append({"root": i + 1, "eV": ev, "nm": 1239.841984 / ev,
                         "fosc": f})
            print(f"   root {i+1}:  {ev:8.4f} eV   {1239.841984/ev:7.1f} nm   "
                  f"f = {f:.4f}", flush=True)
        bright = max(rows, key=lambda r: r["fosc"])
        print(f"   BRIGHT: {bright['nm']:.1f} nm, f = {bright['fosc']:.4f}  "
              f"({time.time()-t0:.0f} s)", flush=True)
        results[level] = {"states": rows, "bright": bright,
                          "wall_seconds": time.time() - t0}
    except Exception as exc:
        print(f"   {level} FAILED: {type(exc).__name__}: {exc}", flush=True)
        results[level] = {"error": f"{type(exc).__name__}: {exc}"}

(HERE / "adc_results.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"\nwrote {HERE/'adc_results.json'}")

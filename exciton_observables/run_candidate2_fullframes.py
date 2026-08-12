#!/usr/bin/env python3
"""Generate full-system frames of the candidate-2 dimer for QM/MM detuning.

The 20 ns production DCD is protein-only, but the site-energy embedding needs
solvent coordinates (full-system embedding, matching the v2 baseline policy).
This continues the candidate-2 system from its final relaxed state for a short
window, writing FULL-system frames:

  200 ps re-equilibration (discarded) + 3.2 ns production,
  full-system DCD every 40 ps -> 80 statistically independent frames
  (the site-energy detuning decorrelates in <= 5 ps).

Start state: classical_relaxed.pdb positions with the box taken from the
production DCD cell (88.2695 A) -- NEVER the PDB's CRYST1, which is the stale
pre-NPT box (see the v3 stale-box trap; same failure mode applies here).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path("/home/robson/PetaChem")
sys.path.insert(0, str(REPO))

from run_nvt import (  # noqa: E402  (configures OPENMM_PLUGIN_DIR for CUDA)
    apply_amber_cr2_parameters,
    ensure_amber_cr2_topology,
    pick_platform,
)

import openmm  # noqa: E402
from openmm import unit  # noqa: E402
from openmm.app import DCDReporter, ForceField, HBonds, PDBFile, PME, Simulation  # noqa: E402

WORKDIR = REPO / "tc_candidate2_ff19sb_opc"
CR2_PRMTOP = REPO / "anionic_build" / "monomer_solv.prmtop"
BOX_A = json.loads((REPO / "exciton_observables/candidate2_box.json").read_text())["box_a"]

EQUIL_STEPS = 100_000        # 200 ps
PROD_STEPS = 1_600_000       # 3.2 ns
DCD_INTERVAL = 20_000        # 40 ps -> 80 frames
SEED = 20260813


def main():
    pdb = PDBFile(str(WORKDIR / "classical_relaxed.pdb"))
    topology = pdb.topology
    topology.setPeriodicBoxVectors(np.diag([BOX_A / 10.0] * 3) * unit.nanometer)

    ensure_amber_cr2_topology(topology, CR2_PRMTOP)
    forcefield = ForceField("amber19-all.xml", "amber14/opc.xml",
                            str(WORKDIR / "nonstandard_residues_generic.xml"))
    system = forcefield.createSystem(topology, nonbondedMethod=PME,
                                     nonbondedCutoff=1.0 * unit.nanometer,
                                     constraints=HBonds)
    apply_amber_cr2_parameters(system, topology, CR2_PRMTOP)
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    net = sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
              for i in range(system.getNumParticles()))
    print(f"[build] {system.getNumParticles()} particles, net charge {net:+.6f} e, "
          f"box {BOX_A:.4f} A", flush=True)
    if abs(net) > 1.0e-4:
        raise RuntimeError("net charge out of tolerance")

    integrator = openmm.LangevinMiddleIntegrator(300.0 * unit.kelvin,
                                                 1.0 / unit.picosecond,
                                                 2.0 * unit.femtoseconds)
    integrator.setRandomNumberSeed(SEED)
    platform = pick_platform("CUDA", strict=True)
    simulation = Simulation(topology, system, integrator, platform,
                            {"DeviceIndex": "0", "Precision": "mixed"})
    simulation.context.setPositions(pdb.positions)
    simulation.context.computeVirtualSites()
    simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, SEED)

    t0 = time.time()
    print(f"[equil] {EQUIL_STEPS*2e-6:.1f} ns discard ...", flush=True)
    simulation.step(EQUIL_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"[equil] done, PE = {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.0f} kJ/mol, "
          f"{(time.time()-t0)/60:.1f} min", flush=True)

    out = WORKDIR / "candidate2_fullsystem.dcd"
    simulation.reporters.append(DCDReporter(str(out), DCD_INTERVAL))
    print(f"[prod] {PROD_STEPS*2e-6:.1f} ns, full-system frame every "
          f"{DCD_INTERVAL*2e-3:.0f} ps -> {PROD_STEPS//DCD_INTERVAL} frames", flush=True)
    for chunk in range(PROD_STEPS // 100_000):
        simulation.step(100_000)
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(pe):
            raise RuntimeError("PE diverged")
        print(f"[prod] {(chunk+1)*0.2:.1f} ns  PE={pe:.0f}", flush=True)
    print(f"[done] {out}  ({(time.time()-t0)/60:.1f} min total)", flush=True)


if __name__ == "__main__":
    main()

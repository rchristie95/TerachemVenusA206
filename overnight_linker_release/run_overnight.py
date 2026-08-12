#!/usr/bin/env python3
"""Overnight test of the linker-tension hypothesis for the tandem register.

Background (notes/dimer_geometry_audit.md): the tandem was built in the 1MYW
crystal register, which forces the 33-residue linker (res 230-262) to span
~54 A against ~38 A natural end-to-end; across all v3 production frames it sits
pinned at 52.3 +/- 0.7 A. Both force fields move the interface angle AWAY from
the 131.3 deg the limiting anisotropy requires, and 1 ns cannot re-dock two
beta-barrels, so an alternative register has never been sampled.

Two arms, same v3 ff19SB/OPC system (rebuilt from the retained solvated box,
CR2 physics transplanted from the retained AMBER prmtop exactly as in
run_nvt.py):

  control : unbiased 300 K NVT continuation from the v3 end state. Null arm --
            does theta stay pinned given >10x the previous sampling time?
  release : steered harmonic on the linker end-to-end distance
            (CA229-CA263): ramp 52.3 -> 40 A over RAMP_NS, hold HOLD_NS,
            then set k=0 and run unbiased for the rest of the night.
            If linker tension holds the crystal register, the inter-axis
            angle should move once slack is granted -- and the direction it
            moves (toward or away from 131 deg) is the result.

Every 10 ps both arms log: inter-CR2-long-axis angle alpha (audit route 3:
OH minus imidazolinone-ring centroid, no dipoles needed), cos(alpha), linker
end-to-end, CR2 centroid separation, chirality triple product
R_AB.(axA x axB), potential energy and temperature.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path("/home/robson/PetaChem")
sys.path.insert(0, str(REPO))

# Importing run_nvt configures OPENMM_PLUGIN_DIR so CUDA is visible.
from run_nvt import (  # noqa: E402
    apply_amber_cr2_parameters,
    ensure_amber_cr2_topology,
    pick_platform,
)

import openmm  # noqa: E402
from openmm import unit  # noqa: E402
from openmm.app import (  # noqa: E402
    CheckpointReporter,
    DCDReporter,
    ForceField,
    HBonds,
    PDBFile,
    PME,
    Simulation,
)

V3_DIR = REPO / "tc_tandem_nvt_v3_ff19sb_opc"
START_PDB = V3_DIR / "classical_relaxed.pdb"  # topology source only
NONSTANDARD_XML = V3_DIR / "nonstandard_residues_generic.xml"
CR2_PRMTOP = REPO / "anionic_build" / "monomer_solv.prmtop"

# True v3 end state: classical_relaxed.pdb carries the pre-NPT CRYST1 box
# (90.373 A) with unwrapped coordinates, which overlaps catastrophically when
# re-imaged. The production box after NPT is 88.198 A; positions and cell are
# taken from the last frame of the production DCD instead.
HERE = Path(__file__).resolve().parent
START_XYZ_NPY = HERE / "v3_last_frame_xyz_ang.npy"
START_BOX_JSON = HERE / "v3_last_frame_box.json"

RING = ("CA2", "C2", "N2", "C1", "N3")  # imidazolinone ring, as in the audit
LINKER_ANCHORS = ("229", "263")  # CA of barrel-1 C-term, CA of barrel-2 N-term

CHUNK_STEPS = 5000          # 10 ps at 2 fs: metrics + bias-schedule update cadence
DCD_INTERVAL = 12500        # 25 ps
CHECKPOINT_INTERVAL = 50000  # 100 ps


def build_simulation(workdir: Path, seed: int):
    pdb = PDBFile(str(START_PDB))
    topology = pdb.topology

    box_info = json.loads(START_BOX_JSON.read_text())
    xyz_ang = np.load(START_XYZ_NPY)
    if xyz_ang.shape[0] != topology.getNumAtoms():
        raise RuntimeError(
            f"DCD frame has {xyz_ang.shape[0]} atoms, topology {topology.getNumAtoms()}")
    box_nm = box_info["box_a"] / 10.0
    topology.setPeriodicBoxVectors(np.diag([box_nm] * 3) * unit.nanometer)
    positions = (xyz_ang / 10.0) * unit.nanometer
    print(f"[build] start = production DCD frame {box_info['frame_index']} "
          f"(box {box_info['box_a']:.3f} A)", flush=True)

    cr2_summary = ensure_amber_cr2_topology(topology, CR2_PRMTOP)
    print(f"[build] CR2 topology: bonds/site={cr2_summary['bonds_per_site']}, "
          f"added={cr2_summary['bonds_added']}", flush=True)

    forcefield = ForceField("amber19-all.xml", "amber14/opc.xml", str(NONSTANDARD_XML))
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    transplant = apply_amber_cr2_parameters(system, topology, CR2_PRMTOP)
    print(f"[build] CR2 transplant charges: {transplant['charge_per_cr2_e']}", flush=True)

    nonbonded = next(f for f in system.getForces()
                     if isinstance(f, openmm.NonbondedForce))
    net = sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
              for i in range(system.getNumParticles()))
    print(f"[build] net charge after transplant: {net:+.6f} e", flush=True)
    if abs(net) > 1.0e-4:
        raise RuntimeError(f"net charge {net:+.6f} e out of tolerance")

    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)
    integrator.setRandomNumberSeed(seed)

    platform = pick_platform("CUDA", strict=True)
    properties = {"DeviceIndex": "0", "Precision": "mixed"}
    simulation = Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    simulation.context.computeVirtualSites()
    return simulation, topology, system


def analysis_indices(topology):
    """Atom indices for the on-the-fly geometry metrics."""
    cr2 = []
    anchors = {}
    for residue in topology.residues():
        if residue.name == "CR2":
            site = {}
            for atom in residue.atoms():
                if atom.name in RING or atom.name == "OH":
                    site[atom.name] = atom.index
            missing = [n for n in (*RING, "OH") if n not in site]
            if missing:
                raise RuntimeError(f"CR2 {residue.id} missing atoms {missing}")
            site["all_heavy"] = [a.index for a in residue.atoms()
                                 if a.element is not None and a.element.symbol != "H"]
            cr2.append(site)
        if residue.id in LINKER_ANCHORS:
            for atom in residue.atoms():
                if atom.name == "CA":
                    anchors[residue.id] = atom.index
    if len(cr2) != 2:
        raise RuntimeError(f"expected 2 CR2 residues, found {len(cr2)}")
    if set(anchors) != set(LINKER_ANCHORS):
        raise RuntimeError(f"linker anchor CAs not found: {anchors}")
    return cr2, (anchors[LINKER_ANCHORS[0]], anchors[LINKER_ANCHORS[1]])


def metrics(pos_ang: np.ndarray, cr2, anchor_pair):
    out = {}
    axes, cents = [], []
    for site in cr2:
        ring_c = pos_ang[[site[n] for n in RING]].mean(axis=0)
        axis = pos_ang[site["OH"]] - ring_c
        axes.append(axis / np.linalg.norm(axis))
        cents.append(pos_ang[site["all_heavy"]].mean(axis=0))
    cos_a = float(np.dot(axes[0], axes[1]))
    out["alpha_deg"] = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
    out["cos_alpha"] = cos_a
    r_ab = cents[1] - cents[0]
    out["cr2_sep_ang"] = float(np.linalg.norm(r_ab))
    out["triple_product_ang"] = float(np.dot(r_ab, np.cross(axes[0], axes[1])))
    out["linker_e2e_ang"] = float(
        np.linalg.norm(pos_ang[anchor_pair[0]] - pos_ang[anchor_pair[1]]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", choices=("control", "release"), required=True)
    ap.add_argument("--hours", type=float, default=8.5, help="wall-clock budget")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--resume", type=Path, default=None,
                    help="Checkpoint to continue from; velocities and box are taken "
                         "from it rather than re-initialised.")
    ap.add_argument("--tag", default=None,
                    help="Output subdirectory name (default: the arm name). Use for "
                         "continuations so the original run's data is never overwritten.")
    ap.add_argument("--t0-ns", type=float, default=0.0,
                    help="Simulated time already accumulated, for continuous logging.")
    ap.add_argument("--ramp-ns", type=float, default=2.0)
    ap.add_argument("--hold-ns", type=float, default=2.0)
    ap.add_argument("--target-ang", type=float, default=40.0)
    ap.add_argument("--k-pull", type=float, default=1500.0,
                    help="kJ/mol/nm^2 on the linker end-to-end distance")
    args = ap.parse_args()

    seed = args.seed if args.seed is not None else (20260811 if args.arm == "control" else 20260812)
    workdir = Path(__file__).resolve().parent / (args.tag or args.arm)
    workdir.mkdir(parents=True, exist_ok=True)
    stem = args.tag or args.arm

    simulation, topology, system = build_simulation(workdir, seed)
    cr2, anchor_pair = analysis_indices(topology)

    # Steering force: harmonic on the linker end-to-end distance, controlled
    # through global parameters so the schedule never rebuilds the context.
    if args.arm == "release":
        bias = openmm.CustomBondForce("0.5*k_pull*(r-r0_pull)^2")
        bias.addGlobalParameter("k_pull", 0.0)
        bias.addGlobalParameter("r0_pull", 5.23)
        bias.addBond(anchor_pair[0], anchor_pair[1], [])
        bias.setForceGroup(15)
        system.addForce(bias)
        simulation.context.reinitialize(preserveState=True)

    if args.resume is not None:
        simulation.loadCheckpoint(str(args.resume))
        print(f"[resume] loaded {args.resume}; continuing with its positions, "
              f"velocities and box (t0 = {args.t0_ns} ns)", flush=True)
    else:
        simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, seed)

    protein_atoms = [a.index for a in topology.atoms()
                     if a.residue.name not in ("HOH", "WAT", "SOL", "NA", "CL", "K")]
    dcd_path = workdir / f"{stem}.dcd"
    simulation.reporters.append(
        DCDReporter(str(dcd_path), DCD_INTERVAL, atomSubset=protein_atoms))
    simulation.reporters.append(
        CheckpointReporter(str(workdir / f"{stem}.chk"), CHECKPOINT_INTERVAL))

    state0 = simulation.context.getState(getPositions=True)
    pos0 = state0.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    m0 = metrics(np.asarray(pos0), cr2, anchor_pair)
    print(f"[start] {json.dumps(m0)}", flush=True)

    ramp_steps = int(args.ramp_ns * 500000)
    hold_steps = int(args.hold_ns * 500000)
    r_start_nm = m0["linker_e2e_ang"] / 10.0
    r_target_nm = args.target_ang / 10.0

    csv_path = workdir / f"{stem}_metrics.csv"
    fields = ["time_ps", "phase", "k_pull", "r0_pull_ang", "alpha_deg", "cos_alpha",
              "linker_e2e_ang", "cr2_sep_ang", "triple_product_ang",
              "potential_kj_mol", "temperature_k"]
    csv_handle = open(csv_path, "w", newline="")
    writer = csv.DictWriter(csv_handle, fieldnames=fields)
    writer.writeheader()

    ndof = system.getNumParticles() * 3 - system.getNumConstraints() - 3
    # OPC virtual sites carry no DOF
    n_virtual = sum(1 for i in range(system.getNumParticles())
                    if system.isVirtualSite(i))
    ndof -= 3 * n_virtual

    deadline = time.time() + args.hours * 3600.0
    step = 0
    wall0 = time.time()
    print(f"[run] arm={args.arm} seed={seed} budget={args.hours} h "
          f"start linker={m0['linker_e2e_ang']:.2f} A alpha={m0['alpha_deg']:.2f} deg",
          flush=True)

    while time.time() < deadline:
        if args.arm == "release":
            if step < ramp_steps:
                phase = "ramp"
                frac = step / max(ramp_steps, 1)
                k_now = args.k_pull * min(1.0, frac * 10.0)  # full k within first 10% of ramp
                r0_now = r_start_nm + (r_target_nm - r_start_nm) * frac
            elif step < ramp_steps + hold_steps:
                phase, k_now, r0_now = "hold", args.k_pull, r_target_nm
            else:
                phase, k_now, r0_now = "released", 0.0, r_target_nm
            simulation.context.setParameter("k_pull", k_now)
            simulation.context.setParameter("r0_pull", r0_now)
        else:
            phase, k_now, r0_now = "control", 0.0, 0.0

        simulation.step(CHUNK_STEPS)
        step += CHUNK_STEPS

        state = simulation.context.getState(getPositions=True, getEnergy=True)
        pos = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
        row = metrics(pos, cr2, anchor_pair)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        temp = 2.0 * ke * 1000.0 / (ndof * unit.MOLAR_GAS_CONSTANT_R.value_in_unit(
            unit.joule / (unit.kelvin * unit.mole)))
        if not np.isfinite(pe):
            raise RuntimeError(f"potential energy diverged at step {step}: {pe}")
        row.update(time_ps=args.t0_ns * 1000.0 + step * 0.002, phase=phase, k_pull=k_now,
                   r0_pull_ang=r0_now * 10.0, potential_kj_mol=pe,
                   temperature_k=temp)
        writer.writerow(row)
        csv_handle.flush()

        if step % 50000 == 0:  # every 100 ps
            ns_per_day = (step * 2.0e-6) / max(time.time() - wall0, 1.0) * 86400.0
            print(f"[{stem}] t={args.t0_ns + step*0.002/1000.0:8.3f} ns  phase={phase:8s}  "
                  f"linker={row['linker_e2e_ang']:6.2f} A  alpha={row['alpha_deg']:7.2f}  "
                  f"cos={row['cos_alpha']:+.3f}  sep={row['cr2_sep_ang']:6.2f} A  "
                  f"T={temp:6.1f} K  {ns_per_day:5.1f} ns/day", flush=True)

    total_ns = step * 2.0e-6
    print(f"[done] arm={args.arm} sampled {total_ns:.2f} ns in "
          f"{(time.time()-wall0)/3600.0:.2f} h", flush=True)
    csv_handle.close()


if __name__ == "__main__":
    main()

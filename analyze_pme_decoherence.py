#!/usr/bin/env python3
"""Compare minimum-image and full PME electrostatic gap fluctuations.

The STEOM difference charges are placed on the 29 actual CR2 atoms.  For the
PME calculation, each site's own CR2 and stacked Tyr MM charges are zeroed and
the difference charges are introduced through an OpenMM global-parameter
offset.  The centred finite-difference derivative with respect to that parameter is the
linear electrostatic contribution to the excitation-energy gap, including the
real- and reciprocal-space PME terms.  A direct minimum-image calculation is
performed on the same positions and charges for a controlled comparison.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
from openmm import Context, NonbondedForce, Platform, Vec3, VerletIntegrator, XmlSerializer, unit
from openmm.app import ForceField, HBonds, PME, PDBFile

import run_nvt
from analyze_solvation_decoherence import (
    ElectrostaticEngine,
    HARTREE_TO_CM,
    analyze_series,
    classify_atoms,
)


KJ_MOL_TO_CM = 83.593472251353


def build_base_system(topology_pdb: Path, generic_xml: Path, amber_prmtop: Path):
    pdb = PDBFile(str(topology_pdb))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml", str(generic_xml))
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=run_nvt.choose_cutoff_from_box(pdb.topology, default_nm=1.0),
        constraints=HBonds,
        ignoreExternalBonds=True,
    )
    run_nvt.apply_amber_cr2_parameters(system, pdb.topology, amber_prmtop)
    atoms = list(pdb.topology.atoms())
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    charges = np.empty(system.getNumParticles(), dtype=np.float32)
    for index in range(system.getNumParticles()):
        charge, _, _ = nonbonded.getParticleParameters(index)
        charges[index] = charge.value_in_unit(unit.elementary_charge)
    return pdb, system, atoms, charges


def cr2_name_maps(atoms: list) -> list[dict[str, int]]:
    residues = []
    seen = set()
    for atom in atoms:
        if atom.residue.name != "CR2":
            continue
        key = (atom.residue.chain.id, atom.residue.id)
        if key not in seen:
            residues.append(atom.residue)
            seen.add(key)
    if len(residues) != 2:
        raise RuntimeError(f"Expected two CR2 residues, found {len(residues)}")
    return [{atom.name: atom.index for atom in residue.atoms()} for residue in residues]


def make_site_system(base_xml: str, exclude_code: np.ndarray, site: int, atom_indices: list[int], dq: np.ndarray):
    system = XmlSerializer.deserialize(base_xml)
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    for index in np.flatnonzero(exclude_code == site):
        _, sigma, epsilon = nonbonded.getParticleParameters(int(index))
        nonbonded.setParticleParameters(int(index), 0.0, sigma, epsilon)
    nonbonded.addGlobalParameter("lambda_gap", 0.0)
    for index, delta_q in zip(atom_indices, dq):
        nonbonded.addParticleParameterOffset("lambda_gap", int(index), float(delta_q), 0.0, 0.0)
    for force in system.getForces():
        force.setForceGroup(0 if isinstance(force, NonbondedForce) else 1)
    return system


def set_periodic_box(context: Context, vectors_nm: np.ndarray) -> None:
    vectors = [Vec3(*row) * unit.nanometer for row in np.asarray(vectors_nm, dtype=float)]
    context.setPeriodicBoxVectors(*vectors)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectories", type=Path, nargs="+", required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--generic-xml", type=Path, required=True)
    parser.add_argument("--amber-cr2-prmtop", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--dt-fs", type=float, required=True)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=-1)
    parser.add_argument("--finite-difference-scale", type=float, default=10.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    probe = np.load(args.probe)
    cr2_count = int(probe["cr2_atom_count"])
    cr2_names = [str(value) for value in probe["cr2_atom_names"][:cr2_count]]
    dq = np.asarray(probe["atom_delta_q_e"][:cr2_count], dtype=float).copy()
    # The omitted Tyr tail carries <0.5% of the absolute difference charge.
    # Remove its tiny compensating monopole uniformly so the PME probe is exact neutral.
    neutralization_per_atom = float(dq.sum() / len(dq))
    dq -= neutralization_per_atom

    pdb, base_system, atoms, mm_q = build_base_system(
        args.topology, args.generic_xml, args.amber_cr2_prmtop
    )
    group_code, exclude_code, topology_audit = classify_atoms(atoms)
    name_maps = cr2_name_maps(atoms)
    site_indices = []
    for mapping in name_maps:
        missing = set(cr2_names) - set(mapping)
        if missing:
            raise RuntimeError(f"Missing CR2 atoms: {sorted(missing)}")
        site_indices.append([mapping[name] for name in cr2_names])

    base_xml = XmlSerializer.serialize(base_system)
    pme_energies = []
    direct_energies = []
    total_frames = None
    platform = Platform.getPlatformByName("OpenCL")

    for site in range(2):
        system = make_site_system(base_xml, exclude_code, site, site_indices[site], dq)
        integrator = VerletIntegrator(1.0 * unit.femtosecond)
        context = Context(
            system,
            integrator,
            platform,
            {"DeviceIndex": "0", "Precision": "mixed"},
        )
        context.setParameter("lambda_gap", 0.0)
        direct_engine = ElectrostaticEngine(dq, mm_q, group_code, exclude_code)
        site_pme = []
        site_direct = []
        frame_count = 0
        for trajectory in args.trajectories:
            for chunk in md.iterload(
                str(trajectory),
                top=str(args.topology),
                chunk=10,
                stride=max(1, args.stride),
            ):
                if chunk.unitcell_vectors is None:
                    raise RuntimeError("Trajectory lacks periodic unit-cell vectors")
                for local in range(chunk.n_frames):
                    if args.max_frames >= 0 and frame_count >= args.max_frames:
                        break
                    xyz_nm = np.asarray(chunk.xyz[local], dtype=np.float32)
                    box_nm = np.asarray(chunk.unitcell_vectors[local], dtype=float)
                    context.setPositions(xyz_nm * unit.nanometer)
                    set_periodic_box(context, box_nm)
                    scale = float(args.finite_difference_scale)
                    context.setParameter("lambda_gap", scale)
                    energy_plus = context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
                    context.setParameter("lambda_gap", -scale)
                    energy_minus = context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
                    context.setParameter("lambda_gap", 0.0)
                    derivative = (
                        energy_plus.value_in_unit(unit.kilojoule_per_mole)
                        - energy_minus.value_in_unit(unit.kilojoule_per_mole)
                    ) / (2.0 * scale)
                    site_pme.append(float(derivative) * KJ_MOL_TO_CM)

                    xyz_ang = xyz_nm * 10.0
                    direct_engine.set_mm_positions(xyz_ang)
                    direct_groups = direct_engine.calculate(
                        xyz_ang[np.asarray(site_indices[site], dtype=int)],
                        np.linalg.norm(box_nm, axis=1) * 10.0,
                        site,
                    )
                    site_direct.append(float(np.sum(direct_groups)))
                    frame_count += 1
                    if frame_count == 1 or frame_count % 200 == 0:
                        print(f"[site {site + 1}] {frame_count} frames", flush=True)
                if args.max_frames >= 0 and frame_count >= args.max_frames:
                    break
            if args.max_frames >= 0 and frame_count >= args.max_frames:
                break
        if total_frames is None:
            total_frames = frame_count
        elif frame_count != total_frames:
            raise RuntimeError("Site frame counts differ")
        pme_energies.append(site_pme)
        direct_energies.append(site_direct)
        del context, integrator, system, direct_engine

    pme = np.asarray(pme_energies, dtype=float).T
    direct = np.asarray(direct_energies, dtype=float).T
    time_fs = np.arange(total_frames, dtype=float) * args.dt_fs * max(1, args.stride)
    pme_summary, pme_arrays = analyze_series(time_fs, pme[:, 0], pme[:, 1])
    direct_summary, direct_arrays = analyze_series(time_fs, direct[:, 0], direct[:, 1])
    correction = pme - direct
    correction_summary, _ = analyze_series(time_fs, correction[:, 0], correction[:, 1])
    summary = {
        "method": "centred OpenMM PME charge-offset derivative on actual CR2 atoms",
        "trajectories": [f"external/{path.name}" for path in args.trajectories],
        "frames": int(total_frames),
        "dt_fs": float(args.dt_fs * max(1, args.stride)),
        "duration_fs": float(time_fs[-1]) if len(time_fs) else 0.0,
        "finite_difference_scale": float(args.finite_difference_scale),
        "probe": {
            "cr2_atoms": cr2_count,
            "net_before_neutralization_e": float(np.asarray(probe["atom_delta_q_e"][:cr2_count]).sum()),
            "neutralization_per_atom_e": neutralization_per_atom,
            "net_after_neutralization_e": float(dq.sum()),
            "omitted_tyr_absolute_charge_e": float(np.abs(probe["atom_delta_q_e"][cr2_count:]).sum()),
        },
        "topology": topology_audit,
        "PME": pme_summary,
        "minimum_image_same_probe": direct_summary,
        "PME_minus_minimum_image": correction_summary,
        "site_trace_correlations": {
            "A_PME_vs_minimum_image": float(np.corrcoef(pme[:, 0], direct[:, 0])[0, 1]),
            "B_PME_vs_minimum_image": float(np.corrcoef(pme[:, 1], direct[:, 1])[0, 1]),
            "difference_PME_vs_minimum_image": float(
                np.corrcoef(pme[:, 0] - pme[:, 1], direct[:, 0] - direct[:, 1])[0, 1]
            ),
        },
        "limitations": [
            "difference density is fixed and restricted to the 29 CR2 atoms",
            "PME uses the same fixed-charge AMBER/TIP3P Hamiltonian as the trajectory",
            "particle-parameter offsets do not restore scaled 1-4 exception interactions across the QM boundary",
            "classical second cumulant; no quantum correction or geometry-dependent excitation energies",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    np.savez_compressed(
        args.out / "pme_gap_timeseries.npz",
        time_fs=time_fs,
        pme_cm=pme,
        minimum_image_cm=direct,
        pme_minus_minimum_image_cm=correction,
        pme_coherence=pme_arrays["coherence_classical"],
        direct_coherence=direct_arrays["coherence_classical"],
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

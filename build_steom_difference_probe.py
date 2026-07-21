#!/usr/bin/env python3
"""Build a quantitative STEOM S1 excited-minus-ground electrostatic probe.

The manuscript difference-density cube is normalized for visualization and is
therefore unsuitable for energies.  This script reconstructs the electron
number-density change from all printed S1 NTO pairs,

    Delta rho = sum_k n_k (|particle_k|^2 - |hole_k|^2),

normalizes every orbital density on the numerical cube, converts the result to
an electronic charge-density change, and partitions that charge onto the 41
physical QM atoms.  Density assigned nearest to one of the three link atoms is
reassigned to its nearest physical atom, preserving the exactly neutral total
charge change.

The resulting atom-centred probe is intended for fast electrostatic energy-gap
tests against an explicit protein/water MM charge trajectory.  It is not a
replacement for state-specific QM/MM excitation energies.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903

# S1 NTO occupations printed by the retained STEOM calculation.
NTO_PAIRS = [
    (92, 93, 0.96558779),
    (91, 94, 0.01476372),
    (90, 95, 0.00315100),
    (89, 96, 0.00237768),
    (88, 97, 0.00196516),
    (87, 98, 0.00149604),
    (86, 99, 0.00111854),
]


def read_cube(path: Path) -> dict[str, np.ndarray | int]:
    with path.open("r", encoding="utf-8") as handle:
        handle.readline()
        handle.readline()
        fields = handle.readline().split()
        signed_natoms = int(fields[0])
        natoms = abs(signed_natoms)
        origin = np.asarray([float(value) for value in fields[1:4]], dtype=float)
        shape = []
        axes = []
        for _ in range(3):
            fields = handle.readline().split()
            shape.append(int(fields[0]))
            axes.append([float(value) for value in fields[1:4]])
        atomic_numbers = []
        atoms_bohr = []
        for _ in range(natoms):
            fields = handle.readline().split()
            atomic_numbers.append(int(float(fields[0])))
            atoms_bohr.append([float(value) for value in fields[2:5]])
        if signed_natoms < 0:
            handle.readline()
        values = np.fromstring(" ".join(handle.read().split()), sep=" ")
    shape_array = np.asarray(shape, dtype=int)
    expected = int(np.prod(shape_array))
    if values.size != expected:
        raise ValueError(f"Cube size mismatch in {path}: {values.size} != {expected}")
    return {
        "natoms": natoms,
        "origin": origin,
        "shape": shape_array,
        "axes": np.asarray(axes, dtype=float),
        "atomic_numbers": np.asarray(atomic_numbers, dtype=int),
        "atoms_bohr": np.asarray(atoms_bohr, dtype=float),
        "values": values.reshape(tuple(shape_array)),
    }


def read_cr2_names(monomer: Path) -> list[str]:
    names = []
    with monomer.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            if line[17:20].strip() == "CR2":
                names.append(line[12:16].strip())
    if len(names) != 29:
        raise ValueError(f"Expected 29 CR2 atoms in {monomer}, found {len(names)}")
    return names


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cube-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="steom_phenol_svpd.s1")
    parser.add_argument("--monomer", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=len(NTO_PAIRS),
        help="Use only the first N occupation-ranked NTO pairs (default: all)",
    )
    parser.add_argument(
        "--match-grid-dipole",
        action="store_true",
        help="Apply the minimum-norm atomic-charge correction that exactly preserves the grid monopole and dipole",
    )
    args = parser.parse_args(argv)

    if args.max_pairs < 1 or args.max_pairs > len(NTO_PAIRS):
        parser.error(f"--max-pairs must be between 1 and {len(NTO_PAIRS)}")
    selected_pairs = NTO_PAIRS[: args.max_pairs]

    reference = read_cube(args.cube_dir / f"{args.prefix}.mo92a.cube")
    shape = np.asarray(reference["shape"], dtype=int)
    origin = np.asarray(reference["origin"], dtype=float)
    axes = np.asarray(reference["axes"], dtype=float)
    atoms_bohr = np.asarray(reference["atoms_bohr"], dtype=float)
    atomic_numbers = np.asarray(reference["atomic_numbers"], dtype=int)
    if int(reference["natoms"]) != 44:
        raise ValueError(f"Expected 44 QM atoms, found {reference['natoms']}")
    physical_atoms_bohr = atoms_bohr[:41]
    dvol_bohr3 = abs(float(np.linalg.det(axes)))

    delta_number_density = np.zeros(tuple(shape), dtype=np.float64)
    pair_audit = []
    for hole, particle, occupation in selected_pairs:
        hole_cube = read_cube(args.cube_dir / f"{args.prefix}.mo{hole}a.cube")
        particle_cube = read_cube(args.cube_dir / f"{args.prefix}.mo{particle}a.cube")
        for cube in (hole_cube, particle_cube):
            if not np.array_equal(cube["shape"], shape):
                raise ValueError("NTO cubes do not share a grid")
            if not np.allclose(cube["origin"], origin) or not np.allclose(cube["axes"], axes):
                raise ValueError("NTO cubes do not share an origin and axes")
        hole_sq = np.asarray(hole_cube["values"], dtype=float) ** 2
        particle_sq = np.asarray(particle_cube["values"], dtype=float) ** 2
        hole_norm = float(hole_sq.sum() * dvol_bohr3)
        particle_norm = float(particle_sq.sum() * dvol_bohr3)
        delta_number_density += occupation * (
            particle_sq / particle_norm - hole_sq / hole_norm
        )
        pair_audit.append(
            {
                "hole": hole,
                "particle": particle,
                "occupation": occupation,
                "hole_grid_norm_before": hole_norm,
                "particle_grid_norm_before": particle_norm,
            }
        )

    # Electron gain is negative charge.  Per-pair orbital normalization makes
    # the total charge change neutral to numerical precision.
    delta_charge_voxels = -delta_number_density.reshape(-1) * dvol_bohr3
    indices = np.indices(tuple(shape), dtype=float).reshape(3, -1).T
    points_bohr = origin + indices @ axes

    nearest_physical = np.zeros(points_bohr.shape[0], dtype=np.int16)
    nearest_d2 = np.full(points_bohr.shape[0], np.inf, dtype=np.float64)
    for atom_index, atom_bohr in enumerate(physical_atoms_bohr):
        d2 = np.sum((points_bohr - atom_bohr) ** 2, axis=1)
        update = d2 < nearest_d2
        nearest_d2[update] = d2[update]
        nearest_physical[update] = atom_index
    atom_delta_q = np.bincount(
        nearest_physical,
        weights=delta_charge_voxels,
        minlength=41,
    ).astype(np.float64)

    grid_delta_mu = np.sum(points_bohr * delta_charge_voxels[:, None], axis=0)
    atom_delta_mu = np.sum(physical_atoms_bohr * atom_delta_q[:, None], axis=0)
    partition_delta_mu_before = atom_delta_mu.copy()
    charge_correction = np.zeros_like(atom_delta_q)
    if args.match_grid_dipole:
        center = physical_atoms_bohr.mean(axis=0)
        constraints = np.vstack(
            [np.ones(len(physical_atoms_bohr)), (physical_atoms_bohr - center).T]
        )
        delta_net = float(delta_charge_voxels.sum() - atom_delta_q.sum())
        target = np.concatenate(
            [[delta_net], grid_delta_mu - atom_delta_mu - center * delta_net]
        )
        charge_correction = constraints.T @ np.linalg.solve(
            constraints @ constraints.T, target
        )
        atom_delta_q += charge_correction
        atom_delta_mu = np.sum(physical_atoms_bohr * atom_delta_q[:, None], axis=0)
    cr2_names = read_cr2_names(args.monomer)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        atom_coords_ang=physical_atoms_bohr * BOHR_TO_ANGSTROM,
        atom_delta_q_e=atom_delta_q,
        atomic_numbers=atomic_numbers[:41],
        cr2_atom_names=np.asarray(cr2_names),
        cr2_atom_count=np.asarray(29),
        physical_qm_atom_count=np.asarray(41),
        link_atom_count=np.asarray(3),
        delta_mu_grid_e_bohr=grid_delta_mu,
        delta_mu_atom_e_bohr=atom_delta_mu,
        nto_occupation_represented=np.asarray(sum(row[2] for row in selected_pairs)),
        nto_pair_count=np.asarray(len(selected_pairs)),
        source_grid_shape=shape,
        source_grid_origin_bohr=origin,
        source_grid_axes_bohr=axes,
        partition=np.asarray("nearest physical QM atom Voronoi partition"),
        grid_dipole_matched=np.asarray(args.match_grid_dipole),
    )

    summary = {
        "output": args.out.name,
        "grid_shape": shape.tolist(),
        "grid_voxels": int(points_bohr.shape[0]),
        "voxel_volume_bohr3": dvol_bohr3,
        "nto_pair_count": len(selected_pairs),
        "nto_occupation_represented": float(sum(row[2] for row in selected_pairs)),
        "grid_net_delta_charge_e": float(delta_charge_voxels.sum()),
        "atom_net_delta_charge_e": float(atom_delta_q.sum()),
        "grid_delta_mu_e_bohr": grid_delta_mu.tolist(),
        "atom_delta_mu_e_bohr": atom_delta_mu.tolist(),
        "atom_delta_mu_before_correction_e_bohr": partition_delta_mu_before.tolist(),
        "delta_mu_partition_error_e_bohr": float(np.linalg.norm(atom_delta_mu - grid_delta_mu)),
        "grid_dipole_matched": bool(args.match_grid_dipole),
        "charge_correction_rms_e": float(np.sqrt(np.mean(charge_correction**2))),
        "charge_correction_max_abs_e": float(np.max(np.abs(charge_correction))),
        "atom_delta_q_min_e": float(atom_delta_q.min()),
        "atom_delta_q_max_e": float(atom_delta_q.max()),
        "pairs": pair_audit,
    }
    summary_path = args.summary or args.out.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

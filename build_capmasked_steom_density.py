#!/usr/bin/env python3
"""Reconstruct the Windows STEOM S1 transition density from its NTO cubes.

All printed S1 NTO pairs are included.  Grid points whose nearest QM atom is
one of the three link caps (one-based atoms 42--44) are removed before the
remaining density is scaled to the ORCA right transition dipole.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


BOHR_A = 0.52917721067
PAIRS = [
    (92, 93, 0.96558779),
    (91, 94, 0.01476372),
    (90, 95, 0.00315100),
    (89, 96, 0.00237768),
    (88, 97, 0.00196516),
    (87, 98, 0.00149604),
    (86, 99, 0.00111854),
]
MU_RIGHT_AU = np.array([1.13353, -2.06197, 3.12464])


def read_cube(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        title = handle.readline().rstrip("\n")
        subtitle = handle.readline().rstrip("\n")
        fields = handle.readline().split()
        signed_natoms = int(fields[0])
        natoms = abs(signed_natoms)
        origin = np.array([float(value) for value in fields[1:4]])
        shape = []
        axes = []
        for _ in range(3):
            fields = handle.readline().split()
            shape.append(int(fields[0]))
            axes.append([float(value) for value in fields[1:4]])
        atoms = []
        for _ in range(natoms):
            fields = handle.readline().split()
            atoms.append([float(value) for value in fields[2:5]])
        if signed_natoms < 0:
            handle.readline()
        values = np.fromstring(" ".join(handle.read().split()), sep=" ")
    shape = np.asarray(shape, dtype=int)
    if values.size != int(np.prod(shape)):
        raise ValueError(f"Cube size mismatch in {path}: {values.size} != {tuple(shape)}")
    return {
        "title": title,
        "subtitle": subtitle,
        "origin": origin,
        "shape": shape,
        "axes": np.asarray(axes, dtype=float),
        "atoms_bohr": np.asarray(atoms, dtype=float),
        "values": values.reshape(tuple(shape)),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cube-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="steom_phenol_svpd.s1")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=1.0e-6)
    args = parser.parse_args(argv)

    reference = read_cube(args.cube_dir / f"{args.prefix}.mo92a.cube")
    shape = reference["shape"]
    axes = reference["axes"]
    origin = reference["origin"]
    dvol = abs(float(np.linalg.det(axes)))
    indices = np.indices(tuple(shape), dtype=float).reshape(3, -1).T
    points_bohr = origin + indices @ axes
    atoms_bohr = reference["atoms_bohr"]
    if atoms_bohr.shape[0] != 44:
        raise ValueError(f"Expected 44 QM atoms, found {atoms_bohr.shape[0]}")

    rho_total = np.zeros(points_bohr.shape[0], dtype=float)
    pair_rows = []
    for hole, particle, occupation in PAIRS:
        hole_cube = read_cube(args.cube_dir / f"{args.prefix}.mo{hole}a.cube")
        particle_cube = read_cube(args.cube_dir / f"{args.prefix}.mo{particle}a.cube")
        for cube in (hole_cube, particle_cube):
            if not np.array_equal(cube["shape"], shape):
                raise ValueError("NTO cubes do not share a grid")
            if not np.allclose(cube["origin"], origin) or not np.allclose(cube["axes"], axes):
                raise ValueError("NTO cube origins/axes do not match")
        rho_pair = (hole_cube["values"] * particle_cube["values"]).reshape(-1)
        q_pair = rho_pair * dvol
        mu_pair = -(points_bohr * q_pair[:, None]).sum(axis=0)
        sign = 1.0 if float(np.dot(mu_pair, MU_RIGHT_AU)) >= 0.0 else -1.0
        coefficient = sign * np.sqrt(occupation)
        rho_total += coefficient * rho_pair
        pair_rows.append((hole, particle, occupation, coefficient, *mu_pair))

    q_unmasked = rho_total * dvol
    mu_unmasked = -(points_bohr * q_unmasked[:, None]).sum(axis=0)

    # Voronoi assignment of every voxel to its nearest QM atom.  Removing the
    # cells assigned to atoms 42--44 excludes cap-associated density without an
    # arbitrary spherical cutoff or changing the 41 physical atoms.
    min_physical_d2 = np.full(points_bohr.shape[0], np.inf)
    for atom in atoms_bohr[:41]:
        min_physical_d2 = np.minimum(
            min_physical_d2, ((points_bohr - atom) ** 2).sum(axis=1)
        )
    min_link_d2 = np.full(points_bohr.shape[0], np.inf)
    for atom in atoms_bohr[41:44]:
        min_link_d2 = np.minimum(min_link_d2, ((points_bohr - atom) ** 2).sum(axis=1))
    physical_mask = min_physical_d2 <= min_link_d2
    q_masked = q_unmasked.copy()
    q_masked[~physical_mask] = 0.0
    mu_masked = -(points_bohr * q_masked[:, None]).sum(axis=0)

    scale = float(np.dot(MU_RIGHT_AU, mu_masked) / np.dot(mu_masked, mu_masked))
    q_scaled = q_masked * scale
    mu_scaled = -(points_bohr * q_scaled[:, None]).sum(axis=0)
    nonzero = np.abs(q_scaled) > args.threshold * np.max(np.abs(q_scaled))
    points_ang = points_bohr[nonzero] * BOHR_A
    charges = q_scaled[nonzero]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        pts_ang=points_ang,
        q=charges,
        mu_au=mu_scaled,
        mu_target_au=MU_RIGHT_AU,
        mu_unmasked_au=mu_unmasked,
        mu_masked_unscaled_au=mu_masked,
        nto_pairs=np.asarray(pair_rows, dtype=float),
        physical_qm_atom_count=np.asarray(41),
        excluded_link_atom_indices_one_based=np.asarray([42, 43, 44]),
        mask_method=np.asarray("nearest-QM-atom Voronoi assignment"),
        retained_grid_fraction=np.asarray(float(physical_mask.mean())),
        retained_point_count=np.asarray(int(nonzero.sum())),
        source_grid_shape=shape,
        source_grid_origin_bohr=origin,
        source_grid_axes_bohr=axes,
    )

    print(f"S1 NTO occupation represented: {sum(row[2] for row in PAIRS):.8f}")
    print(f"Unmasked mu (au): {mu_unmasked} |mu|={np.linalg.norm(mu_unmasked):.6f}")
    print(f"Cap-masked mu (au): {mu_masked} |mu|={np.linalg.norm(mu_masked):.6f}")
    print(f"Spectroscopic scale: {scale:.8f}")
    print(f"Scaled mu (au): {mu_scaled} |mu|={np.linalg.norm(mu_scaled):.6f}")
    print(f"Target mu (au): {MU_RIGHT_AU} |mu|={np.linalg.norm(MU_RIGHT_AU):.6f}")
    print(
        f"Retained {nonzero.sum()} density points; cap Voronoi mask retained "
        f"{physical_mask.mean():.6%} of the source grid"
    )
    print(f"Net retained transition charge: {charges.sum():+.6e}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

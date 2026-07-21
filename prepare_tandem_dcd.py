#!/usr/bin/env python3
"""Make the protein-only tandem DCD PBC-whole along its covalent residue order.

The solvated NVT trajectory legitimately wraps atoms into the primary unit cell.
For a long, single-chain tandem this can put a chromophore in a different image
from its own beta barrel.  Each residue is first made internally whole and is
then translated by whole box vectors to sit next to the preceding residue.
Coordinates are never fitted or otherwise deformed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from mdtraj.formats import DCDTrajectoryFile
from openmm import unit
from openmm.app import PDBFile


def _nearest_image_shift(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Return a lattice-vector shift that puts *delta* in its nearest image."""
    fractional = _right_multiply(np.asarray(delta), _inverse3(box))
    return _right_multiply(-np.rint(fractional), box)


def _right_multiply(values: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Compute values @ matrix explicitly, avoiding BLAS overhead for 3 columns."""
    v = np.asarray(values)
    return np.stack(
        (
            v[..., 0] * matrix[0, 0] + v[..., 1] * matrix[1, 0] + v[..., 2] * matrix[2, 0],
            v[..., 0] * matrix[0, 1] + v[..., 1] * matrix[1, 1] + v[..., 2] * matrix[2, 1],
            v[..., 0] * matrix[0, 2] + v[..., 1] * matrix[1, 2] + v[..., 2] * matrix[2, 2],
        ),
        axis=-1,
    )


def _inverse3(matrix: np.ndarray) -> np.ndarray:
    """Small dependency-free 3x3 inverse (avoids a BLAS call per residue)."""
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    cof = np.array(
        [
            [e * i - f * h, c * h - b * i, b * f - c * e],
            [f * g - d * i, a * i - c * g, c * d - a * f],
            [d * h - e * g, b * g - a * h, a * e - b * d],
        ]
    )
    determinant = a * cof[0, 0] + b * cof[1, 0] + c * cof[2, 0]
    if abs(determinant) < 1.0e-12:
        raise ValueError("singular periodic box")
    return cof / determinant


def _box_from_lengths_angles(lengths: np.ndarray, angles_deg: np.ndarray) -> np.ndarray:
    """Build row-vector triclinic boxes from DCD lengths/angles (all in Angstrom)."""
    a, b, c = map(float, lengths)
    alpha, beta, gamma = np.deg2rad(angles_deg)
    cos_a, cos_b, cos_g = np.cos(alpha), np.cos(beta), np.cos(gamma)
    sin_g = np.sin(gamma)
    av = np.array([a, 0.0, 0.0])
    bv = np.array([b * cos_g, b * sin_g, 0.0])
    cx = c * cos_b
    cy = c * (cos_a - cos_b * cos_g) / sin_g
    cz = np.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.vstack((av, bv, np.array([cx, cy, cz])))


def unwrap_frame(xyz: np.ndarray, box: np.ndarray, residue_atoms: list[np.ndarray]) -> np.ndarray:
    out = np.array(xyz, dtype=np.float64, copy=True)
    previous_center = None
    for atom_indices in residue_atoms:
        coords = out[atom_indices]

        # A residue (especially CR2) can itself straddle a periodic boundary.
        anchor = coords[0]
        deltas = coords - anchor
        fractional = _right_multiply(deltas, _inverse3(box))
        coords -= _right_multiply(np.rint(fractional), box)

        center = coords.mean(axis=0)
        if previous_center is not None:
            shift = _nearest_image_shift(center - previous_center, box)
            coords += shift
            center += shift

        out[atom_indices] = coords
        previous_center = center
    return out


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traj", type=Path, required=True)
    ap.add_argument("--topology", type=Path, required=True)
    ap.add_argument("--out-dcd", type=Path, required=True)
    ap.add_argument("--out-first-pdb", type=Path, default=None)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    # Use MDTraj's low-level DCD reader.  It retains native Angstrom units and
    # avoids its high-level PDB loader, which is unnecessary for this operation.
    with DCDTrajectoryFile(str(args.traj), "r") as dcd:
        xyz, cell_lengths, cell_angles = dcd.read()
    if cell_lengths is None or cell_angles is None:
        raise ValueError("DCD has no unit-cell vectors; periodic unwrapping is undefined")

    pdb = PDBFile(str(args.topology))
    residue_atoms = [np.array([a.index for a in r.atoms()], dtype=int) for r in pdb.topology.residues()]
    if xyz.shape[1] != pdb.topology.getNumAtoms():
        raise ValueError(f"DCD/PDB atom mismatch: {xyz.shape[1]} versus {pdb.topology.getNumAtoms()}")

    corrected = np.empty_like(xyz)
    for frame in range(xyz.shape[0]):
        box = _box_from_lengths_angles(cell_lengths[frame], cell_angles[frame])
        corrected[frame] = unwrap_frame(xyz[frame], box, residue_atoms)

    args.out_dcd.parent.mkdir(parents=True, exist_ok=True)
    with DCDTrajectoryFile(str(args.out_dcd), "w") as dcd:
        dcd.write(corrected, cell_lengths=cell_lengths, cell_angles=cell_angles)
    if args.out_first_pdb is not None:
        args.out_first_pdb.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_first_pdb, "w") as handle:
            PDBFile.writeFile(pdb.topology, corrected[0] * unit.angstrom, handle, keepIds=True)

    print(
        f"[unwrap] wrote {args.out_dcd}: {xyz.shape[0]} frames, "
        f"{xyz.shape[1]} atoms, {len(residue_atoms)} sequentially unwrapped residues"
    )


if __name__ == "__main__":
    main()

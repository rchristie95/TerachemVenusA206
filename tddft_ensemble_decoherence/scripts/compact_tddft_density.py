#!/usr/bin/env python3
"""Cap-mask and compact a TeraChem transition density for tracking and TDC.

The 44-atom validated QM region contains 41 physical atoms followed by three
link hydrogens. Grid points assigned by nearest-atom Voronoi partition to a
link atom are excluded. The retained physical density is phase-aligned to a
reference atomic transition-charge vector when supplied; otherwise it is
phase-aligned to the printed TeraChem transition dipole. The cap-masked density
is scaled only to restore that same calculation's printed dipole magnitude.
All phase and normalization diagnostics are retained in the output NPZ.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
import coupling_core as cc  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(errors="replace").splitlines()
    count = int(lines[0].split()[0])
    rows = [line.split() for line in lines[2 : 2 + count]]
    if len(rows) != count:
        raise ValueError(f"Incomplete XYZ: {path}")
    return [row[0] for row in rows], np.asarray([[float(x) for x in row[1:4]] for row in rows])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dx", type=Path, required=True)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--root", type=int, required=True)
    parser.add_argument("--energy-ev", type=float, required=True)
    parser.add_argument("--oscillator-strength", type=float, required=True)
    parser.add_argument("--target-mu-au", nargs=3, type=float, required=True)
    parser.add_argument("--reference", type=Path, help="Reference compact NPZ for atomic-overlap phase alignment")
    parser.add_argument("--physical-atoms", type=int, default=41)
    parser.add_argument("--absolute-threshold", type=float, default=1.0e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    symbols, atoms = read_xyz(args.xyz)
    if len(atoms) != 44 or args.physical_atoms != 41:
        raise ValueError(f"Validated mapping requires 44 total/41 physical atoms; got {len(atoms)}/{args.physical_atoms}")

    points, charges = cc.read_dx(args.dx, threshold=args.absolute_threshold, stride=1)
    if not len(charges):
        raise ValueError("Transition density is empty after thresholding")
    nearest_physical_distance, nearest_physical = cKDTree(atoms[: args.physical_atoms]).query(points)
    nearest_link_distance, _ = cKDTree(atoms[args.physical_atoms :]).query(points)
    keep = nearest_physical_distance <= nearest_link_distance
    points = np.asarray(points[keep], dtype=np.float64)
    charges = np.asarray(charges[keep], dtype=np.float64)
    nearest_physical = np.asarray(nearest_physical[keep], dtype=np.int32)

    # A transition density is exactly neutral. Voronoi removal of link-cap
    # regions and numerical thresholding can break that identity. Restore it
    # with the least invasive sign-wise normalization: positive and negative
    # lobes are scaled to their mean absolute charge before spectroscopic
    # dipole normalization. Retain the pre-correction error and both factors.
    net_charge_before_neutralization = float(charges.sum())
    positive_sum = float(charges[charges > 0.0].sum())
    negative_magnitude = float(-charges[charges < 0.0].sum())
    if positive_sum <= 0.0 or negative_magnitude <= 0.0:
        raise ValueError("Transition density does not contain both phases")
    mean_lobe_charge = 0.5 * (positive_sum + negative_magnitude)
    positive_scale = mean_lobe_charge / positive_sum
    negative_scale = mean_lobe_charge / negative_magnitude
    charges[charges > 0.0] *= positive_scale
    charges[charges < 0.0] *= negative_scale

    atom_q_raw = np.bincount(nearest_physical, weights=charges, minlength=args.physical_atoms)
    origin = atoms[: args.physical_atoms].mean(axis=0)
    mu_raw = cc.transition_dipole_au(points, charges, origin_angstrom=origin)
    target_mu = np.asarray(args.target_mu_au, dtype=np.float64)

    overlap_signed = np.nan
    overlap_score = np.nan
    phase_source = "printed_transition_dipole"
    phase = 1.0 if float(np.dot(mu_raw, target_mu)) >= 0.0 else -1.0
    if args.reference is not None:
        reference = np.load(args.reference)
        reference_q = np.asarray(reference["atom_transition_q"], dtype=np.float64)
        if reference_q.shape != atom_q_raw.shape:
            raise ValueError("Reference physical-atom transition vector has the wrong shape")
        denominator = float(np.linalg.norm(atom_q_raw) * np.linalg.norm(reference_q))
        overlap_signed = float(np.dot(atom_q_raw, reference_q) / denominator) if denominator else np.nan
        overlap_score = abs(overlap_signed)
        phase = 1.0 if overlap_signed >= 0.0 else -1.0
        phase_source = "physical_atom_transition_density_overlap"

    charges *= phase
    atom_q_raw *= phase
    mu_phase = mu_raw * phase
    denominator = float(np.dot(mu_phase, mu_phase))
    if args.reference is not None:
        # TeraChem's transition-density and printed transition-dipole signs are
        # each arbitrary.  Reference overlap defines the physical phase across
        # jobs; only the positive spectroscopic magnitude is normalized here.
        scale = float(np.linalg.norm(target_mu) / np.linalg.norm(mu_phase)) if denominator else 1.0
    else:
        scale = float(np.dot(target_mu, mu_phase) / denominator) if denominator else 1.0
    if scale < 0.0:
        raise ValueError("Phase alignment left a negative spectroscopic normalization")
    charges *= scale
    atom_q = atom_q_raw * scale
    mu_scaled = mu_phase * scale

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        pts_ang=points,
        q=charges,
        atom_transition_q=atom_q,
        atom_symbols=np.asarray(symbols[: args.physical_atoms]),
        physical_atom_coords_ang=atoms[: args.physical_atoms],
        excluded_link_atom_coords_ang=atoms[args.physical_atoms :],
        root=np.asarray(args.root),
        energy_eV=np.asarray(args.energy_ev),
        oscillator_strength=np.asarray(args.oscillator_strength),
        target_mu_au=target_mu,
        masked_mu_before_phase_au=mu_raw,
        final_mu_au=mu_scaled,
        phase=np.asarray(phase),
        phase_source=np.asarray(phase_source),
        assignment_overlap_signed=np.asarray(overlap_signed),
        assignment_score=np.asarray(overlap_score),
        spectroscopic_scale=np.asarray(scale),
        retained_point_count=np.asarray(len(charges)),
        retained_grid_charge=np.asarray(float(charges.sum())),
        net_charge_before_neutralization=np.asarray(net_charge_before_neutralization),
        positive_lobe_neutralization_scale=np.asarray(positive_scale),
        negative_lobe_neutralization_scale=np.asarray(negative_scale),
        physical_atom_count=np.asarray(args.physical_atoms),
        excluded_link_atom_indices_one_based=np.arange(args.physical_atoms + 1, len(atoms) + 1),
        mask_method=np.asarray("nearest-QM-atom Voronoi assignment"),
        source_dx_sha256=np.asarray(sha256(args.dx)),
        source_xyz_sha256=np.asarray(sha256(args.xyz)),
    )
    print(
        f"root={args.root} points={len(charges)} phase={phase:+.0f} scale={scale:.8f} "
        f"overlap={overlap_score} net_q={charges.sum():+.3e} "
        f"raw_net_q={net_charge_before_neutralization:+.3e} "
        f"|mu|={np.linalg.norm(mu_scaled):.6f} au"
    )


if __name__ == "__main__":
    main()

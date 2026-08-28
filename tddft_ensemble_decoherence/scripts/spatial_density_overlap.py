#!/usr/bin/env python3
"""Bidirectional spatial overlap of aligned compact transition densities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target_center - source_center @ rotation.T
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def interpolate(source_points: np.ndarray, source_q: np.ndarray, target_points: np.ndarray, sigma: float, cutoff: float) -> tuple[np.ndarray, float]:
    tree = cKDTree(source_points)
    distance, index = tree.query(target_points, k=4, workers=-1)
    weight = np.exp(-0.5 * (distance / sigma) ** 2)
    weight[distance > cutoff] = 0.0
    denominator = weight.sum(axis=1)
    values = np.zeros(len(target_points), dtype=float)
    valid = denominator > 0.0
    values[valid] = np.sum(weight[valid] * source_q[index[valid]], axis=1) / denominator[valid]
    return values, float(np.mean(valid))


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator else float("nan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sigma", type=float, default=0.12)
    parser.add_argument("--cutoff", type=float, default=0.35)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reference = np.load(args.reference)
    candidate = np.load(args.candidate)
    ref_atoms = np.asarray(reference["physical_atom_coords_ang"], float)
    cand_atoms = np.asarray(candidate["physical_atom_coords_ang"], float)
    ref_symbols = np.asarray(reference["atom_symbols"]).astype(str)
    cand_symbols = np.asarray(candidate["atom_symbols"]).astype(str)
    heavy = np.asarray([index for index in range(29) if ref_symbols[index].upper() != "H"])
    if not np.array_equal(ref_symbols[:29], cand_symbols[:29]) or len(heavy) != 19:
        raise ValueError("Reference/candidate CR2 mappings are not the validated 19-heavy-atom mapping")
    rotation, translation, fit_rmsd = kabsch(cand_atoms[heavy], ref_atoms[heavy])
    ref_points = np.asarray(reference["pts_ang"], float)
    ref_q = np.asarray(reference["q"], float)
    cand_points = np.asarray(candidate["pts_ang"], float) @ rotation.T + translation
    cand_q = np.asarray(candidate["q"], float)
    cand_on_ref, fraction_ref = interpolate(cand_points, cand_q, ref_points, args.sigma, args.cutoff)
    ref_on_cand, fraction_cand = interpolate(ref_points, ref_q, cand_points, args.sigma, args.cutoff)
    forward = cosine(ref_q, cand_on_ref)
    reverse = cosine(cand_q, ref_on_cand)
    result = {
        "reference": str(args.reference.resolve()),
        "candidate": str(args.candidate.resolve()),
        "root": int(candidate["root"]),
        "energy_eV": float(candidate["energy_eV"]),
        "oscillator_strength": float(candidate["oscillator_strength"]),
        "heavy_atom_fit_rmsd_A": fit_rmsd,
        "sigma_A": args.sigma,
        "cutoff_A": args.cutoff,
        "reference_points_covered_fraction": fraction_ref,
        "candidate_points_covered_fraction": fraction_cand,
        "spatial_overlap_forward": forward,
        "spatial_overlap_reverse": reverse,
        "spatial_overlap_mean": 0.5 * (forward + reverse),
        "atomic_assignment_score": float(candidate["assignment_score"]),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

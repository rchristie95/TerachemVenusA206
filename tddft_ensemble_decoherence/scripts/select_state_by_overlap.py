#!/usr/bin/env python3
"""Select a physical state by spatial-density overlap with an ambiguity guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from spatial_density_overlap import cosine, interpolate, kabsch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidates", nargs="+", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sigma", type=float, default=0.12)
    parser.add_argument("--cutoff", type=float, default=0.35)
    parser.add_argument("--minimum-margin", type=float, default=0.10)
    return parser.parse_args()


def overlap(reference, candidate, sigma: float, cutoff: float) -> dict:
    ref_atoms = np.asarray(reference["physical_atom_coords_ang"], float)
    cand_atoms = np.asarray(candidate["physical_atom_coords_ang"], float)
    symbols = np.asarray(reference["atom_symbols"]).astype(str)
    heavy = np.asarray([index for index in range(29) if symbols[index].upper() != "H"])
    rotation, translation, fit_rmsd = kabsch(cand_atoms[heavy], ref_atoms[heavy])
    ref_points, ref_q = np.asarray(reference["pts_ang"], float), np.asarray(reference["q"], float)
    cand_points = np.asarray(candidate["pts_ang"], float) @ rotation.T + translation
    cand_q = np.asarray(candidate["q"], float)
    cand_on_ref, ref_coverage = interpolate(cand_points, cand_q, ref_points, sigma, cutoff)
    ref_on_cand, cand_coverage = interpolate(ref_points, ref_q, cand_points, sigma, cutoff)
    forward, reverse = cosine(ref_q, cand_on_ref), cosine(cand_q, ref_on_cand)
    return {
        "root": int(candidate["root"]),
        "energy_eV": float(candidate["energy_eV"]),
        "oscillator_strength": float(candidate["oscillator_strength"]),
        "spatial_overlap_forward": forward,
        "spatial_overlap_reverse": reverse,
        "spatial_overlap_mean": 0.5 * (forward + reverse),
        "spatial_overlap_absolute": abs(0.5 * (forward + reverse)),
        "atomic_assignment_score": float(candidate["assignment_score"]),
        "heavy_atom_fit_rmsd_A": fit_rmsd,
        "reference_points_covered_fraction": ref_coverage,
        "candidate_points_covered_fraction": cand_coverage,
    }


def main() -> None:
    args = parse_args()
    reference = np.load(args.reference)
    candidates = [overlap(reference, np.load(path), args.sigma, args.cutoff) for path in args.candidates]
    candidates.sort(key=lambda item: item["spatial_overlap_absolute"], reverse=True)
    margin = (
        candidates[0]["spatial_overlap_absolute"] - candidates[1]["spatial_overlap_absolute"]
        if len(candidates) > 1 else float("inf")
    )
    result = {
        "status": "selected" if margin >= args.minimum_margin else "ambiguous_overlap_margin",
        "selected_root": candidates[0]["root"],
        "selected_energy_eV": candidates[0]["energy_eV"],
        "selected_oscillator_strength": candidates[0]["oscillator_strength"],
        "overlap_margin": margin,
        "minimum_margin": args.minimum_margin,
        "sigma_A": args.sigma,
        "cutoff_A": args.cutoff,
        "candidates": candidates,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["status"] != "selected":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

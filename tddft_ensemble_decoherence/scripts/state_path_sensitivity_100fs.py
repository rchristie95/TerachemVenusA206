#!/usr/bin/env python3
"""Quantify decoherence sensitivity to globally plausible state paths."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from summarize_100fs50 import cumulant


ROOT = Path(__file__).resolve().parents[1]
WINDOW = "short100fs50_fixedmm"
OUTDIR = ROOT / "rapid_100fs50_fixedmm"


def root_energies(site: str) -> list[dict[int, float]]:
    values = []
    for frame in range(50):
        path = (
            ROOT
            / f"frames/{WINDOW}/frame_{frame:04d}/site_{site}/density_summary.json"
        )
        parsed = json.loads(path.read_text())
        values.append(
            {int(root["root"]): float(root["energy_eV"]) for root in parsed["roots"]}
        )
    return values


def evaluate(
    a_roots: list[int],
    b_roots: list[int],
    a_values: list[dict[int, float]],
    b_values: list[dict[int, float]],
) -> dict:
    ea = np.asarray([a_values[i][root] for i, root in enumerate(a_roots)])
    eb = np.asarray([b_values[i][root] for i, root in enumerate(b_roots)])
    gap = ea - eb
    result = cumulant(gap, 2.0)
    return {
        "A_roots": a_roots,
        "B_roots": b_roots,
        "T2star_fs": result["T2star_fs"],
        "gap_mean_eV": float(np.mean(gap)),
        "gap_sample_sd_eV": float(np.std(gap, ddof=1)),
        "maximum_adjacent_gap_jump_eV": float(np.max(np.abs(np.diff(gap)))),
    }


def main() -> None:
    audit = json.loads((OUTDIR / "global_state_tracking.json").read_text())
    a_values, b_values = root_energies("A"), root_energies("B")
    canonical_a = [int(root) for root in audit["sites"]["A"]["global_roots"]]
    canonical_b = [int(root) for root in audit["sites"]["B"]["global_roots"]]
    variants = []
    for index, item in enumerate(audit["sensitivity"]["B"]):
        record = evaluate(
            canonical_a,
            [int(root) for root in item["roots"]],
            a_values,
            b_values,
        )
        record.update(
            {
                "variant_index": index,
                "density_weight": item["density_weight"],
                "energy_weight": item["energy_weight"],
            }
        )
        variants.append(record)
    finite = np.asarray(
        [item["T2star_fs"] for item in variants if item["T2star_fs"] is not None]
    )
    result = {
        "status": "complete",
        "method": (
            "Twelve predeclared Viterbi weight combinations using exact adjacent "
            "B-site spatial transition-density matrices"
        ),
        "canonical": evaluate(
            canonical_a, canonical_b, a_values, b_values
        ),
        "variants": variants,
        "finite_variant_count": int(len(finite)),
        "censored_variant_count": int(len(variants) - len(finite)),
        "finite_T2star_range_fs": (
            [float(np.min(finite)), float(np.max(finite))]
            if len(finite)
            else [None, None]
        ),
        "finite_T2star_median_fs": (
            float(np.median(finite)) if len(finite) else None
        ),
    }
    (OUTDIR / "state_path_sensitivity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

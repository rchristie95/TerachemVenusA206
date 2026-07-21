#!/usr/bin/env python3
"""Validate timescales extracted from high-cadence tandem-Venus gap traces.

This script operates only on the compact ``gap_timeseries.npz`` files written
by ``analyze_solvation_decoherence.py``.  It joins contiguous segments, repeats
the analysis for trajectory blocks and sampling strides, and compares the
ion-inclusive result with a protein+water trace.  The latter is useful here
because the exact -1 CR2 charge transplant left the pre-solvated box at -2e.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


LIGHT_CM_PER_S = 2.99792458e10
KB_OVER_HC_CM_PER_K = 0.695034800


def portable_path(path: Path) -> str:
    """Keep provenance useful without recording a workstation root."""
    project_root = Path(__file__).resolve().parent
    try:
        return path.resolve().relative_to(project_root).as_posix()
    except ValueError:
        return f"external/{path.name}"


def correlation(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float) - np.mean(values)
    n = len(values)
    return np.correlate(values, values, mode="full")[n - 1 :] / np.arange(n, 0, -1)


def cross_correlation(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    return float(np.mean(a * b) / np.sqrt(np.mean(a * a) * np.mean(b * b)))


def first_crossing(time_fs: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    hits = np.flatnonzero(values <= threshold)
    if not len(hits):
        return None
    index = int(hits[0])
    if index == 0:
        return float(time_fs[0])
    x0, x1 = time_fs[index - 1 : index + 1]
    y0, y1 = values[index - 1 : index + 1]
    if y1 == y0:
        return float(x1)
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def positive_integral(time_fs: np.ndarray, normalized: np.ndarray) -> float:
    zeros = np.flatnonzero(normalized <= 0.0)
    stop = int(zeros[0]) + 1 if len(zeros) else len(normalized)
    return float(np.trapezoid(normalized[:stop], time_fs[:stop]))


def cumulative_trapezoid(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros(len(values), dtype=float)
    if len(values) > 1:
        out[1:] = np.cumsum(0.5 * (values[1:] + values[:-1]) * np.diff(x))
    return out


def metrics(a: np.ndarray, b: np.ndarray, dt_fs: float, temperature_k: float) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    difference = a - b
    time_fs = np.arange(len(a), dtype=float) * dt_fs
    ca = correlation(a)
    cb = correlation(b)
    cd = correlation(difference)
    sa = ca / ca[0]
    sb = cb / cb[0]
    sd = cd / cd[0]

    time_s = time_fs * 1.0e-15
    c_omega = (2.0 * np.pi * LIGHT_CM_PER_S) ** 2 * cd
    g = cumulative_trapezoid(cumulative_trapezoid(c_omega, time_s), time_s)
    coherence = np.exp(-np.clip(g, 0.0, 700.0))
    t2_fs = first_crossing(time_fs, coherence, 1.0 / np.e)
    sigma_d = float(np.std(difference, ddof=1))
    kbt_cm = KB_OVER_HC_CM_PER_K * temperature_k
    return {
        "frames": int(len(a)),
        "duration_fs": float(time_fs[-1]) if len(time_fs) else 0.0,
        "dt_fs": float(dt_fs),
        "sigma_A_cm-1": float(np.std(a, ddof=1)),
        "sigma_B_cm-1": float(np.std(b, ddof=1)),
        "sigma_difference_cm-1": sigma_d,
        "rho_AB_zero_lag": cross_correlation(a, b),
        "site_A_1e_fs": first_crossing(time_fs, sa, 1.0 / np.e),
        "site_B_1e_fs": first_crossing(time_fs, sb, 1.0 / np.e),
        "difference_1e_fs": first_crossing(time_fs, sd, 1.0 / np.e),
        "site_A_positive_integral_fs": positive_integral(time_fs, sa),
        "site_B_positive_integral_fs": positive_integral(time_fs, sb),
        "difference_positive_integral_fs": positive_integral(time_fs, sd),
        "classical_cumulant_T2_1e_fs": t2_fs,
        "effective_differential_reorganization_cm-1": sigma_d**2 / (2.0 * kbt_cm),
    }


def trace_for_groups(energies: np.ndarray, groups: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    return energies[:, 0, groups].sum(axis=-1), energies[:, 1, groups].sum(axis=-1)


def aggregate(rows: list[dict], keys: list[str]) -> dict:
    result = {"blocks": len(rows)}
    for key in keys:
        values = np.asarray([row[key] for row in rows if row.get(key) is not None], dtype=float)
        result[key] = {
            "median": float(np.median(values)),
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--block-fs", type=float, default=1000.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args()

    loaded = [np.load(path) for path in args.inputs]
    names = [str(item) for item in loaded[0]["group_names"]]
    for item in loaded[1:]:
        if [str(value) for value in item["group_names"]] != names:
            raise RuntimeError("Group names do not match across inputs")
    dt_fs = float(np.median(np.diff(loaded[0]["time_fs"])))
    energies = np.concatenate([np.asarray(item["energies_cm"], dtype=float) for item in loaded], axis=0)
    group_sets = {
        "total_with_ions": tuple(range(len(names))),
        "protein_plus_water": tuple(index for index, name in enumerate(names) if name != "ions"),
        **{name: (index,) for index, name in enumerate(names)},
    }

    full = {}
    for label, groups in group_sets.items():
        a, b = trace_for_groups(energies, groups)
        full[label] = metrics(a, b, dt_fs, args.temperature_k)

    segments = []
    start = 0
    for path, item in zip(args.inputs, loaded):
        count = len(item["time_fs"])
        row = {"input": portable_path(path), "frames": count}
        for label in ("total_with_ions", "protein_plus_water"):
            a, b = trace_for_groups(energies[start : start + count], group_sets[label])
            row[label] = metrics(a, b, dt_fs, args.temperature_k)
        segments.append(row)
        start += count

    block_frames = max(16, int(round(args.block_fs / dt_fs)))
    block_rows = []
    for block_start in range(0, len(energies) - block_frames + 1, block_frames):
        block = energies[block_start : block_start + block_frames]
        row = {
            "block": len(block_rows) + 1,
            "start_fs": float(block_start * dt_fs),
            "stop_fs": float((block_start + block_frames - 1) * dt_fs),
        }
        for label in ("total_with_ions", "protein_plus_water", "protein", "water"):
            a, b = trace_for_groups(block, group_sets[label])
            row[label] = metrics(a, b, dt_fs, args.temperature_k)
        block_rows.append(row)

    validation_keys = [
        "sigma_difference_cm-1",
        "rho_AB_zero_lag",
        "site_A_1e_fs",
        "site_B_1e_fs",
        "difference_1e_fs",
        "classical_cumulant_T2_1e_fs",
    ]
    block_summary = {
        label: aggregate([row[label] for row in block_rows], validation_keys)
        for label in ("total_with_ions", "protein_plus_water", "protein", "water")
    }

    stride_rows = []
    for stride in (1, 2, 4, 8):
        if len(energies) // stride < 100:
            continue
        sampled = energies[::stride]
        row = {"stride": stride, "dt_fs": dt_fs * stride}
        for label in ("total_with_ions", "protein_plus_water"):
            a, b = trace_for_groups(sampled, group_sets[label])
            row[label] = metrics(a, b, dt_fs * stride, args.temperature_k)
        stride_rows.append(row)

    output = {
        "inputs": [portable_path(path) for path in args.inputs],
        "frames": int(len(energies)),
        "dt_fs": dt_fs,
        "duration_fs": float((len(energies) - 1) * dt_fs),
        "temperature_k": args.temperature_k,
        "group_names": names,
        "full_trajectory": full,
        "segments": segments,
        "block_size_fs": args.block_fs,
        "block_metrics": block_rows,
        "block_summary": block_summary,
        "stride_sensitivity": stride_rows,
        "interpretation_guardrails": [
            "protein_plus_water excludes explicit ions but not the PME neutralizing background used in the MD propagation",
            "block estimates remove each block mean and therefore target fast local fluctuations rather than slow conformational drift",
            "the cumulant is classical and uses only the fixed electrostatic STEOM difference-density probe",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "validation_summary.json").write_text(json.dumps(output, indent=2), encoding="utf-8")

    with (args.out / "block_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["block", "start_fs", "stop_fs", "environment", *validation_keys])
        for row in block_rows:
            for label in ("total_with_ions", "protein_plus_water", "protein", "water"):
                writer.writerow([
                    row["block"], row["start_fs"], row["stop_fs"], label,
                    *(row[label][key] for key in validation_keys),
                ])
    print(json.dumps({
        "full_trajectory": {
            key: full[key] for key in ("total_with_ions", "protein_plus_water", "protein", "water")
        },
        "block_summary": block_summary,
        "stride_sensitivity": stride_rows,
    }, indent=2))


if __name__ == "__main__":
    main()

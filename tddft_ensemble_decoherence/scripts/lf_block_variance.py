#!/usr/bin/env python3
"""How much of the gap variance is slow, by block averaging rather than by FFT.

WHY NOT THE SPECTRUM. The stride-20 campaign (100 frames at 80 fs over 7.92 ps)
was launched to resolve the spectral density below v8dense's 33 cm^-1 floor. It
cannot: the gap is not band-limited below the 208 cm^-1 Nyquist that an 80 fs
cadence implies -- 92% of its variance lies above 208 -- so that power folds
into the low-frequency band. Subsampling v8dense itself to 80 fs reproduces the
artefact exactly, inflating the 33-208 cm^-1 spread from 136 to ~430 cm^-1 on
identical underlying data. Any spectral density read off this record below
208 cm^-1 is aliased high-frequency power. Resolving 4 cm^-1 without aliasing
needs 8 ps at the ORIGINAL cadence, about 2000 frames, which is ~9 days.

WHAT THE RECORD CAN STILL ANSWER. Block averaging is immune to the aliasing,
because aliasing only relabels which frequency the fast power appears at -- it
does not make consecutive samples correlated. With tau_c ~ 5 fs and a sampling
interval of 80 fs, the fast component of successive samples is independent, so

    Var(mean of M samples) = sigma_slow^2(T = M dt) + sigma_fast^2 / M

and the question "is there variance slower than a picosecond?" is answered by
whether that curve plateaus above zero or decays as 1/M. That is exactly the
question the cumulant left open, and it needs the 8 ps SPAN rather than a
resolved spectrum.

Read-only apart from results/lf_block_variance.json.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def series(files, key):
    rows = [r for f in files for r in csv.DictReader((ROOT / "manifests" / f).open())
            if r["convergence_status"] == "complete" and r[key]]
    rows.sort(key=lambda r: int(r["frame_index"]))
    return np.array([float(r[key]) for r in rows])


def block_variance(x, m):
    """Variance of the means of consecutive blocks of m samples, and its se."""
    n = x.size // m
    if n < 4:
        return None
    means = x[: n * m].reshape(n, m).mean(axis=1)
    v = float(means.var(ddof=1))
    return v, float(v * np.sqrt(2.0 / (n - 1))), n


def curve(x, dt_fs, label):
    print(f"  --- {label}: {x.size} samples at {dt_fs:g} fs, "
          f"span {x.size * dt_fs / 1000:.2f} ps, sd {x.std(ddof=1):.0f} cm^-1")
    rows = []
    print(f"      {'block':>8s} {'blocks':>7s} {'sd of block means':>19s}"
          f" {'sigma_fast/sqrt(M) alone':>26s}")
    sd_tot = x.std(ddof=1)
    for m in (1, 2, 4, 8, 16, 32, 64, 128):
        r = block_variance(x, m)
        if r is None:
            continue
        v, se, n = r
        pure_fast = sd_tot / np.sqrt(m)      # what an all-fast signal would give
        rows.append(dict(block_samples=m, block_fs=m * dt_fs, n_blocks=n,
                         sd_block_mean_cm=float(np.sqrt(v)),
                         se_variance=se,
                         all_fast_prediction_cm=float(pure_fast)))
        print(f"      {m*dt_fs:7.0f}fs {n:7d} {np.sqrt(v):15.0f} cm^-1"
              f" {pure_fast:22.0f} cm^-1")
    return rows


def main() -> None:
    out = {}
    print("=== variance of block means: does anything survive averaging? ===")
    print("    If all the disorder were fast, the block-mean spread would fall")
    print("    as sigma/sqrt(M) (last column). Any excess over that is slow.\n")
    for site in ("A", "B"):
        key = f"site_{site}_energy_cm-1"
        out[f"site_{site}"] = {
            "lf_8ps": curve(series([f"v4full_lf{w}_frames.csv" for w in "abc"], key),
                            80.0, f"site {site}  8 ps record (NEW)"),
            "v8dense_1ps": curve(series([f"v8dense_{w}_frames.csv" for w in "abc"], key),
                                 2.0, f"site {site}  1 ps record"),
            "v5camp_1ns": curve(series(["v5camp_frames.csv"], key),
                                2000.0, f"site {site}  1 ns record"),
        }
        print()
    (ROOT / "results" / "lf_block_variance.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print("  wrote results/lf_block_variance.json")


if __name__ == "__main__":
    main()

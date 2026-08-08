#!/usr/bin/env python3
"""Compare a regenerated per-frame coupling series with the retained manuscript series."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def load_csv(path: Path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def stats(values):
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(len(arr)),
        "mean": float(arr.mean()),
        "std_sample": float(arr.std(ddof=1)),
        "min": float(arr.min()),
        "max": float(arr.max()),
        "median": float(np.median(arr)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manuscript", type=Path, required=True)
    ap.add_argument("--regenerated", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    old, new = load_csv(args.manuscript), load_csv(args.regenerated)
    if len(old) != len(new):
        raise ValueError(f"Frame count mismatch: {len(old)} versus {len(new)}")
    fields = ["J_cm", "J_pda_cm", "angle_deg", "separation_A", "aln_A_rms", "aln_B_rms"]
    combined = []
    for old_row, new_row in zip(old, new):
        if int(old_row["frame"]) != int(new_row["frame"]):
            raise ValueError("Frame numbering mismatch")
        row = {"frame": int(old_row["frame"])}
        for field in fields:
            row[f"manuscript_{field}"] = float(old_row[field])
            row[f"regenerated_{field}"] = float(new_row[field])
        combined.append(row)

    out_fields = ["frame"] + [f"{prefix}_{field}" for prefix in ("manuscript", "regenerated") for field in fields]
    with open(args.out / "coupling_frame_comparison.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(combined)

    summary = {
        "spread_convention": "sample standard deviation (ddof=1), matching JPCB_tandem_round_2.tex",
        "optical_dielectric": 1.77,
        "manuscript_retained": stats([row["manuscript_J_cm"] for row in combined]),
        "regenerated_generic_CR2_fallback": stats([row["regenerated_J_cm"] for row in combined]),
        "regenerated_PDA": stats([row["regenerated_J_pda_cm"] for row in combined]),
        "regenerated_separation_A": stats([row["regenerated_separation_A"] for row in combined]),
        "method_note": (
            "The regenerated trajectory log explicitly reports the approximate generic CR2 fallback force field. "
            "It is suitable for the requested visualization but is not a like-for-like reproduction of the "
            "published RESP/AMBER-parameterized ensemble."
        ),
    }
    with open(args.out / "coupling_comparison_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = np.asarray([row["frame"] for row in combined])
    j_old = np.asarray([row["manuscript_J_cm"] for row in combined])
    j_new = np.asarray([row["regenerated_J_cm"] for row in combined])
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.8), constrained_layout=True)
    axes[0].plot(frame, j_old, lw=1.0, label="retained manuscript ensemble")
    axes[0].plot(frame, j_new, lw=1.0, label="regenerated (generic CR2 fallback)")
    axes[0].set(xlabel="NVT frame", ylabel=r"$J$ (cm$^{-1}$)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    bins = np.linspace(min(j_old.min(), j_new.min()), max(j_old.max(), j_new.max()), 42)
    axes[1].hist(j_old, bins=bins, alpha=0.68, label=f"retained: {j_old.mean():.1f} +/- {j_old.std(ddof=1):.1f}")
    axes[1].hist(j_new, bins=bins, alpha=0.58, label=f"regenerated: {j_new.mean():.1f} +/- {j_new.std(ddof=1):.1f}")
    axes[1].set(xlabel=r"$J$ (cm$^{-1}$)", ylabel="Frames")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    fig.savefig(args.out / "coupling_comparison.png", dpi=220)
    fig.savefig(args.out / "coupling_comparison.pdf")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

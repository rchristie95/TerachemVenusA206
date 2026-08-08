#!/usr/bin/env python3
"""Migrate archived TDC outputs after the inverse-Angstrom unit correction.

Old GPU results multiplied a Coulomb sum evaluated with distances in Angstrom
by ANGSTROM_TO_BOHR.  Reciprocal distances require BOHR_TO_ANGSTROM instead,
so every TDC ``J_cm`` value is multiplied by BOHR_TO_ANGSTROM**2.  PDA values
were already evaluated after converting coordinates to bohr and are unchanged.

The migration is idempotent: a JSON audit marker prevents double application.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import date
from pathlib import Path

import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903
FACTOR = BOHR_TO_ANGSTROM**2
SCALAR_KEYS = ("mean", "std", "min", "max", "median")


def correction_metadata() -> dict[str, object]:
    return {
        "status": "corrected",
        "date": date(2026, 8, 7).isoformat(),
        "factor_old_to_new": FACTOR,
        "old_reciprocal_distance_factor": 1.0 / BOHR_TO_ANGSTROM,
        "new_reciprocal_distance_factor": BOHR_TO_ANGSTROM,
        "reason": "Coulomb kernel distances are in Angstrom; 1/Angstrom to 1/bohr uses BOHR_TO_ANGSTROM",
        "pda_values_changed": False,
    }


def is_corrected(payload: dict[str, object]) -> bool:
    native = payload.get("tdc_units", {})
    if (
        isinstance(native, dict)
        and native.get("status") == "corrected"
        and native.get("pair_distance_unit") == "angstrom"
        and math.isclose(
            float(native.get("reciprocal_distance_to_atomic_units", float("nan"))),
            BOHR_TO_ANGSTROM,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    ):
        return True
    correction = payload.get("unit_correction", {})
    return bool(
        isinstance(correction, dict)
        and correction.get("status") == "corrected"
        and math.isclose(
            float(correction.get("new_reciprocal_distance_factor", float("nan"))),
            BOHR_TO_ANGSTROM,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    )


def migrate_directory(directory: Path) -> None:
    csv_path = directory / "coupling_samples.csv"
    json_path = directory / "coupling_distribution.json"
    if not csv_path.is_file() or not json_path.is_file():
        return

    summary = json.loads(json_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as handle:
        existing_rows = list(csv.DictReader(handle))

    rows = existing_rows
    already_corrected = is_corrected(summary)
    if not already_corrected:
        fieldnames = list(rows[0]) if rows else []
        if "J_cm" not in fieldnames:
            raise ValueError(f"{csv_path} has no J_cm column")
        for row in rows:
            row["J_cm"] = repr(float(row["J_cm"]) * FACTOR)
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        for key in SCALAR_KEYS:
            if key in summary:
                summary[key] = float(summary[key]) * FACTOR
        if "samples" in summary:
            summary["samples"] = [float(value) * FACTOR for value in summary["samples"]]
        summary["unit_correction"] = correction_metadata()
        json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"[corrected] {directory} by factor {FACTOR:.15f}")
    else:
        print(f"[skip] {directory} already corrected")

    # Geometry archives duplicate J_cm for downstream spectral generation.
    # Compare them with the final authoritative CSV (after any correction), so
    # a single migration call cannot leave an uncorrected NPZ beside a corrected
    # CSV. This synchronization is also safe on idempotent repeat calls.
    geometry_path = directory / "coupling_geometry.npz"
    if geometry_path.is_file() and rows:
        csv_j = np.asarray([float(row["J_cm"]) for row in rows])
        with np.load(geometry_path) as archive:
            arrays = {key: np.asarray(archive[key]) for key in archive.files}
        if "J_cm" in arrays:
            if np.allclose(arrays["J_cm"], csv_j, rtol=0.0, atol=1.0e-10):
                return
            if np.allclose(arrays["J_cm"] * FACTOR, csv_j, rtol=0.0, atol=1.0e-10):
                arrays["J_cm"] = arrays["J_cm"] * FACTOR
                np.savez(geometry_path, **arrays)
                print(f"[synchronised] {geometry_path}")
                return
            raise ValueError(f"{geometry_path}: J_cm matches neither raw nor corrected CSV")


def _rewrite_csv_column(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    fieldnames = list(rows[0]) if rows else []
    missing = [column for column in columns if column not in fieldnames]
    if missing:
        raise ValueError(f"{path}: missing expected TDC columns {missing}")
    for row in rows:
        for column in columns:
            row[column] = repr(float(row[column]) * FACTOR)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _rewrite_npz_j(path: Path) -> None:
    with np.load(path) as archive:
        arrays = {key: np.asarray(archive[key]) for key in archive.files}
    if "J_cm" not in arrays:
        raise ValueError(f"{path}: no J_cm array")
    arrays["J_cm"] = arrays["J_cm"] * FACTOR
    np.savez(path, **arrays)


def _migrate_prefixed_distribution(directory: Path, suffix: str) -> None:
    json_path = directory / f"coupling_distribution_{suffix}.json"
    csv_path = directory / f"coupling_samples_{suffix}.csv"
    npz_path = directory / f"coupling_geometry_{suffix}.npz"
    if not (json_path.is_file() and csv_path.is_file() and npz_path.is_file()):
        return
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    if is_corrected(summary):
        print(f"[skip] {json_path} already corrected")
        return
    _rewrite_csv_column(csv_path, ("J_cm",))
    _rewrite_npz_j(npz_path)
    for key in SCALAR_KEYS:
        if key in summary:
            summary[key] = float(summary[key]) * FACTOR
    if "samples" in summary:
        summary["samples"] = [float(value) * FACTOR for value in summary["samples"]]
    summary["unit_correction"] = correction_metadata()
    summary["archive_status"] = "superseded generic-CR2 diagnostic; not manuscript evidence"
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[corrected] {json_path}")


def _migrate_prefixed_comparison(directory: Path, suffix: str) -> None:
    stem = "coupling_comparison_summary" + (f"_{suffix}" if suffix else "")
    csv_stem = "coupling_frame_comparison" + (f"_{suffix}" if suffix else "")
    json_path = directory / f"{stem}.json"
    csv_path = directory / f"{csv_stem}.csv"
    if not (json_path.is_file() and csv_path.is_file()):
        return
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    if is_corrected(summary):
        print(f"[skip] {json_path} already corrected")
        return
    _rewrite_csv_column(csv_path, ("manuscript_J_cm", "regenerated_J_cm"))
    for block_name in ("manuscript_retained", "regenerated_generic_CR2_fallback"):
        block = summary[block_name]
        for key in ("mean", "std_sample", "min", "max", "median"):
            block[key] = float(block[key]) * FACTOR
    summary["unit_correction"] = correction_metadata()
    summary["archive_status"] = "superseded generic-CR2 diagnostic; not manuscript evidence"
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[corrected] {json_path}")


def _regenerate_prefixed_plots(directory: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for suffix in ("", "full_gpu"):
        suffix_part = f"_{suffix}" if suffix else ""
        path = directory / f"coupling_frame_comparison{suffix_part}.csv"
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        frame = np.asarray([int(row["frame"]) for row in rows])
        retained = np.asarray([float(row["manuscript_J_cm"]) for row in rows])
        regenerated = np.asarray([float(row["regenerated_J_cm"]) for row in rows])
        fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.8), constrained_layout=True)
        fig.suptitle("Superseded generic-CR2 diagnostic (unit-corrected)")
        axes[0].plot(frame, retained, lw=1.0, label="legacy retained ensemble")
        axes[0].plot(frame, regenerated, lw=1.0, label="generic CR2 fallback")
        axes[0].set(xlabel="NVT frame", ylabel=r"$J$ (cm$^{-1}$)")
        axes[0].grid(alpha=0.25)
        axes[0].legend(frameon=False)
        bins = np.linspace(min(retained.min(), regenerated.min()), max(retained.max(), regenerated.max()), 42)
        axes[1].hist(retained, bins=bins, alpha=0.68, label=f"legacy: {retained.mean():.1f} +/- {retained.std(ddof=1):.1f}")
        axes[1].hist(regenerated, bins=bins, alpha=0.58, label=f"fallback: {regenerated.mean():.1f} +/- {regenerated.std(ddof=1):.1f}")
        axes[1].set(xlabel=r"$J$ (cm$^{-1}$)", ylabel="Frames")
        axes[1].grid(alpha=0.25)
        axes[1].legend(frameon=False)
        fig.savefig(directory / f"coupling_comparison{suffix_part}.png", dpi=220)
        fig.savefig(directory / f"coupling_comparison{suffix_part}.pdf")
        plt.close(fig)

    sample_path = directory / "coupling_samples_full_gpu.csv"
    if sample_path.is_file():
        with sample_path.open(newline="", encoding="utf-8") as handle:
            values = np.asarray([float(row["J_cm"]) for row in csv.DictReader(handle)])
        fig, ax = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)
        ax.hist(values, bins=42, color="#3b7ddd", alpha=0.82)
        ax.axvline(values.mean(), color="black", ls="--", lw=1.2, label=f"{values.mean():.2f} +/- {values.std(ddof=1):.2f}")
        ax.set(title="Superseded generic-CR2 diagnostic (unit-corrected)", xlabel=r"$J$ (cm$^{-1}$)", ylabel="Frames")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
        fig.savefig(directory / "Fig_Coupling_Histogram_full_gpu.png", dpi=220)
        fig.savefig(directory / "Fig_Coupling_Histogram_full_gpu.pdf")
        plt.close(fig)


def migrate_prefixed_legacy(directory: Path) -> None:
    """Migrate the explicitly superseded, non-production 2026-07-21 archive."""
    if not any(directory.glob("coupling_distribution_*.json")):
        return
    for suffix in ("regenerated", "full_gpu"):
        _migrate_prefixed_distribution(directory, suffix)
    for suffix in ("", "full_gpu"):
        _migrate_prefixed_comparison(directory, suffix)
    _regenerate_prefixed_plots(directory)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    for root in args.directories:
        root = root.resolve()
        migrate_directory(root)
        migrate_prefixed_legacy(root)
        for json_path in sorted(root.rglob("coupling_distribution.json")):
            migrate_directory(json_path.parent)


if __name__ == "__main__":
    main()

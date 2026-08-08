#!/usr/bin/env python3
"""Cross-check every retained coupling archive, including locally ignored ones."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RECIPROCAL_FACTOR = 0.529177210903


def is_corrected(payload: dict) -> bool:
    native = payload.get("tdc_units", {})
    migrated = payload.get("unit_correction", {})
    return bool(
        (
            native.get("status") == "corrected"
            and native.get("pair_distance_unit") == "angstrom"
            and abs(
                float(native.get("reciprocal_distance_to_atomic_units", np.nan))
                - RECIPROCAL_FACTOR
            )
            < 1.0e-12
        )
        or (
            migrated.get("status") == "corrected"
            and abs(
                float(migrated.get("new_reciprocal_distance_factor", np.nan))
                - RECIPROCAL_FACTOR
            )
            < 1.0e-12
        )
    )


def distribution_paths() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("coupling_distribution*.json")
        if "tmp" not in path.relative_to(ROOT).parts
    )


def test_all_coupling_archives_are_unit_marked_and_internally_consistent() -> None:
    paths = distribution_paths()
    assert paths
    for json_path in paths:
        payload = json.loads(json_path.read_text())
        assert is_corrected(payload), json_path

        suffix = json_path.stem.removeprefix("coupling_distribution")
        csv_path = json_path.with_name(f"coupling_samples{suffix}.csv")
        assert csv_path.is_file(), csv_path
        with csv_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        values = np.asarray([float(row["J_cm"]) for row in rows], dtype=float)
        assert len(values) == int(payload["n"])
        assert np.array_equal(values, np.asarray(payload["samples"], dtype=float))
        assert abs(float(payload["mean"]) - float(values.mean())) < 1.0e-12
        expected_std = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        assert abs(float(payload["std"]) - expected_std) < 1.0e-12
        assert abs(float(payload["min"]) - float(values.min())) < 1.0e-12
        assert abs(float(payload["max"]) - float(values.max())) < 1.0e-12
        assert abs(float(payload["median"]) - float(np.median(values))) < 1.0e-12

        geometry_path = json_path.with_name(f"coupling_geometry{suffix}.npz")
        if geometry_path.is_file():
            with np.load(geometry_path, allow_pickle=False) as archive:
                assert np.allclose(
                    np.asarray(archive["J_cm"], dtype=float),
                    values,
                    rtol=0.0,
                    atol=1.0e-10,
                ), geometry_path


def test_all_prefixed_comparison_archives_are_consistent() -> None:
    paths = sorted(
        path
        for path in ROOT.rglob("coupling_comparison_summary*.json")
        if "tmp" not in path.relative_to(ROOT).parts
    )
    for json_path in paths:
        payload = json.loads(json_path.read_text())
        assert is_corrected(payload), json_path
        suffix = json_path.stem.removeprefix("coupling_comparison_summary")
        csv_path = json_path.with_name(f"coupling_frame_comparison{suffix}.csv")
        assert csv_path.is_file(), csv_path
        with csv_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        for block_name, column in (
            ("manuscript_retained", "manuscript_J_cm"),
            ("regenerated_generic_CR2_fallback", "regenerated_J_cm"),
        ):
            values = np.asarray([float(row[column]) for row in rows], dtype=float)
            block = payload[block_name]
            assert int(block["n"]) == len(values)
            assert abs(float(block["mean"]) - float(values.mean())) < 1.0e-12
            assert abs(float(block["std_sample"]) - float(values.std(ddof=1))) < 1.0e-12
            assert abs(float(block["min"]) - float(values.min())) < 1.0e-12
            assert abs(float(block["max"]) - float(values.max())) < 1.0e-12
            assert abs(float(block["median"]) - float(np.median(values))) < 1.0e-12


def test_decoherence_note_uses_the_corrected_production_coupling() -> None:
    audit = json.loads((ROOT / "notes/decoherence_note_figure_audit.json").read_text())
    with (
        ROOT / "coupling_nvt_production_cr2_1000_20260721/coupling_samples.csv"
    ).open(newline="") as handle:
        values = np.asarray(
            [float(row["J_cm"]) for row in csv.DictReader(handle)], dtype=float
        )

    mean_j = float(values.mean())
    assert audit["coupling_unit_status"] == "corrected reciprocal-distance conversion"
    assert abs(float(audit["mean_J_cm-1"]) - mean_j) < 1.0e-12
    assert abs(float(audit["mean_splitting_cm-1"]) - 2.0 * abs(mean_j)) < 1.0e-12
    assert abs(float(audit["peak_ratio_24_over_60"]) - 0.1733481453374537) < 1.0e-12

    prose = "\n".join(
        (ROOT / relative).read_text()
        for relative in (
            "notes/solvation_decoherence_note.tex",
            "solvation_decoherence_test/HANDOFF.md",
            "solvation_decoherence_test/NUMERICAL_TEST_REPORT.md",
        )
    )
    for stale in ("117.19", "234.38", "520.69", "527.12", "45.3"):
        assert stale not in prose
    assert "32.82" in prose and "65.63" in prose and "161.7" in prose
    assert (ROOT / "notes/Fig_T2_CommonScale.pdf").is_file()

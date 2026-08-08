"""Exercise the one-off migration used for superseded prefixed archives."""

from __future__ import annotations

import csv
import json

import numpy as np

from correct_coupling_unit_outputs import FACTOR, migrate_directory, migrate_prefixed_legacy


def test_prefixed_legacy_migration_is_exact_and_idempotent(tmp_path):
    distribution = {
        "n": 2,
        "mean": 15.0,
        "std": np.sqrt(50.0),
        "min": 10.0,
        "max": 20.0,
        "median": 15.0,
        "samples": [10.0, 20.0],
    }
    (tmp_path / "coupling_distribution_full_gpu.json").write_text(
        json.dumps(distribution), encoding="utf-8"
    )
    sample_fields = ["frame", "J_cm", "J_pda_cm"]
    with (tmp_path / "coupling_samples_full_gpu.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields)
        writer.writeheader()
        writer.writerows(
            [
                {"frame": 0, "J_cm": 10.0, "J_pda_cm": 3.0},
                {"frame": 1, "J_cm": 20.0, "J_pda_cm": 4.0},
            ]
        )
    np.savez(
        tmp_path / "coupling_geometry_full_gpu.npz",
        frame=np.array([0, 1]),
        J_cm=np.array([10.0, 20.0]),
        epsilon=np.array(1.77),
    )

    comparison = {
        "manuscript_retained": {
            "n": 2,
            "mean": 15.0,
            "std_sample": np.sqrt(50.0),
            "min": 10.0,
            "max": 20.0,
            "median": 15.0,
        },
        "regenerated_generic_CR2_fallback": {
            "n": 2,
            "mean": 35.0,
            "std_sample": np.sqrt(50.0),
            "min": 30.0,
            "max": 40.0,
            "median": 35.0,
        },
        "regenerated_PDA": {"mean": 3.5},
    }
    (tmp_path / "coupling_comparison_summary_full_gpu.json").write_text(
        json.dumps(comparison), encoding="utf-8"
    )
    comparison_fields = [
        "frame",
        "manuscript_J_cm",
        "manuscript_J_pda_cm",
        "regenerated_J_cm",
        "regenerated_J_pda_cm",
    ]
    with (tmp_path / "coupling_frame_comparison_full_gpu.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_fields)
        writer.writeheader()
        writer.writerows(
            [
                {
                    "frame": 0,
                    "manuscript_J_cm": 10.0,
                    "manuscript_J_pda_cm": 3.0,
                    "regenerated_J_cm": 30.0,
                    "regenerated_J_pda_cm": 3.0,
                },
                {
                    "frame": 1,
                    "manuscript_J_cm": 20.0,
                    "manuscript_J_pda_cm": 4.0,
                    "regenerated_J_cm": 40.0,
                    "regenerated_J_pda_cm": 4.0,
                },
            ]
        )

    migrate_prefixed_legacy(tmp_path)
    tracked = [
        tmp_path / "coupling_distribution_full_gpu.json",
        tmp_path / "coupling_samples_full_gpu.csv",
        tmp_path / "coupling_geometry_full_gpu.npz",
        tmp_path / "coupling_comparison_summary_full_gpu.json",
        tmp_path / "coupling_frame_comparison_full_gpu.csv",
    ]
    first_bytes = {path: path.read_bytes() for path in tracked}

    with (tmp_path / "coupling_samples_full_gpu.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert np.allclose([float(row["J_cm"]) for row in rows], [10 * FACTOR, 20 * FACTOR])
    assert [float(row["J_pda_cm"]) for row in rows] == [3.0, 4.0]
    with np.load(tmp_path / "coupling_geometry_full_gpu.npz") as archive:
        assert np.allclose(archive["J_cm"], [10 * FACTOR, 20 * FACTOR])
        assert float(archive["epsilon"]) == 1.77

    migrate_prefixed_legacy(tmp_path)
    assert {path: path.read_bytes() for path in tracked} == first_bytes


def test_native_corrected_distribution_is_never_rescaled(tmp_path):
    summary = {
        "n": 1,
        "mean": 30.453973388943076,
        "samples": [30.453973388943076],
        "tdc_units": {
            "status": "corrected",
            "pair_distance_unit": "angstrom",
            "reciprocal_distance_to_atomic_units": 0.529177210903,
        },
    }
    csv_path = tmp_path / "coupling_samples.csv"
    json_path = tmp_path / "coupling_distribution.json"
    csv_path.write_text("frame,J_cm,J_pda_cm\n0,30.453973388943076,24.814298716448004\n")
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    before = {path: path.read_bytes() for path in (csv_path, json_path)}

    migrate_directory(tmp_path)

    assert {path: path.read_bytes() for path in (csv_path, json_path)} == before


def test_standard_migration_synchronises_geometry_on_first_call(tmp_path):
    summary = {
        "n": 2,
        "mean": 15.0,
        "std": np.sqrt(50.0),
        "min": 10.0,
        "max": 20.0,
        "median": 15.0,
        "samples": [10.0, 20.0],
    }
    (tmp_path / "coupling_distribution.json").write_text(
        json.dumps(summary), encoding="utf-8"
    )
    (tmp_path / "coupling_samples.csv").write_text(
        "frame,J_cm,J_pda_cm\n0,10.0,3.0\n1,20.0,4.0\n", encoding="utf-8"
    )
    np.savez(
        tmp_path / "coupling_geometry.npz",
        frame=np.array([0, 1]),
        J_cm=np.array([10.0, 20.0]),
        J_pda_cm=np.array([3.0, 4.0]),
    )

    migrate_directory(tmp_path)

    with np.load(tmp_path / "coupling_geometry.npz") as archive:
        assert np.allclose(archive["J_cm"], [10 * FACTOR, 20 * FACTOR])
        assert np.array_equal(archive["J_pda_cm"], [3.0, 4.0])
    first_bytes = {
        path: path.read_bytes()
        for path in (
            tmp_path / "coupling_distribution.json",
            tmp_path / "coupling_samples.csv",
            tmp_path / "coupling_geometry.npz",
        )
    }

    migrate_directory(tmp_path)

    assert {
        path: path.read_bytes()
        for path in (
            tmp_path / "coupling_distribution.json",
            tmp_path / "coupling_samples.csv",
            tmp_path / "coupling_geometry.npz",
        )
    } == first_bytes

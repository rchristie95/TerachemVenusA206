#!/usr/bin/env python3
"""Verify that archived coupling outputs received only the intended unit fix."""

from __future__ import annotations

import csv
import io
import json
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FACTOR = 0.529177210903**2
# Immutable repository state immediately before the reciprocal-distance fix.
# Pinning this commit keeps the migration audit valid after the corrected files
# themselves are committed; using HEAD would then compare the files to themselves.
PRE_CORRECTION_COMMIT = "34775aa82847f654d77f8720182ddc36782095c7"
DIRECTORIES = (
    "coupling_paper_steom_thermal",
    "coupling_tandem_1000",
    "coupling_nvt_production_cr2_1000_20260721",
    "coupling_nvt_production_cr2_1000_20260721/precision_check_f64/frame_0000",
    "coupling_nvt_production_cr2_1000_20260721/precision_check_f64/frame_0499",
    "coupling_nvt_production_cr2_1000_20260721/precision_check_f64/frame_0999",
)


def head_bytes(relative: str) -> bytes:
    return subprocess.check_output(
        ["git", "show", f"{PRE_CORRECTION_COMMIT}:{relative}"],
        cwd=ROOT,
        stderr=subprocess.DEVNULL,
    )


def csv_rows(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode("utf-8"))))


def verify_csv(relative: str) -> tuple[float, int]:
    before = csv_rows(head_bytes(relative))
    after = csv_rows((ROOT / relative).read_bytes())
    assert len(before) == len(after)
    old_j = np.asarray([float(row["J_cm"]) for row in before])
    new_j = np.asarray([float(row["J_cm"]) for row in after])
    ratio_error = float(np.max(np.abs(new_j / old_j - FACTOR)))
    assert ratio_error < 5.0e-16, (relative, ratio_error)
    for key in before[0]:
        if key == "J_cm":
            continue
        assert [row[key] for row in after] == [row[key] for row in before], (
            relative,
            key,
        )
    return ratio_error, len(before)


def verify_npz(relative: str) -> None:
    with np.load(io.BytesIO(head_bytes(relative))) as archive:
        before = {key: np.asarray(archive[key]) for key in archive.files}
    with np.load(ROOT / relative) as archive:
        after = {key: np.asarray(archive[key]) for key in archive.files}
    assert before.keys() == after.keys()
    for key in before:
        if key == "J_cm":
            assert np.allclose(after[key], before[key] * FACTOR, rtol=0.0, atol=1.0e-12)
        else:
            assert np.array_equal(after[key], before[key]), (relative, key)


def main() -> None:
    maximum_error = 0.0
    rows = 0
    for directory in DIRECTORIES:
        error, count = verify_csv(f"{directory}/coupling_samples.csv")
        maximum_error = max(maximum_error, error)
        rows += count
        npz_relative = f"{directory}/coupling_geometry.npz"
        if (ROOT / npz_relative).is_file():
            verify_npz(npz_relative)
    print(
        f"PASS: {rows} archived rows; TDC factor={FACTOR:.15f}; "
        f"maximum ratio error={maximum_error:.3e}; all non-TDC fields unchanged"
    )


def test_migrated_archives() -> None:
    main()


def test_static_recomputed_with_production_density() -> None:
    """The static control was recomputed, not blindly rescaled.

    It must use the same definitive seven-pair, cap-masked density as the
    production ensemble and the independently regenerated multipole table.
    """
    relative = Path("coupling_paper_steom_static/coupling_samples.csv")
    rows = csv_rows((ROOT / relative).read_bytes())
    assert len(rows) == 1
    row = rows[0]
    assert float(row["J_cm"]) == 30.453973388943076
    assert float(row["J_pda_cm"]) == 24.814298716448004
    assert float(row["separation_A"]) == 25.209016247697047
    assert float(row["angle_deg"]) == 103.68594446779686

    distribution = json.loads(
        (ROOT / "coupling_paper_steom_static/coupling_distribution.json").read_text()
    )
    provenance = distribution["density_provenance"]
    assert distribution["mean"] == float(row["J_cm"])
    assert distribution["tdc_units"]["reciprocal_distance_to_atomic_units"] == 0.529177210903
    assert provenance["source"].endswith("steom_transdens_capmasked_oldframe.npz")
    assert provenance["retained_points"] == 259277
    assert provenance["included_nto_pairs"] == 7
    assert abs(provenance["represented_nto_occupation"] - 0.99045993) < 1.0e-12


if __name__ == "__main__":
    main()

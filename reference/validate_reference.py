#!/usr/bin/env python3
"""Validate compact published data without ORCA, OpenMM, PyMOL, or a GPU."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def validate_checksums() -> None:
    for raw in (ROOT / "reference/SHA256SUMS").read_text().splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        expected, rel = raw.split(maxsplit=1)
        path = ROOT / rel
        if not path.exists():
            print(f"SKIP external/generated: {rel}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"checksum mismatch: {rel}"


def validate_tandem_statistics() -> None:
    ref = json.loads((ROOT / "reference/orca_validation.json").read_text())
    with (ROOT / "coupling_tandem_1000/coupling_samples.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = [float(row["J_cm"]) for row in rows]
    expected = ref["tandem_ensemble"]
    assert len(values) == expected["frames"] == 1000
    assert abs(statistics.fmean(values) - expected["mean_J_cm-1"]) < 1e-12
    # coupling_ensemble.py reports the sample standard deviation (ddof=1).
    assert abs(statistics.stdev(values) - expected["std_J_cm-1"]) < 1e-12
    assert min(values) == expected["min_J_cm-1"]
    assert max(values) == expected["max_J_cm-1"]
    assert [int(row["frame"]) for row in rows] == list(range(1000))


def validate_qm_inputs() -> None:
    geom = (ROOT / "neo_model/orca_steom/geom_cthrp.xyz").read_text().splitlines()
    field = (ROOT / "neo_model/orca_steom/field.pc").read_text().splitlines()
    assert int(geom[0]) == 44
    assert int(field[0]) == 2350
    assert len(geom[2:]) == 44
    assert len(field[1:]) == 2350


if __name__ == "__main__":
    validate_checksums()
    validate_tandem_statistics()
    validate_qm_inputs()
    print("reference validation: PASS")

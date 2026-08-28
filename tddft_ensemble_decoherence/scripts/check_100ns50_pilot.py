#!/usr/bin/env python3
"""Require a complete corrected-cap pilot before the remaining 49 frames."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    calibration = json.loads(
        (ROOT / "rapid_100ns50/calibration.json").read_text()
    )
    expected_delta = float(calibration["steom_constant_shift_eV"])
    manifest = ROOT / "manifests/long100ns50_frames.csv"
    with manifest.open(newline="") as handle:
        row = next(csv.DictReader(handle))
    if row["convergence_status"] != "complete":
        raise RuntimeError(
            f"Pilot frame 0 is not complete: {row['convergence_status']} "
            f"{row['failure_class']} {row['notes']}"
        )
    prep = json.loads(
        (ROOT / "frames/long100ns50/frame_0000/preparation.json").read_text()
    )
    caps = [
        cap
        for site in ("A", "B")
        for cap in prep["sites"][site]["link_caps"]
    ]
    if len(caps) != 6:
        raise RuntimeError(f"Expected six site/link diagnostics, found {len(caps)}")
    if any(abs(cap["cap_direction_dot_parent_to_partner"] - 1.0) > 1.0e-12 for cap in caps):
        raise RuntimeError("A pilot link cap is not directed toward its MM partner")
    delta_a = float(row["site_A_corrected_eV"]) - float(row["site_A_energy_eV"])
    delta_b = float(row["site_B_corrected_eV"]) - float(row["site_B_energy_eV"])
    if abs(delta_a - expected_delta) > 1.0e-8 or abs(delta_b - expected_delta) > 1.0e-8:
        raise RuntimeError("Pilot STEOM reference shift mismatch")
    if int(row["site_A_selected_root"]) == 4 or int(row["site_B_selected_root"]) == 4:
        raise RuntimeError("Pilot selected a boundary root")
    print(
        "[pilot-valid] frame=0 time=2.0 ns caps=6 dot=1.0 "
        f"roots={row['site_A_selected_root']}/{row['site_B_selected_root']} "
        f"Delta={delta_a:+.8f} eV"
    )


if __name__ == "__main__":
    main()

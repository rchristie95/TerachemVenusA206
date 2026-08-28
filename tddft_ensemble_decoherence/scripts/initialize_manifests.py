#!/usr/bin/env python3
"""Initialize missing-input manifests without overwriting classified rows."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

COLUMNS = [
    "window", "frame_index", "segment_frame_index", "time_fs", "time_ps",
    "source_trajectory", "coordinate_sha256", "site_A_energy_eV",
    "site_A_energy_cm-1", "site_A_wavelength_nm", "site_B_energy_eV",
    "site_B_energy_cm-1", "site_B_wavelength_nm", "site_A_corrected_eV",
    "site_B_corrected_eV", "site_A_oscillator_strength",
    "site_B_oscillator_strength", "site_A_selected_root",
    "site_B_selected_root", "site_A_assignment_score",
    "site_B_assignment_score", "J_vacuum_cm-1", "J_screened_cm-1",
    "J_pda_vacuum_cm-1", "J_pda_screened_cm-1", "separation_A",
    "dipole_angle_deg", "grid_normalization_error_A",
    "grid_normalization_error_B", "convergence_status", "failure_class",
    "retry_count", "runtime_seconds", "gpu", "compact_output_path", "notes",
]


def existing_rows(path: Path) -> dict[int, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {int(row["frame_index"]): row for row in csv.DictReader(handle)}


def write_manifest(path: Path, planned: list[dict[str, object]]) -> None:
    retained = existing_rows(path)
    rows = []
    for row in planned:
        old = retained.get(int(row["frame_index"]))
        if old:
            rows.append(old)
        else:
            rows.append({column: row.get(column, "") for column in COLUMNS})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def short_rows() -> list[dict[str, object]]:
    rows = []
    for frame in range(2000):
        segment = 0 if frame < 1000 else 1
        local = frame if segment == 0 else frame - 1000
        time_fs = float(frame * 4)
        rows.append({
            "window": "short",
            "frame_index": frame,
            "segment_frame_index": local,
            "time_fs": f"{time_fs:.1f}",
            "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": (
                "external/full_highres_0to4ps_4fs.dcd" if segment == 0
                else "external/full_highres_4to8ps_4fs.dcd"
            ),
            "convergence_status": "blocked_missing_source",
            "failure_class": "missing_source_trajectory",
            "retry_count": 0,
            "notes": "No coordinate substitution or interpolation permitted",
        })
    return rows


def long_rows() -> list[dict[str, object]]:
    rows = []
    for frame in range(1000):
        time_ps = float(frame + 1)
        rows.append({
            "window": "long",
            "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_ps * 1000.0:.1f}",
            "time_ps": f"{time_ps:.6f}",
            "source_trajectory": "external/tandem_nvt_production_cr2_1000.dcd",
            "convergence_status": "blocked_missing_source",
            "failure_class": "missing_source_trajectory",
            "retry_count": 0,
            "notes": "Compact coupling exists, but full-system coordinates/MM embedding do not",
        })
    return rows


def long100ns50_rows() -> list[dict[str, object]]:
    rows = []
    for frame in range(50):
        time_ps = float((frame + 1) * 2000)
        rows.append({
            "window": "long100ns50",
            "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_ps * 1000.0:.1f}",
            "time_ps": f"{time_ps:.6f}",
            "source_trajectory": (
                "tc_tandem_nvt_100ns_50frames_20260723/"
                "tandem_nvt_100ns_50frames.dcd"
            ),
            "convergence_status": "blocked_missing_source",
            "failure_class": "missing_source_trajectory",
            "retry_count": 0,
            "notes": (
                "Planned intact-MM snapshots at 2 ns cadence over 2--100 ns; "
                "no interpolation"
            ),
        })
    return rows


def short100fs50_rows(window: str = "short100fs50") -> list[dict[str, object]]:
    rows = []
    for frame in range(50):
        time_fs = float((frame + 1) * 2)
        rows.append({
            "window": window,
            "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_fs:.1f}",
            "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": (
                "tc_tandem_nvt_100fs_50frames_20260723/"
                "tandem_nvt_100fs_50frames.dcd"
            ),
            "convergence_status": "blocked_missing_source",
            "failure_class": "missing_source_trajectory",
            "retry_count": 0,
            "notes": (
                "Fifty consecutive intact-MM frames at 2 fs cadence over 2--100 fs; "
                "no interpolation"
            ),
        })
    return rows


def short100fs100_rows() -> list[dict[str, object]]:
    rows = []
    for frame in range(100):
        time_fs = float(frame + 1)
        rows.append({
            "window": "short100fs100_wb97xd3",
            "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_fs:.1f}",
            "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": (
                "tc_tandem_nvt_100fs_100frames_20260724/"
                "tandem_nvt_100fs_100frames.dcd"
            ),
            "convergence_status": "blocked_missing_source",
            "failure_class": "missing_source_trajectory",
            "retry_count": 0,
            "notes": (
                "One hundred consecutive intact-MM frames at 1 fs cadence over "
                "1--100 fs; wB97X-D3/6-311G** paired-site campaign; no interpolation"
            ),
        })
    return rows


def main() -> None:
    manifests = ROOT / "manifests"
    manifests.mkdir(parents=True, exist_ok=True)
    write_manifest(manifests / "short_frames.csv", short_rows())
    write_manifest(manifests / "long_frames.csv", long_rows())
    write_manifest(manifests / "long100ns50_frames.csv", long100ns50_rows())
    write_manifest(manifests / "short100fs50_frames.csv", short100fs50_rows())
    write_manifest(
        manifests / "short100fs50_fixedmm_frames.csv",
        short100fs50_rows("short100fs50_fixedmm"),
    )
    write_manifest(
        manifests / "short100fs100_wb97xd3_frames.csv",
        short100fs100_rows(),
    )
    print(
        "Initialized 2000 short, 1000 long, 50 long100ns50, 50 short100fs50, "
        "50 short100fs50_fixedmm, and 100 short100fs100_wb97xd3 "
        "restart rows; "
        "existing rows preserved."
    )


if __name__ == "__main__":
    main()

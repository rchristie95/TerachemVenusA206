#!/usr/bin/env python3
"""Preserve seven-root boundary attempts and retry them with nine CAM-B3LYP roots."""

from __future__ import annotations

import fcntl
import json
from pathlib import Path

import run_rapid_adaptive_100ns50 as adaptive
from run_ensemble import ROOT, load_manifest, save_manifest


WINDOW = "short100fs50"
MANIFEST = ROOT / "manifests/short100fs50_frames.csv"
OUTDIR = ROOT / "rapid_100fs50"


def archive_canonical(frame: int) -> Path:
    canonical = ROOT / f"frames/{WINDOW}/frame_{frame:04d}"
    if not canonical.exists():
        raise RuntimeError(f"Missing retained boundary attempt: {canonical}")
    attempt = 0
    while True:
        target = ROOT / (
            f"frames/{WINDOW}/frame_{frame:04d}_attempt{attempt}_rapid_cam_7roots"
        )
        if not target.exists():
            canonical.rename(target)
            return target
        attempt += 1


def main() -> None:
    lock_handle = (OUTDIR / "campaign.lock").open("w")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        raise RuntimeError("Another campaign driver still holds the lock") from error
    calibration = json.loads((OUTDIR / "calibration.json").read_text())
    delta = float(calibration["steom_constant_shift_eV"])
    reference = Path(calibration["state_tracking_reference"])
    columns, rows = load_manifest(MANIFEST)
    targets = [
        index
        for index, row in enumerate(rows)
        if row["convergence_status"] != "complete"
        and (
            row["failure_class"] in {"", "missing_target_state_or_root_flipping"}
        )
    ]
    if not targets:
        print("[nine-root-retry] no boundary frames remain", flush=True)
        return
    archives = {}
    for frame in targets:
        archives[str(frame)] = str(archive_canonical(frame).resolve())
        rows[frame]["convergence_status"] = "pending_tddft"
        rows[frame]["failure_class"] = ""
        rows[frame]["retry_count"] = str(
            int(rows[frame].get("retry_count") or 0) + 1
        )
        rows[frame]["notes"] = (
            "Seven-root physical state landed at the computed-root boundary; "
            "original retained and queued for nine-root CAM-B3LYP retry"
        )
    save_manifest(MANIFEST, columns, rows)

    adaptive.WINDOW = WINDOW
    adaptive.MANIFEST = MANIFEST
    adaptive.OUTDIR = OUTDIR
    print(
        f"[nine-root-retry] frames={targets} workers=3 archives={archives}",
        flush=True,
    )
    result = adaptive.run_batch(
        targets,
        3,
        rows,
        columns,
        "0",
        9,
        "rapid_cam",
        delta,
        reference,
    )
    final_failed = [
        int(row["frame_index"])
        for row in rows
        if row["convergence_status"] == "failed"
    ]
    audit = {
        "status": "complete" if not final_failed else "complete_with_failures",
        "reason": "selected state at seven-root boundary",
        "roots": 9,
        "profile": "rapid_cam",
        "workers": 3,
        "target_frames": targets,
        "archives": archives,
        "batch": result,
        "final_failed_frames": final_failed,
    }
    (OUTDIR / "nine_root_retry_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print(
        f"[nine-root-retry-complete] failed={final_failed}",
        flush=True,
    )


if __name__ == "__main__":
    main()

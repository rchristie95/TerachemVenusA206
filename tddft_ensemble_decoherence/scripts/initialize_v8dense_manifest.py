#!/usr/bin/env python3
"""Manifest for the 1 ps / 2 fs dense segment (v8dense).

500 frames at 2 fs. Unlike the TDX run, --npt-steps is 0 and the segment is
seeded directly from a v5 frame, so the DCD holds only production frames and
frame_index runs 0-499 with no equilibration offset. That is verified below
against the frame count rather than assumed.

time_fs = (frame + 1) * 2, because OpenMM writes after each interval and never
at step 0.
"""
from __future__ import annotations
import csv, hashlib
from pathlib import Path
import numpy as np
from initialize_manifests import write_manifest
from register_nvt_trajectories import dcd_frames

ROOT = Path(__file__).resolve().parents[1]
WINDOW, FRAMES, CADENCE_FS = "v8dense", 500, 2.0
SEED_FRAME = 1756
TRAJ = ROOT.parent / "tc_tandem_v8_dense" / "tandem_v8_1ps_2fs.dcd"


def main() -> None:
    if not TRAJ.exists():
        raise FileNotFoundError(TRAJ)
    rows = [{
        "window": WINDOW, "frame_index": i, "segment_frame_index": i,
        "time_fs": f"{(i + 1) * CADENCE_FS:.1f}",
        "time_ps": f"{(i + 1) * CADENCE_FS / 1000.0:.6f}",
        "source_trajectory": str(TRAJ),
        "convergence_status": "pending_tddft", "failure_class": "", "retry_count": 0,
        "notes": (f"dense segment for the cumulant lineshape; 1 ps at 2 fs, seeded "
                  f"from v5camp frame {SEED_FRAME} (median |Delta|); svp profile; "
                  f"full-system embedding, link-only exclusion"),
    } for i in range(FRAMES)]
    path = ROOT / "manifests" / f"{WINDOW}_frames.csv"
    write_manifest(path, rows)

    with path.open(newline="") as h:
        r = csv.DictReader(h); cols = r.fieldnames; out = list(r)
    n = 0
    for i, xyz in enumerate(dcd_frames(TRAJ)):
        if i >= len(out):
            raise RuntimeError(f"{TRAJ.name} has more frames than the manifest")
        if out[i]["convergence_status"] == "complete":
            continue
        out[i]["coordinate_sha256"] = hashlib.sha256(
            np.asarray(xyz, dtype="<f4").tobytes(order="C")).hexdigest()
        n += 1
    if n != len(out):
        raise RuntimeError(f"stamped {n} of {len(out)}; trajectory is short")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cols); w.writeheader(); w.writerows(out)
    tmp.replace(path)
    print(f"wrote {path.name}: {FRAMES} rows at {CADENCE_FS} fs, checksums bound")


if __name__ == "__main__":
    main()

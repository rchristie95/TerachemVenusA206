#!/usr/bin/env python3
"""Manifests for the 500 fs low-frequency windows.

Only frames BEYOND 100 fs are listed: DCD frames 25-124, i.e. t = 104-500 fs.
The first 100 fs is already covered at 1 fs by v6t2a/b/c, and recomputing it at
4 fs would buy nothing.

Rows are 0-99; the DCD frame to read is each row's own frame_index. run_ensemble
handles that split (`dcd_frame = int(row.get("frame_index") or frame)`), which is
the same pattern v5camp used for its strided production frames.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

from initialize_manifests import write_manifest
from register_nvt_trajectories import dcd_frames

ROOT = Path(__file__).resolve().parents[1]
# DCD frame 24 is t = 100 fs. This was 25 on the assumption that the 500 fs MD
# would reproduce the 100 fs run and its first 25 frames would duplicate v6t2.
# check_lf_overlap.py showed it does NOT -- OpenMM CUDA is not bit-reproducible
# across runs even with identical seeds -- so those frames are a genuine hole in
# each v6lf trajectory, and --first-frame 0 closes it.
FIRST_NEW = 25
N_FRAMES = 125
WINDOWS = {"v6lfa": (900, 20260901), "v6lfb": (1400, 20260902), "v6lfc": (1900, 20260903)}


def rows_for(window: str) -> list[dict[str, object]]:
    seed_frame, vseed = WINDOWS[window]
    rows = []
    for frame in range(FIRST_NEW, N_FRAMES):
        time_fs = float((frame + 1) * 4)
        rows.append({
            "window": window,
            "frame_index": frame,
            "segment_frame_index": frame - FIRST_NEW,
            "time_fs": f"{time_fs:.1f}",
            "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": str(
                ROOT.parent / f"tc_tandem_v6_lf_w{seed_frame}"
                / f"tandem_v6_w{seed_frame}_500fs_4fs.dcd"
            ),
            "convergence_status": "pending_tddft",
            "failure_class": "",
            "retry_count": 0,
            "notes": (
                f"low-frequency extension; 4 fs cadence over 104--500 fs; seeded "
                f"from v5camp frame {seed_frame}, velocity seed {vseed}; the "
                f"0--100 fs part of this window is v6t2 at 1 fs; svp profile; "
                f"full-system embedding; no interpolation"
            ),
        })
    return rows


def stamp(path: Path, trajectory: Path) -> None:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns, rows = reader.fieldnames, list(reader)
    by_index = {int(r["frame_index"]): r for r in rows}
    for index, xyz in enumerate(dcd_frames(trajectory)):
        row = by_index.get(index)
        if row is None or row["convergence_status"] == "complete":
            continue
        row["coordinate_sha256"] = hashlib.sha256(
            np.asarray(xyz, dtype="<f4").tobytes(order="C")).hexdigest()
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def main() -> None:
    import argparse
    global FIRST_NEW          # must precede any use of the name in this scope
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--first-frame", type=int, default=25,
                    help="lowest DCD frame to include (0 closes the 0-100 fs hole)")
    args = ap.parse_args()
    FIRST_NEW = args.first_frame
    manifests = ROOT / "manifests"
    for window, (seed_frame, _) in WINDOWS.items():
        traj = (ROOT.parent / f"tc_tandem_v6_lf_w{seed_frame}"
                / f"tandem_v6_w{seed_frame}_500fs_4fs.dcd")
        if not traj.exists():
            raise FileNotFoundError(f"{window}: {traj}")
        path = manifests / f"{window}_frames.csv"
        write_manifest(path, rows_for(window))
        stamp(path, traj)
        print(f"wrote {path.name}: {N_FRAMES - FIRST_NEW} rows at 4 fs "
              f"(DCD frames {FIRST_NEW}-{N_FRAMES - 1}, "
              f"t = {(FIRST_NEW + 1) * 4}-{N_FRAMES * 4} fs)")


if __name__ == "__main__":
    main()

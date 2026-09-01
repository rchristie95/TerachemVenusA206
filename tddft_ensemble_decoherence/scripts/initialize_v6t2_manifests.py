#!/usr/bin/env python3
"""Initialize the three T2-campaign manifests.

Each window is 100 frames at 1 fs, seeded from the v5 frame whose STEOM
reference already exists. Frame 0 of a window is at step 4, i.e. 1 fs, NOT
0 fs: OpenMM's DCDReporter writes after each interval, never at step 0. The
window's true t = 0 is the v5 frame it was seeded from, which is already
complete in manifests/v5camp_frames.csv -- the analysis attaches it there
rather than recomputing it.

Note the time convention. The other short windows use
`time_fs = (frame + 1) * cadence`, which is correct. v5camp instead labels
`time_ps = frame_index * 0.5`, one save interval early; harmless there because
nothing differences its timestamps, but a 1 fs cadence leaves no room for it.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

from initialize_manifests import write_manifest
from register_nvt_trajectories import dcd_frames

ROOT = Path(__file__).resolve().parents[1]

# window -> (v5 seed frame, velocity seed, |Delta| at the seed frame in cm^-1)
WINDOWS = {
    "v6t2a": (900, 20260901, 826.6),
    "v6t2b": (1400, 20260902, 58.6),
    "v6t2c": (1900, 20260903, 514.2),
}


def rows_for(window: str) -> list[dict[str, object]]:
    seed_frame, vseed, detuning = WINDOWS[window]
    rows = []
    for frame in range(100):
        time_fs = float(frame + 1)
        rows.append({
            "window": window,
            "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_fs:.1f}",
            "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": str(
                ROOT.parent / f"tc_tandem_v6_t2_w{seed_frame}"
                / f"tandem_v6_w{seed_frame}_100fs_1fs.dcd"
            ),
            "convergence_status": "pending_tddft",
            "failure_class": "",
            "retry_count": 0,
            "notes": (
                f"T2 campaign; 1 fs cadence over 1--100 fs; seeded from v5camp frame "
                f"{seed_frame} (|Delta| = {detuning:.1f} cm^-1, STEOM reference on "
                f"disk), fresh MB velocities seed {vseed}; svp profile; full-system "
                f"embedding; no interpolation"
            ),
        })
    return rows


def stamp_checksums(path: Path, trajectory: Path) -> None:
    """Bind each manifest row to the exact coordinates it will be run on.

    run_ensemble.py refuses any frame whose prepared coordinates do not hash to
    the value recorded here, so a trajectory silently regenerated under a
    half-finished manifest cannot be mistaken for the original. The hash must be
    built the same way prepare_production_frame.py builds it: canonical
    little-endian float32 interleaved XYZ.

    register_nvt_trajectories.register() does this for the older windows but
    also overwrites the notes column with its own text, so it is not reused here.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = list(reader)

    stamped = 0
    for index, xyz in enumerate(dcd_frames(trajectory)):
        if index >= len(rows):
            raise RuntimeError(f"{trajectory.name} has more frames than {path.name}")
        row = rows[index]
        if int(row["frame_index"]) != index:
            raise RuntimeError(f"Nonsequential manifest frame at row {index}")
        if row["convergence_status"] == "complete":
            continue
        canonical = np.asarray(xyz, dtype="<f4").tobytes(order="C")
        row["coordinate_sha256"] = hashlib.sha256(canonical).hexdigest()
        stamped += 1
    if stamped and stamped != len(rows):
        print(f"    note: {len(rows) - stamped} row(s) already complete, left alone")

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    manifests = ROOT / "manifests"
    for window, (seed_frame, _, _) in WINDOWS.items():
        path = manifests / f"{window}_frames.csv"
        trajectory = (
            ROOT.parent / f"tc_tandem_v6_t2_w{seed_frame}"
            / f"tandem_v6_w{seed_frame}_100fs_1fs.dcd"
        )
        if not trajectory.exists():
            raise FileNotFoundError(f"{window}: trajectory not built yet: {trajectory}")
        write_manifest(path, rows_for(window))
        stamp_checksums(path, trajectory)
        print(f"wrote {path.name}: 100 rows at 1 fs, checksums bound to {trajectory.name}")


if __name__ == "__main__":
    main()

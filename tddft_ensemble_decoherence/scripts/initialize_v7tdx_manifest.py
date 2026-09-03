#!/usr/bin/env python3
"""Initialize the dVenus-TDX campaign manifest.

500 frames at 2 ps over 1 ns, the same sampling as v5camp, on the construct with
one chromophore.

The DCD holds 600 frames, not 500: the 200 ps NPT equilibration is written to it
as well, at the same 2000-step interval. Production is DCD frames 100-599 and
the manifest indexes them directly. This was verified rather than assumed, from
the unit cell: frames 0-98 have a moving box (91.600 -> 90.457 A, barostat live)
and from frame 99 it is pinned at 90.56254 A, which is the box NVT inherited.
Taking frame 0 as production would start sampling in a box 1.04 A too large and
still relaxing -- the same class of error as the v3 stale-box trap.

Time convention: `time_ps = (frame + 1) * 2`, because OpenMM writes after each
interval and never at step 0. v5camp instead labels `time_ps = frame_index * 0.5`,
one save interval early; that is not reproduced here.

The site-B and J columns exist in the shared schema and stay empty for this
window. That is deliberate: TDX has no second chromophore, so a detuning or a
coupling written here would be indistinguishable downstream from a measured one.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

from initialize_manifests import write_manifest
from register_nvt_trajectories import dcd_frames

ROOT = Path(__file__).resolve().parents[1]
WINDOW = "v7tdx"
FRAMES = 500
FIRST_PRODUCTION_FRAME = 100     # DCD 0-99 are the NPT equilibration
CADENCE_PS = 2.0
NVT_BOX_A = 90.56254             # constant from DCD frame 99 onward
TRAJECTORY = ROOT.parent / "tc_tandem_tdx_dark2" / "tandem_tdx_dark2_1ns.dcd"


def rows_for() -> list[dict[str, object]]:
    rows = []
    for frame in range(FIRST_PRODUCTION_FRAME, FIRST_PRODUCTION_FRAME + FRAMES):
        time_ps = (frame - FIRST_PRODUCTION_FRAME + 1) * CADENCE_PS
        rows.append({
            "window": WINDOW,
            "frame_index": frame,
            "segment_frame_index": frame - FIRST_PRODUCTION_FRAME,
            "time_fs": f"{time_ps * 1000.0:.1f}",
            "time_ps": f"{time_ps:.6f}",
            "source_trajectory": str(TRAJECTORY),
            "convergence_status": "pending_tddft",
            "failure_class": "",
            "retry_count": 0,
            "notes": (
                "dVenus-TDX (Addgene #250939, Y339 deletion); 1 ns, 2 ps cadence, "
                "500 frames, matching v5camp sampling; single chromophore, so no "
                "site B, no detuning and no TDC; svp profile; full-system "
                "embedding, link-only exclusion; box neutral at 0.000000 e"
            ),
        })
    return rows


def stamp_checksums(path: Path) -> None:
    """Bind every row to the coordinates it will actually be run on.

    run_ensemble.py refuses any frame whose prepared coordinates do not hash to
    the recorded value, so a regenerated trajectory cannot be mistaken for this
    one. Built exactly as prepare_production_frame.py builds it: canonical
    little-endian float32, interleaved XYZ.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = list(reader)

    by_frame = {int(r["frame_index"]): r for r in rows}
    stamped = 0
    for index, xyz in enumerate(dcd_frames(TRAJECTORY)):
        row = by_frame.get(index)
        if row is None:                       # NPT equilibration frame, not sampled
            continue
        if row["convergence_status"] == "complete":
            continue
        row["coordinate_sha256"] = hashlib.sha256(
            np.asarray(xyz, dtype="<f4").tobytes(order="C")).hexdigest()
        stamped += 1
    if stamped != len(rows):
        raise RuntimeError(
            f"stamped {stamped} of {len(rows)} rows; the trajectory is short "
            f"(needs {FIRST_PRODUCTION_FRAME + FRAMES} frames)")

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    if not TRAJECTORY.exists():
        raise FileNotFoundError(f"trajectory not built yet: {TRAJECTORY}")
    path = ROOT / "manifests" / f"{WINDOW}_frames.csv"
    write_manifest(path, rows_for())
    stamp_checksums(path)
    print(f"wrote {path.name}: {FRAMES} rows at {CADENCE_PS} ps, "
          f"checksums bound to {TRAJECTORY.name}")


if __name__ == "__main__":
    main()

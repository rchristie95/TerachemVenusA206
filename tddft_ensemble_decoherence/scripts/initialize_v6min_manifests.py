#!/usr/bin/env python3
"""Manifests for the minimum-detuning window (v6min, v5 frame 1832).

Frame 1832 is the most degenerate configuration in the whole 500-frame v5
ensemble: |Delta| = 0.170 cm^-1 against 2|J| = 49.1. It is also the frame the
manuscript already features as its extreme delocalised case.

Full 125 frames, t = 4-500 fs at 4 fs -- nothing is skipped here, unlike v6lf,
because there is no 1 fs run covering the first 100 fs of this trajectory.

Writes a master manifest plus three BLOCK-split slices for parallel workers.
Blocks rather than interleaving: run_ensemble threads previous["A"]/["B"]
transition dipoles frame to frame for density-overlap tracking, so contiguous
blocks keep each worker's reference genuinely adjacent in time and break the
chain at only two seams. (v5camp ran interleaved slices with zero failures, so
splitting a window across workers is fine either way -- an earlier note in this
tree calling one-worker-per-window "required" overstated it.)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np

from initialize_manifests import write_manifest
from register_nvt_trajectories import dcd_frames

ROOT = Path(__file__).resolve().parents[1]
SEED_FRAME = 1832
VSEED = 20260904
N_FRAMES = 125
N_BLOCKS = 3
TRAJ = (ROOT.parent / f"tc_tandem_v6_min_w{SEED_FRAME}"
        / f"tandem_v6_w{SEED_FRAME}_500fs_4fs.dcd")


def rows_for(frames) -> list[dict[str, object]]:
    rows = []
    for frame in frames:
        time_fs = float((frame + 1) * 4)
        rows.append({
            "window": "v6min", "frame_index": frame,
            "segment_frame_index": frame,
            "time_fs": f"{time_fs:.1f}", "time_ps": f"{time_fs / 1000.0:.6f}",
            "source_trajectory": str(TRAJ),
            "convergence_status": "pending_tddft", "failure_class": "",
            "retry_count": 0,
            "notes": (
                f"minimum-detuning window; 4 fs cadence over 4--500 fs; seeded from "
                f"v5camp frame {SEED_FRAME} (|Delta| = 0.170 cm^-1, the most "
                f"degenerate configuration in the 500-frame ensemble), velocity "
                f"seed {VSEED}; STEOM reference computed for this seed; svp "
                f"profile; full-system embedding; no interpolation"
            ),
        })
    return rows


def stamp(path: Path) -> None:
    with path.open(newline="") as h:
        reader = csv.DictReader(h)
        columns, rows = reader.fieldnames, list(reader)
    by_index = {int(r["frame_index"]): r for r in rows}
    for index, xyz in enumerate(dcd_frames(TRAJ)):
        row = by_index.get(index)
        if row is None or row["convergence_status"] == "complete":
            continue
        row["coordinate_sha256"] = hashlib.sha256(
            np.asarray(xyz, dtype="<f4").tobytes(order="C")).hexdigest()
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=columns); w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    if not TRAJ.exists():
        raise FileNotFoundError(TRAJ)
    manifests = ROOT / "manifests"

    master = manifests / "v6min_frames.csv"
    write_manifest(master, rows_for(range(N_FRAMES))); stamp(master)
    print(f"wrote {master.name}: {N_FRAMES} rows at 4 fs (t = 4-500 fs)")

    edges = np.linspace(0, N_FRAMES, N_BLOCKS + 1).astype(int)
    for tag, lo, hi in zip("abc", edges[:-1], edges[1:]):
        path = manifests / f"v6min_{tag}_frames.csv"
        write_manifest(path, rows_for(range(int(lo), int(hi)))); stamp(path)
        print(f"wrote {path.name}: frames {lo}-{hi - 1} "
              f"(t = {(lo + 1) * 4}-{hi * 4} fs)")


if __name__ == "__main__":
    main()

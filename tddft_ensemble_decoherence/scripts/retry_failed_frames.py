#!/usr/bin/env python3
"""Archive seven-root failures and retry them with the nine-root TDDFT template."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--window",
        choices=("short", "long", "long100ns50"),
        required=True,
    )
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--gpu", default="0")
    args = parser.parse_args()

    manifest = ROOT / "manifests" / f"{args.window}_frames.csv"
    rows = load_rows(manifest)
    unresolved = [
        int(row["frame_index"])
        for row in rows[args.start : args.stop]
        if row["convergence_status"] not in {"complete", "failed"}
    ]
    if unresolved:
        raise RuntimeError(
            "Refusing failed-frame retries before the capped first pass is classified; "
            f"unresolved frames begin {unresolved[:10]} (total={len(unresolved)})"
        )
    failed = [
        int(row["frame_index"])
        for row in rows[args.start : args.stop]
        if row["convergence_status"] == "failed"
    ]
    print(f"[retry-plan] {args.window} failed_frames={failed}", flush=True)

    for frame in failed:
        frame_dir = ROOT / "frames" / args.window / f"frame_{frame:04d}"
        archive = ROOT / "frames" / args.window / f"frame_{frame:04d}_attempt0_7roots"
        if archive.exists() and frame_dir.exists():
            # A prior retry was interrupted; retain and resume its canonical directory.
            print(f"[retry-resume] frame={frame} archive={archive}", flush=True)
        elif archive.exists():
            raise RuntimeError(
                f"Archived attempt exists but retry directory is absent for frame {frame}: {archive}"
            )
        else:
            if not frame_dir.exists():
                raise RuntimeError(f"Missing failed-frame directory: {frame_dir}")
            frame_dir.rename(archive)
            preparation = json.loads((archive / "preparation.json").read_text())
            command = [
                sys.executable,
                str(ROOT / "scripts/prepare_production_frame.py"),
                "--topology",
                preparation["topology"],
                "--trajectory",
                preparation["trajectory"],
                "--frame-index",
                str(frame),
                "--output-dir",
                str(frame_dir),
                "--embedding-cache",
                str(ROOT / "cache/embedding_charges.npz"),
            ]
            subprocess.run(command, check=True)
            if json.loads((frame_dir / "preparation.json").read_text())["coordinate_sha256"] != preparation[
                "coordinate_sha256"
            ]:
                raise RuntimeError(f"Coordinate checksum changed while preparing retry for frame {frame}")
            print(f"[retry-archived] frame={frame} archive={archive}", flush=True)

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/run_ensemble.py"),
                "--window",
                args.window,
                "--start",
                str(frame),
                "--stop",
                str(frame + 1),
                "--gpu",
                args.gpu,
                "--roots",
                "9",
                "--retry-failed-only",
            ],
            check=True,
        )

    final_rows = load_rows(manifest)
    final = {
        int(row["frame_index"]): {
            "status": row["convergence_status"],
            "failure_class": row["failure_class"],
            "retry_count": row["retry_count"],
            "notes": row["notes"],
        }
        for row in final_rows
        if int(row["frame_index"]) in failed
    }
    print("[retry-result] " + json.dumps(final, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

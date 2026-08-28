#!/usr/bin/env python3
"""Register completed DCD sources and per-frame coordinate checksums."""

from __future__ import annotations

import argparse
import csv
import hashlib
import struct
from pathlib import Path

import numpy as np


def dcd_frames(path: Path):
    """Yield canonical interleaved little-endian float32 XYZ arrays."""
    with path.open("rb") as handle:
        raw = handle.read(4)
        if len(raw) != 4:
            raise RuntimeError(f"Empty DCD: {path}")
        little, big = struct.unpack("<i", raw)[0], struct.unpack(">i", raw)[0]
        endian = "<" if little == 84 else ">" if big == 84 else None
        if endian is None:
            raise RuntimeError(f"Unrecognized DCD header marker: {little}/{big}")
        header = handle.read(84)
        if header[:4] != b"CORD" or struct.unpack(endian + "i", handle.read(4))[0] != 84:
            raise RuntimeError(f"Malformed DCD header: {path}")
        expected_frames = struct.unpack(endian + "i", header[4:8])[0]
        title_size = struct.unpack(endian + "i", handle.read(4))[0]
        handle.seek(title_size, 1)
        if struct.unpack(endian + "i", handle.read(4))[0] != title_size:
            raise RuntimeError("Malformed DCD title record")
        if struct.unpack(endian + "i", handle.read(4))[0] != 4:
            raise RuntimeError("Malformed DCD atom-count record")
        atom_count = struct.unpack(endian + "i", handle.read(4))[0]
        if struct.unpack(endian + "i", handle.read(4))[0] != 4:
            raise RuntimeError("Malformed DCD atom-count trailer")

        frame_count = 0
        while True:
            marker_raw = handle.read(4)
            if not marker_raw:
                break
            marker = struct.unpack(endian + "i", marker_raw)[0]
            if marker in (48, 56):
                handle.seek(marker, 1)
                if struct.unpack(endian + "i", handle.read(4))[0] != marker:
                    raise RuntimeError(f"Malformed unit-cell record at frame {frame_count}")
                marker = struct.unpack(endian + "i", handle.read(4))[0]
            if marker != 4 * atom_count:
                raise RuntimeError(f"Bad coordinate record at frame {frame_count}: {marker}")
            xyz = np.empty((atom_count, 3), dtype=np.float32)
            for axis in range(3):
                if axis and struct.unpack(endian + "i", handle.read(4))[0] != marker:
                    raise RuntimeError(f"Bad coordinate marker at frame {frame_count}")
                values = np.fromfile(handle, dtype=endian + "f4", count=atom_count)
                if values.size != atom_count:
                    raise RuntimeError(f"Truncated coordinates at frame {frame_count}")
                xyz[:, axis] = values
                if struct.unpack(endian + "i", handle.read(4))[0] != marker:
                    raise RuntimeError(f"Bad coordinate trailer at frame {frame_count}")
            frame_count += 1
            yield np.asarray(xyz, dtype="<f4")
        if frame_count != expected_frames:
            raise RuntimeError(
                f"DCD header/data frame mismatch: {expected_frames} versus {frame_count}"
            )


def register(manifest: Path, trajectory: Path, expected_frames: int) -> None:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = list(reader)
    if columns is None or len(rows) != expected_frames:
        raise RuntimeError(f"Unexpected manifest length for {manifest}: {len(rows)}")

    seen = 0
    for index, xyz in enumerate(dcd_frames(trajectory)):
        if index >= len(rows):
            raise RuntimeError(f"Trajectory has too many frames for {manifest}")
        row = rows[index]
        if int(row["frame_index"]) != index:
            raise RuntimeError(f"Nonsequential manifest frame at row {index}")
        if row["convergence_status"] in {"blocked_missing_source", "pending_tddft"}:
            row["source_trajectory"] = str(trajectory.resolve())
            row["coordinate_sha256"] = hashlib.sha256(xyz.tobytes(order="C")).hexdigest()
            row["convergence_status"] = "pending_tddft"
            row["failure_class"] = ""
            row["notes"] = (
                "New exact-AMBER CR2 NVT; coordinate hash is canonical "
                "little-endian float32 interleaved XYZ"
            )
        seen += 1
    if seen != expected_frames:
        raise RuntimeError(f"Unexpected trajectory length for {trajectory}: {seen}")

    temporary = manifest.with_suffix(manifest.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--short-dcd", type=Path)
    parser.add_argument("--long-dcd", type=Path)
    parser.add_argument("--long100ns50-dcd", type=Path)
    parser.add_argument("--short100fs50-dcd", type=Path)
    parser.add_argument("--short100fs50-fixedmm-dcd", type=Path)
    parser.add_argument("--short100fs100-wb97xd3-dcd", type=Path)
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "manifests",
    )
    args = parser.parse_args()
    selected = [
        ("short", args.short_dcd, 2000),
        ("long", args.long_dcd, 1000),
        ("long100ns50", args.long100ns50_dcd, 50),
        ("short100fs50", args.short100fs50_dcd, 50),
        ("short100fs50_fixedmm", args.short100fs50_fixedmm_dcd, 50),
        ("short100fs100_wb97xd3", args.short100fs100_wb97xd3_dcd, 100),
    ]
    selected = [item for item in selected if item[1] is not None]
    if not selected:
        parser.error("provide at least one trajectory argument")
    for window, trajectory, expected in selected:
        register(
            args.manifest_dir / f"{window}_frames.csv",
            trajectory,
            expected,
        )
        print(f"Registered {expected} {window} trajectory frames.")


if __name__ == "__main__":
    main()

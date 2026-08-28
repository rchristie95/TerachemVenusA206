#!/usr/bin/env python3
"""Stream selected zero-based models from a large multi-model PDB."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("frames", nargs="+", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    wanted = set(args.frames)
    if not wanted or min(wanted) < 0:
        raise ValueError("Frame indices must be non-negative")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    headers: list[str] = []
    current = -1
    handle = None
    found: set[int] = set()
    try:
        with args.source.open(errors="replace") as source:
            for line in source:
                if line.startswith("MODEL"):
                    current += 1
                    if current in wanted:
                        path = args.output_dir / f"frame_{current:04d}.pdb"
                        if path.exists():
                            raise FileExistsError(f"Refusing to overwrite {path}")
                        handle = path.open("w")
                        handle.writelines(headers)
                        handle.write(line)
                        found.add(current)
                    continue
                if current < 0:
                    headers.append(line)
                elif handle is not None:
                    handle.write(line)
                    if line.startswith("ENDMDL"):
                        handle.write("END\n")
                        handle.close()
                        handle = None
                if current > max(wanted) and handle is None:
                    break
    finally:
        if handle is not None:
            handle.close()
    missing = wanted - found
    if missing:
        raise RuntimeError(f"Missing requested frames: {sorted(missing)}")
    print(f"extracted {len(found)} frames: {sorted(found)}")


if __name__ == "__main__":
    main()

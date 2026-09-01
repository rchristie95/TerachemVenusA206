#!/usr/bin/env python3
"""Did the 500 fs MD reproduce the first 100 fs of the matching 1 fs run?

Same seed frame, same velocity and integrator seeds, same system -- so the
trajectories should be identical and only the DCD save interval differs. That is
worth checking rather than assuming: if it holds, the 4 fs and 1 fs records
splice into one continuous 500 fs window per seed; if not, each stands alone and
only the >100 fs frames of the long run are used.

Compares DCD frame 24 of the 500 fs run (t = 100 fs) with frame 99 of the 100 fs
run (t = 100 fs), by the same canonical little-endian float32 hash the manifests
are bound to.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from register_nvt_trajectories import dcd_frames

REPO = Path("/home/robson/PetaChem")
SEEDS = (900, 1400, 1900)


def frame_hash(path: Path, index: int) -> str:
    for i, xyz in enumerate(dcd_frames(path)):
        if i == index:
            return hashlib.sha256(np.asarray(xyz, dtype="<f4").tobytes(order="C")).hexdigest()
    raise IndexError(f"{path} has no frame {index}")


def main() -> None:
    ok = True
    for seed in SEEDS:
        short = REPO / f"tc_tandem_v6_t2_w{seed}/tandem_v6_w{seed}_100fs_1fs.dcd"
        long = REPO / f"tc_tandem_v6_lf_w{seed}/tandem_v6_w{seed}_500fs_4fs.dcd"
        a, b = frame_hash(short, 99), frame_hash(long, 24)
        match = a == b
        ok &= match
        print(f"  seed {seed}: t=100 fs  1fs-run {a[:16]}  4fs-run {b[:16]}  "
              f"{'MATCH -- splice' if match else 'DIFFER -- keep separate'}")
    print("  => trajectories are reproducible; the 1 fs and 4 fs records are one window"
          if ok else
          "  => not reproducible; treat the 500 fs runs as independent windows")


if __name__ == "__main__":
    main()

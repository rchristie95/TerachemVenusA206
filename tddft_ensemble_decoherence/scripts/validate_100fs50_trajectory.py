#!/usr/bin/env python3
"""Validate the 50-frame, 2 fs cadence NVT trajectory before QM extraction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from openmm.app import PDBFile

from prepare_production_frame import read_dcd_frame
from register_nvt_trajectories import dcd_frames


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
TRAJECTORY = (
    REPO / "tc_tandem_nvt_100fs_50frames_20260723/tandem_nvt_100fs_50frames.dcd"
)
TOPOLOGY = REPO / "tc_tandem_nvt_100fs_50frames_20260723/solvated_protonated.pdb"
OUT = ROOT / "rapid_100fs50/trajectory_provenance.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with TRAJECTORY.open("rb") as handle:
        if struct.unpack("<i", handle.read(4))[0] != 84:
            raise RuntimeError("Expected little-endian DCD")
        raw = handle.read(84)
    frames, start, interval, final = struct.unpack("<4i", raw[4:20])
    if (frames, start, interval, final) != (50, 1, 1, 50):
        raise RuntimeError(
            f"Unexpected DCD header: {(frames, start, interval, final)}"
        )
    topology_atoms = PDBFile(str(TOPOLOGY)).topology.getNumAtoms()
    if topology_atoms != 69_456:
        raise RuntimeError(f"Unexpected topology atom count: {topology_atoms}")
    counted = sum(1 for _ in dcd_frames(TRAJECTORY))
    if counted != 50:
        raise RuntimeError(f"Expected 50 coordinate records, found {counted}")
    first, first_box, first_meta = read_dcd_frame(TRAJECTORY, 0)
    last, last_box, last_meta = read_dcd_frame(TRAJECTORY, 49)
    if first_meta["atom_count"] != topology_atoms or last_meta["atom_count"] != topology_atoms:
        raise RuntimeError("Topology/DCD atom-count mismatch")
    if not np.isfinite(first).all() or not np.isfinite(last).all():
        raise RuntimeError("Nonfinite endpoint coordinates")
    if np.any(first_box <= 0.0) or np.any(last_box <= 0.0):
        raise RuntimeError("Invalid periodic box")
    times_fs = [float(2 * (frame + 1)) for frame in range(50)]
    result = {
        "status": "validated",
        "trajectory": str(TRAJECTORY.resolve()),
        "trajectory_sha256": sha256(TRAJECTORY),
        "topology": str(TOPOLOGY.resolve()),
        "topology_sha256": sha256(TOPOLOGY),
        "topology_atoms": topology_atoms,
        "dcd": {
            "frames": frames,
            "start_step": start,
            "save_interval_steps": interval,
            "final_step": final,
        },
        "absolute_source_checkpoint_step": 500_000,
        "relative_times_fs": times_fs,
        "coverage_fs": [2.0, 100.0],
        "sampling_interval_fs": 2.0,
        "periodic_box_first_A": first_box.tolist(),
        "periodic_box_last_A": last_box.tolist(),
        "mm_system": "intact full tandem; no QM cuts or link atoms in MD",
        "source_checkpoint": str(
            (REPO / "tc_tandem_nvt_rerun_20260722_retry2/nvt.chk").resolve()
        ),
        "steom_reference_energy_eV": 2.366633,
        "cap_construction": (
            "instantaneous severed-bond vectors; link atoms 42--44 excluded from TDC"
        ),
        "no_interpolation_or_frame_replacement": True,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate the intact-MM 50-frame, 2--100 ns DCD before QM extraction."""

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
    REPO
    / "tc_tandem_nvt_100ns_50frames_20260723"
    / "tandem_nvt_100ns_50frames.dcd"
)
TOPOLOGY = (
    REPO
    / "tc_tandem_nvt_100ns_50frames_20260723"
    / "solvated_protonated.pdb"
)
OUT = ROOT / "rapid_100ns50" / "trajectory_provenance.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def header(path: Path) -> dict[str, int | float]:
    with path.open("rb") as handle:
        if struct.unpack("<i", handle.read(4))[0] != 84:
            raise RuntimeError("Expected little-endian DCD")
        raw = handle.read(84)
    frames, start_step, interval, final_step = struct.unpack("<4i", raw[4:20])
    return {
        "frames": frames,
        "start_step": start_step,
        "save_interval_steps": interval,
        "final_step": final_step,
        "dcd_delta": struct.unpack("<f", raw[40:44])[0],
    }


def main() -> None:
    metadata = header(TRAJECTORY)
    expected = {
        "frames": 50,
        "start_step": 1_000_000,
        "save_interval_steps": 1_000_000,
        "final_step": 50_000_000,
    }
    for key, value in expected.items():
        if metadata[key] != value:
            raise RuntimeError(f"DCD {key}: expected {value}, got {metadata[key]}")

    topology = PDBFile(str(TOPOLOGY))
    topology_atoms = topology.topology.getNumAtoms()
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
    if np.any(first_box <= 0) or np.any(last_box <= 0):
        raise RuntimeError("Invalid endpoint periodic box")

    times_ps = [
        (metadata["start_step"] + frame * metadata["save_interval_steps"]) * 0.002
        for frame in range(50)
    ]
    if times_ps != [float(value) for value in range(2000, 100001, 2000)]:
        raise RuntimeError("Unexpected physical-time grid")

    result = {
        "status": "validated",
        "trajectory": str(TRAJECTORY.resolve()),
        "trajectory_sha256": sha256(TRAJECTORY),
        "topology": str(TOPOLOGY.resolve()),
        "topology_sha256": sha256(TOPOLOGY),
        "topology_atoms": topology_atoms,
        "dcd": metadata,
        "timestep_fs": 2.0,
        "times_ps": times_ps,
        "coverage_ns": [2.0, 100.0],
        "sampling_interval_ns": 2.0,
        "periodic_box_first_A": first_box.tolist(),
        "periodic_box_last_A": last_box.tolist(),
        "mm_system": "intact full tandem; no QM cuts or link atoms in MD",
        "source_checkpoint": str(
            (
                REPO
                / "tc_tandem_nvt_rerun_20260722_retry2"
                / "nvt.chk"
            ).resolve()
        ),
        "source_checkpoint_step": 500_000,
        "steom_reference_energy_eV": 2.366633,
        "tddft_calibration": (
            "Calculated after MM completion with the selected rapid TDDFT profile "
            "on the identical retained reference geometry and embedding"
        ),
        "cap_construction": "instantaneous severed-bond vectors; link atoms 42--44 excluded from TDC",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

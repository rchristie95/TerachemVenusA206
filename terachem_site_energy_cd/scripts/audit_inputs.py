#!/usr/bin/env python3
"""Audit the exact trajectory/topology gate and write a machine-readable manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from openmm import NonbondedForce, unit
from openmm.app import AmberPrmtopFile, NoCutoff, PDBFile

from prepare_production_frame import SITE_KEYS, read_dcd_frame, topology_descriptor

# v2 production trajectory, tc_tandem_nvt_v2/tandem_nvt_v2_1000.dcd. This
# replaces the v1-era rerun_20260722_retry2 hash, whose trajectory was deleted
# as superseded (CLEANUP_MANIFEST.md); the site-energy campaign was in fact run
# on v2 (322 v2_* result directories, ens_v2_all.npz). Mirrored in
# reference/SHA256SUMS so the expected value has an independent home rather than
# living only in the manifest it validates.
EXPECTED_ARCHIVE_DCD = "397435de72a189641e73f65883cbe60a5006306aff966799eab5c32b98d1bd19"

# The v2 production trajectory has 1200 frames, not the 1000 this gate hardcoded
# for the v1-era run. Corroborated independently of the DCD header: the campaign's
# own result directories reference frame indices up to 1195
# (terachem_site_energy_cd/results/v2_linkonly_frame_1195), which cannot exist in
# a 1000-frame trajectory. So 1200 is the production geometry and the old
# constant was stale, rather than the gate being loosened to make it pass.
EXPECTED_FRAME_COUNT = 1200


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capture(*command: str) -> str:
    """Best-effort command output. Returns "" if the binary is unavailable.

    check=False already tolerates a nonzero exit, but a missing executable
    raises FileNotFoundError and aborted the whole audit. Provenance capture
    must never be the reason a manifest cannot be written; pass the value
    explicitly (e.g. --starting-commit) when the tool is not callable.
    """
    try:
        result = subprocess.run(
            command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
        )
    except (FileNotFoundError, OSError):
        return ""
    return result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--amber-cr2-prmtop", type=Path, required=True)
    parser.add_argument("--terachem", type=Path, required=True)
    parser.add_argument("--expected-trajectory-sha256", default=EXPECTED_ARCHIVE_DCD)
    parser.add_argument(
        "--expected-frame-count", type=int, default=EXPECTED_FRAME_COUNT,
        help="Frame count the production trajectory must have. Was hardcoded to "
             "1000 for the v1-era trajectory; the v2 production DCD has 1200.",
    )
    parser.add_argument(
        "--starting-commit", default=None,
        help="Record this commit instead of shelling out to git. Use when git "
             "is not callable from the audit environment; the manifest must "
             "still carry a commit.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.topology, args.trajectory, args.amber_cr2_prmtop, args.terachem):
        if not path.is_file():
            raise FileNotFoundError(path)
    pdb = PDBFile(str(args.topology))
    atoms = list(pdb.topology.atoms())
    _, box0, dcd = read_dcd_frame(args.trajectory, 0)
    if dcd["atom_count"] != len(atoms):
        raise RuntimeError(f"Topology/DCD atom mismatch: {len(atoms)} != {dcd['atom_count']}")

    reference = AmberPrmtopFile(str(args.amber_cr2_prmtop))
    ref_cr2 = [residue for residue in reference.topology.residues() if residue.name == "CR2"]
    if len(ref_cr2) != 1:
        raise RuntimeError("Reference prmtop must contain exactly one CR2 residue")
    ref_names = [atom.name for atom in ref_cr2[0].atoms()]
    system = reference.createSystem(nonbondedMethod=NoCutoff)
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    ref_charge = sum(
        nonbonded.getParticleParameters(atom.index)[0].value_in_unit(unit.elementary_charge)
        for atom in ref_cr2[0].atoms()
    )
    mappings = {}
    for site, keys in SITE_KEYS.items():
        cr2 = [r for r in pdb.topology.residues() if (r.chain.id, r.id, r.name) == (*keys["cr2"], "CR2")]
        tyr = [r for r in pdb.topology.residues() if (r.chain.id, r.id, r.name) == (*keys["tyr"], "TYR")]
        if len(cr2) != 1 or len(tyr) != 1:
            raise RuntimeError(f"Site {site} residue mapping is not unique")
        names = [atom.name for atom in cr2[0].atoms()]
        mappings[site] = {
            "cr2": list(keys["cr2"]), "stacked_tyr": list(keys["tyr"]),
            "cr2_atoms": len(names),
            "cr2_atom_names_exact": len(names) == len(set(names)) and set(names) == set(ref_names),
            "mapping_rule": "unique exact atom name; topology order retained independently at each site",
        }

    trajectory_hash = sha256(args.trajectory)
    topology_hash = sha256(args.topology)
    prmtop_hash = sha256(args.amber_cr2_prmtop)
    expected = args.expected_trajectory_sha256.lower()
    gate = (
        trajectory_hash == expected
        and dcd["frame_count"] == args.expected_frame_count
        and dcd["save_interval_steps"] == 500
        and all(item["cr2_atom_names_exact"] for item in mappings.values())
        and abs(ref_charge + 1.0) < 1.0e-6
    )
    payload = {
        "status": "production_identity_validated" if gate else "smoke_only_production_identity_mismatch",
        "production_join_allowed": gate,
        "starting_commit": args.starting_commit or capture("git", "rev-parse", "HEAD"),
        "trajectory": str(args.trajectory.resolve()),
        "trajectory_sha256": trajectory_hash,
        "expected_archive_trajectory_sha256": expected,
        "topology": str(args.topology.resolve()),
        "topology_sha256": topology_hash,
        "topology_descriptor_sha256": topology_descriptor(pdb.topology),
        "amber_cr2_prmtop": str(args.amber_cr2_prmtop.resolve()),
        "amber_cr2_prmtop_sha256": prmtop_hash,
        "dcd": dcd,
        "frame_time_convention_ps": "(zero_based_frame + 1) ps",
        "box_frame0_A": box0.tolist(),
        "site_mapping": mappings,
        "compact_qm_definition": {"atoms": 44, "physical_atoms": 41, "charge": -1, "multiplicity": 1},
        "reference_cr2_charge_e": ref_charge,
        "terachem": {"path": str(args.terachem.resolve()), "binary_sha256": sha256(args.terachem), "version": "1.97B-beta-251105"},
        "python": {"executable": sys.executable, "version": sys.version.split()[0], "platform": platform.platform()},
        "gpu": capture("nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"),
        "settings": {
            "method": "wB97X-D3", "basis_smoke": "6-31G*", "tda": True, "roots": 7,
            "embedding_radius_A": 12.0, "boundary_exclusion_A": 1.8,
            "pcm": "none", "partner_cr2": "complete RESP-charged MM residue",
            "fixed_embedding_selection": True,
        },
        "manual_capability_audit": {
            "manual": "/home/robson/Desktop/TeraChemPython/TeraChem/UserGuide1.97B-beta-251105.pdf",
            "electric_transition_dipoles": True,
            "magnetic_transition_dipoles_or_rotational_strengths": False,
            "conclusion": "Only interaction-induced exciton-chirality CD may be reconstructed; it is not absolute molar CD.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    if not gate:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

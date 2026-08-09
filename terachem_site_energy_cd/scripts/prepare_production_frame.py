#!/usr/bin/env python3
"""Prepare both instantaneous 44-atom QM/MM site jobs from one production DCD frame."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import tempfile

import numpy as np
from scipy.spatial import cKDTree
from openmm import NonbondedForce, unit
from openmm.app import AmberPrmtopFile, ForceField, Modeller, NoCutoff, PDBFile

from mapping import (
    image_near,
    kabsch,
    minimum_image,
    reference_mapping,
    residue_by_key,
    unwrap,
)


SITE_KEYS = {
    "A": {"cr2": ("A", "66"), "tyr": ("A", "202")},
    "B": {"cr2": ("A", "328"), "tyr": ("A", "464")},
}

EMBEDDING_POLICY_VERSION = 2
BOUNDARY_REDISTRIBUTION_POLICY = (
    "equal_increment_over_retained_non_qm_atoms_in_each_partially_selected_residue"
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_sha256(payload: dict) -> str:
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(canonical).hexdigest()


def atomic_write_json(path: Path, payload: dict) -> None:
    """Durably replace a JSON object without exposing a partial marker."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            json.dump(payload, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
        temporary_name = None
        try:
            directory_fd = os.open(path.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except OSError:
            # The file itself is fsynced and replacement remains atomic on
            # filesystems that do not implement directory fsync.
            pass
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def embedding_policy(
    *,
    retain_partner_cr2_charges: bool,
    conserve_boundary_residue_charge: bool,
) -> dict:
    """Return the complete, hashable electrostatic-embedding policy."""

    return {
        "version": EMBEDDING_POLICY_VERSION,
        "retain_partner_cr2_charges": bool(retain_partner_cr2_charges),
        "partner_cr2_selection": (
            "retain_complete_non_qm_partner_residue_regardless_of_radius"
            if retain_partner_cr2_charges
            else "legacy_radius_and_boundary_selection"
        ),
        "partner_cr2_charge_source": (
            "configured_amber_prmtop_exact_atom_name_mapping"
            if retain_partner_cr2_charges
            else "zeroed_as_in_legacy_embedding"
        ),
        "conserve_boundary_residue_charge": bool(
            conserve_boundary_residue_charge
        ),
        "boundary_charge_redistribution": (
            BOUNDARY_REDISTRIBUTION_POLICY
            if conserve_boundary_residue_charge
            else "disabled"
        ),
    }


def validate_fixed_embedding_selection(
    selection: dict,
    *,
    topology_descriptor_sha256: str,
    embedding_policy_payload: dict,
    embedding_policy_sha256: str,
    amber_cr2_prmtop: Path | None,
    amber_cr2_prmtop_sha256: str | None,
    radius_A: float,
    boundary_exclusion_A: float,
    atom_count: int,
) -> None:
    """Fail closed when reusing the frame-zero MM atom inventory."""

    if not isinstance(selection, dict):
        raise RuntimeError("Fixed MM selection is not a JSON object")
    if selection.get("status") != "fixed_embedding_selection":
        raise RuntimeError("Fixed MM selection has an invalid status")
    if selection.get("source_frame_index") != 0:
        raise RuntimeError(
            "Fixed MM selection must originate from logical frame 0"
        )
    if (
        selection.get("topology_descriptor_sha256")
        != topology_descriptor_sha256
    ):
        raise RuntimeError("Fixed MM selection topology does not match this frame")

    saved_policy = selection.get("embedding_policy")
    if saved_policy is None:
        if (
            embedding_policy_payload["retain_partner_cr2_charges"]
            or embedding_policy_payload[
                "conserve_boundary_residue_charge"
            ]
        ):
            raise RuntimeError(
                "Fixed MM selection predates the requested embedding policy"
            )
    elif (
        saved_policy != embedding_policy_payload
        or selection.get("embedding_policy_sha256")
        != embedding_policy_sha256
    ):
        raise RuntimeError("Fixed MM-selection embedding policy mismatch")

    saved_prmtop = selection.get("amber_cr2_prmtop")
    saved_prmtop_hash = selection.get("amber_cr2_prmtop_sha256")
    if amber_cr2_prmtop is not None and (
        saved_prmtop is None
        or Path(saved_prmtop).resolve() != amber_cr2_prmtop.resolve()
        or saved_prmtop_hash != amber_cr2_prmtop_sha256
    ):
        raise RuntimeError(
            "Fixed MM-selection CR2 AMBER prmtop provenance mismatch"
        )

    scalar_checks = (
        ("embedding radius", "radius_A_at_selection", radius_A),
        (
            "embedding boundary exclusion",
            "boundary_exclusion_A_at_selection",
            boundary_exclusion_A,
        ),
    )
    for label, key, expected in scalar_checks:
        try:
            saved = float(selection[key])
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(
                f"Fixed MM selection lacks a valid {label}"
            ) from error
        if not np.isfinite(saved) or not np.isclose(
            saved, expected, atol=1.0e-12, rtol=0.0
        ):
            raise RuntimeError(
                f"Fixed MM-selection {label} mismatch: "
                f"{saved!r} != {expected!r}"
            )

    sites = selection.get("sites")
    if not isinstance(sites, dict) or set(sites) != {"A", "B"}:
        raise RuntimeError(
            "Fixed MM selection must contain exactly sites A and B"
        )
    for site, record in sites.items():
        if not isinstance(record, dict):
            raise RuntimeError(
                f"Fixed MM selection site {site} is not an object"
            )
        indices = record.get("atom_indices")
        if (
            not isinstance(indices, list)
            or not indices
            or any(
                isinstance(index, bool) or not isinstance(index, int)
                for index in indices
            )
        ):
            raise RuntimeError(
                f"Fixed MM selection site {site} has invalid atom indices"
            )
        if len(indices) != len(set(indices)):
            raise RuntimeError(
                f"Fixed MM selection site {site} has duplicate atom indices"
            )
        if any(index < 0 or index >= atom_count for index in indices):
            raise RuntimeError(
                f"Fixed MM selection site {site} has out-of-range atom indices"
            )
        if record.get("atom_count") != len(indices):
            raise RuntimeError(
                f"Fixed MM selection site {site} atom count is inconsistent"
            )
        try:
            net_charge = float(record["net_charge_e"])
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(
                f"Fixed MM selection site {site} has invalid net charge"
            ) from error
        if not np.isfinite(net_charge):
            raise RuntimeError(
                f"Fixed MM selection site {site} has nonfinite net charge"
            )


def topology_descriptor(topology) -> str:
    digest = hashlib.sha256()
    for atom in topology.atoms():
        digest.update(
            f"{atom.index}|{atom.residue.chain.id}|{atom.residue.id}|{atom.residue.name}|"
            f"{atom.name}|{atom.element.symbol if atom.element else '?'}\n".encode()
        )
    return digest.hexdigest()


def read_dcd_frame(path: Path, frame_index: int) -> tuple[np.ndarray, np.ndarray, dict]:
    with path.open("rb") as handle:
        raw = handle.read(4)
        little, big = struct.unpack("<i", raw)[0], struct.unpack(">i", raw)[0]
        endian = "<" if little == 84 else ">" if big == 84 else None
        if endian is None:
            raise RuntimeError(f"Unrecognized DCD header: {little}/{big}")
        header = handle.read(84)
        if header[:4] != b"CORD" or struct.unpack(endian + "i", handle.read(4))[0] != 84:
            raise RuntimeError("Malformed DCD header")
        frame_count, start_step, save_interval = struct.unpack(endian + "3i", header[4:16])
        if not 0 <= frame_index < frame_count:
            raise IndexError(f"Frame {frame_index} outside 0..{frame_count - 1}")
        title_size = struct.unpack(endian + "i", handle.read(4))[0]
        handle.seek(title_size, 1)
        if struct.unpack(endian + "i", handle.read(4))[0] != title_size:
            raise RuntimeError("Malformed DCD title")
        if struct.unpack(endian + "i", handle.read(4))[0] != 4:
            raise RuntimeError("Malformed DCD atom count")
        atom_count = struct.unpack(endian + "i", handle.read(4))[0]
        if struct.unpack(endian + "i", handle.read(4))[0] != 4:
            raise RuntimeError("Malformed DCD atom trailer")
        data_start = handle.tell()
        first_marker = struct.unpack(endian + "i", handle.read(4))[0]
        unit_cell_size = first_marker if first_marker in (48, 56) else 0
        frame_size = (8 + unit_cell_size if unit_cell_size else 0) + 3 * (8 + 4 * atom_count)
        handle.seek(data_start + frame_index * frame_size)
        marker = struct.unpack(endian + "i", handle.read(4))[0]
        if unit_cell_size:
            if marker != unit_cell_size:
                raise RuntimeError(f"Unit-cell marker mismatch at frame {frame_index}")
            cell_raw = np.frombuffer(handle.read(marker), dtype=endian + "f8").copy()
            if struct.unpack(endian + "i", handle.read(4))[0] != marker:
                raise RuntimeError("Malformed unit-cell trailer")
            marker = struct.unpack(endian + "i", handle.read(4))[0]
        else:
            raise RuntimeError("Production trajectory must contain periodic box records")
        if marker != 4 * atom_count:
            raise RuntimeError("Malformed X coordinate record")
        xyz = np.empty((atom_count, 3), dtype=np.float64)
        for axis in range(3):
            if axis and struct.unpack(endian + "i", handle.read(4))[0] != marker:
                raise RuntimeError("Malformed coordinate marker")
            xyz[:, axis] = np.fromfile(handle, dtype=endian + "f4", count=atom_count)
            if struct.unpack(endian + "i", handle.read(4))[0] != marker:
                raise RuntimeError("Malformed coordinate trailer")
    if len(cell_raw) != 6 or not np.allclose(cell_raw[[1, 3, 4]], 0.0, atol=1.0e-8):
        raise RuntimeError(f"Only orthorhombic production boxes are supported: {cell_raw}")
    box = np.asarray([cell_raw[0], cell_raw[2], cell_raw[5]], dtype=float)
    metadata = {
        "frame_count": int(frame_count),
        "start_step": int(start_step),
        "save_interval_steps": int(save_interval),
        "atom_count": int(atom_count),
    }
    return xyz, box, metadata


def load_cr2_resp_charges(prmtop_path: Path) -> dict[str, float]:
    """Read the configured CR2 RESP charges using exact atom-name identity."""

    prmtop = AmberPrmtopFile(str(prmtop_path))
    residues = [
        residue
        for residue in prmtop.topology.residues()
        if residue.name == "CR2"
    ]
    if len(residues) != 1:
        raise RuntimeError(
            f"Expected exactly one CR2 template in {prmtop_path}; "
            f"found {len(residues)}"
        )
    atoms = list(residues[0].atoms())
    names = [atom.name for atom in atoms]
    if len(names) != len(set(names)):
        raise RuntimeError(
            f"CR2 atom names are not unique in {prmtop_path}: {names}"
        )
    system = prmtop.createSystem(nonbondedMethod=NoCutoff)
    nonbonded = next(
        force
        for force in system.getForces()
        if isinstance(force, NonbondedForce)
    )
    charges = {
        atom.name: float(
            nonbonded.getParticleParameters(atom.index)[0].value_in_unit(
                unit.elementary_charge
            )
        )
        for atom in atoms
    }
    values = np.asarray(list(charges.values()), dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError(f"Nonfinite CR2 RESP charge in {prmtop_path}")
    return charges


def load_or_build_charges(
    pdb: PDBFile,
    cache: Path,
    *,
    retain_partner_cr2_charges: bool = False,
    conserve_boundary_residue_charge: bool = False,
    amber_cr2_prmtop: Path | None = None,
) -> np.ndarray:
    descriptor = topology_descriptor(pdb.topology)
    policy = embedding_policy(
        retain_partner_cr2_charges=retain_partner_cr2_charges,
        conserve_boundary_residue_charge=conserve_boundary_residue_charge,
    )
    policy_hash = json_sha256(policy)
    legacy_policy = not (
        retain_partner_cr2_charges or conserve_boundary_residue_charge
    )
    if retain_partner_cr2_charges and amber_cr2_prmtop is None:
        raise ValueError(
            "--amber-cr2-prmtop is required when retaining partner CR2 charges"
        )
    prmtop_hash = (
        file_sha256(amber_cr2_prmtop)
        if retain_partner_cr2_charges and amber_cr2_prmtop is not None
        else None
    )
    forcefield_label = (
        "amber14-all.xml + amber14/tip3p.xml; CR2 RESP from configured AMBER "
        "prmtop by exact atom name"
        if retain_partner_cr2_charges
        else "amber14-all.xml + amber14/tip3p.xml; CR2 wholly QM"
    )
    if cache.exists():
        with np.load(cache, allow_pickle=False) as saved:
            if str(saved["topology_descriptor_sha256"]) != descriptor:
                raise RuntimeError(
                    "Embedding-cache topology descriptor mismatch"
                )
            charges = np.asarray(saved["charges_e"], dtype=float)
            if len(charges) != pdb.topology.getNumAtoms():
                raise RuntimeError("Embedding-cache atom-count mismatch")
            if str(saved["forcefield"]) != forcefield_label:
                raise RuntimeError(
                    "Embedding-cache force field does not match the production "
                    f"model: {saved['forcefield']!s} versus {forcefield_label}"
                )
            if "embedding_policy_sha256" in saved.files:
                if str(saved["embedding_policy_sha256"]) != policy_hash:
                    raise RuntimeError(
                        "Embedding-cache electrostatic policy mismatch"
                    )
            elif not legacy_policy:
                raise RuntimeError(
                    "Embedding cache predates the requested electrostatic "
                    "policy; configure a fresh embedding-cache path"
                )
            if retain_partner_cr2_charges:
                if (
                    "amber_cr2_prmtop_sha256" not in saved.files
                    or str(saved["amber_cr2_prmtop_sha256"]) != prmtop_hash
                ):
                    raise RuntimeError(
                        "Embedding-cache CR2 AMBER prmtop hash mismatch"
                    )
        return charges

    original_atoms = list(pdb.topology.atoms())
    retained_atoms = [atom for atom in original_atoms if atom.residue.name != "CR2"]
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.delete([atom for atom in original_atoms if atom.residue.name == "CR2"])
    modelled_atoms = list(modeller.topology.atoms())
    if len(modelled_atoms) != len(retained_atoms):
        raise RuntimeError("CR2 deletion changed the retained MM inventory")
    for old, new in zip(retained_atoms, modelled_atoms):
        old_key = (old.residue.chain.id, old.residue.id, old.residue.name, old.name)
        new_key = (new.residue.chain.id, new.residue.id, new.residue.name, new.name)
        if old_key != new_key:
            raise RuntimeError(f"MM atom order changed: {old_key} versus {new_key}")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=NoCutoff, ignoreExternalBonds=True
    )
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    retained_charges = np.asarray(
        [
            nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
            for index in range(system.getNumParticles())
        ],
        dtype=float,
    )
    charges = np.zeros(len(original_atoms), dtype=float)
    charges[[atom.index for atom in retained_atoms]] = retained_charges
    if retain_partner_cr2_charges:
        assert amber_cr2_prmtop is not None
        cr2_resp = load_cr2_resp_charges(amber_cr2_prmtop)
        production_cr2 = [
            residue
            for residue in pdb.topology.residues()
            if residue.name == "CR2"
        ]
        if len(production_cr2) != 2:
            raise RuntimeError(
                "Partner-CR2 embedding requires exactly two CR2 residues in "
                f"the production topology; found {len(production_cr2)}"
            )
        expected_names = set(cr2_resp)
        for residue in production_cr2:
            atoms_by_name = {atom.name: atom for atom in residue.atoms()}
            if (
                len(atoms_by_name) != len(list(residue.atoms()))
                or set(atoms_by_name) != expected_names
            ):
                raise RuntimeError(
                    f"CR2 atom inventory at {residue.chain.id}:{residue.id} "
                    "does not exactly match the configured AMBER prmtop"
                )
            for name, charge in cr2_resp.items():
                charges[atoms_by_name[name].index] = charge
    cache.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b",
            dir=cache.parent,
            prefix=f".{cache.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            np.savez_compressed(
                handle,
                charges_e=charges,
                topology_descriptor_sha256=np.asarray(descriptor),
                forcefield=np.asarray(forcefield_label),
                embedding_policy_sha256=np.asarray(policy_hash),
                amber_cr2_prmtop_sha256=np.asarray(prmtop_hash or ""),
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, cache)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return charges


def redistribute_residue_boundary_charge(
    charges: np.ndarray,
    keep: np.ndarray,
    *,
    enabled: bool,
) -> tuple[np.ndarray, dict]:
    """Return selected charges, optionally conserving a partial residue's sum."""

    source = np.asarray(charges, dtype=float)
    mask = np.asarray(keep, dtype=bool)
    if source.ndim != 1 or mask.shape != source.shape:
        raise ValueError("Residue charge and selection arrays must be 1-D peers")
    retained_count = int(mask.sum())
    excluded_count = int((~mask).sum())
    if retained_count == 0:
        return source[mask], {
            "retained_atoms": 0,
            "excluded_atoms": excluded_count,
            "excluded_charge_e": float(source.sum()),
            "equal_increment_per_retained_atom_e": 0.0,
            "charge_conservation_error_e": 0.0,
            "redistributed": False,
        }
    selected = source[mask].copy()
    excluded_charge = float(source[~mask].sum())
    increment = 0.0
    redistributed = bool(enabled and excluded_count)
    if redistributed:
        increment = excluded_charge / retained_count
        selected += increment
    expected = float(source.sum()) if redistributed else float(source[mask].sum())
    return selected, {
        "retained_atoms": retained_count,
        "excluded_atoms": excluded_count,
        "excluded_charge_e": excluded_charge,
        "equal_increment_per_retained_atom_e": increment,
        "charge_conservation_error_e": float(selected.sum() - expected),
        "redistributed": redistributed,
    }


def write_field(
    path: Path,
    charges: np.ndarray,
    coords: np.ndarray,
    *,
    policy_description: str,
) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(charges)}\n")
        handle.write(f"{policy_description}\n")
        for charge, coord in zip(charges, coords):
            handle.write(f"{charge: .8f} {coord[0]: .8f} {coord[1]: .8f} {coord[2]: .8f}\n")


def severed_bond(
    topology,
    qm_indices: set[int],
    qm_atom_name: str,
    qm_residue,
):
    """Return the unique QM/MM bond cut at a named QM atom."""
    candidates = []
    for bond in topology.bonds():
        atom1, atom2 = bond.atom1, bond.atom2
        if atom1.index in qm_indices and atom2.index not in qm_indices:
            qm_atom, mm_atom = atom1, atom2
        elif atom2.index in qm_indices and atom1.index not in qm_indices:
            qm_atom, mm_atom = atom2, atom1
        else:
            continue
        if qm_atom.residue is qm_residue and qm_atom.name == qm_atom_name:
            candidates.append((qm_atom, mm_atom))
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one severed bond at {qm_residue.name} {qm_residue.id} "
            f"{qm_atom_name}; found {len(candidates)}"
        )
    return candidates[0]


def cap_along_severed_bond(
    qm_atom,
    mm_atom,
    positions: np.ndarray,
    box: np.ndarray,
    length_A: float,
    qm_parent_coord: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Place a link H from the QM parent toward its actual bonded MM partner."""
    vector = minimum_image(positions[mm_atom.index] - positions[qm_atom.index], box)
    partner_distance = float(np.linalg.norm(vector))
    if partner_distance < 1.0e-8:
        raise RuntimeError(f"Zero-length severed bond at atom {qm_atom.index}")
    direction = vector / partner_distance
    cap = qm_parent_coord + direction * length_A
    diagnostic = {
        "qm_parent_index_zero_based": qm_atom.index,
        "qm_parent": f"{qm_atom.residue.name}:{qm_atom.residue.id}:{qm_atom.name}",
        "mm_partner_index_zero_based": mm_atom.index,
        "mm_partner": f"{mm_atom.residue.name}:{mm_atom.residue.id}:{mm_atom.name}",
        "parent_partner_distance_A": partner_distance,
        "cap_distance_A": length_A,
        "cap_direction_dot_parent_to_partner": 1.0,
        "construction": "instantaneous minimum-image severed-bond vector",
    }
    return cap, diagnostic


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--frame-index", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--embedding-cache", type=Path, required=True)
    parser.add_argument(
        "--fixed-embedding-selection",
        type=Path,
        help=(
            "JSON atom selection shared by every frame. If absent at frame 0, "
            "create it from the validated radius/boundary selection; later "
            "frames must reuse exactly those MM atoms."
        ),
    )
    parser.add_argument("--radius", type=float, default=12.0)
    parser.add_argument("--boundary-exclusion", type=float, default=1.8)
    parser.add_argument("--amber-cr2-prmtop", type=Path)
    parser.add_argument(
        "--retain-partner-cr2-charges",
        action="store_true",
        help=(
            "Keep the complete non-QM partner CR2 at every site and assign its "
            "RESP charges from --amber-cr2-prmtop by exact atom name."
        ),
    )
    parser.add_argument(
        "--conserve-boundary-residue-charge",
        action="store_true",
        help=(
            "Redistribute the charge of boundary-excluded non-QM atoms equally "
            "over retained atoms in the same partially selected residue."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    policy = embedding_policy(
        retain_partner_cr2_charges=args.retain_partner_cr2_charges,
        conserve_boundary_residue_charge=(
            args.conserve_boundary_residue_charge
        ),
    )
    policy_hash = json_sha256(policy)
    active_prmtop = (
        args.amber_cr2_prmtop
        if args.retain_partner_cr2_charges
        else None
    )
    if active_prmtop is not None and not active_prmtop.is_file():
        raise FileNotFoundError(f"CR2 AMBER prmtop not found: {active_prmtop}")
    active_prmtop_hash = (
        file_sha256(active_prmtop) if active_prmtop is not None else None
    )
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    pdb = PDBFile(str(args.topology))
    positions, box, dcd_metadata = read_dcd_frame(args.trajectory, args.frame_index)
    atoms = list(pdb.topology.atoms())
    if len(atoms) != len(positions):
        raise RuntimeError(f"Topology/DCD mismatch: {len(atoms)} versus {len(positions)}")
    if not np.isfinite(positions).all():
        raise RuntimeError("Nonfinite coordinates")
    charges = load_or_build_charges(
        pdb,
        args.embedding_cache,
        retain_partner_cr2_charges=args.retain_partner_cr2_charges,
        conserve_boundary_residue_charge=(
            args.conserve_boundary_residue_charge
        ),
        amber_cr2_prmtop=active_prmtop,
    )
    fixed_selection = None
    if args.fixed_embedding_selection is not None and args.fixed_embedding_selection.exists():
        try:
            fixed_selection = json.loads(
                args.fixed_embedding_selection.read_text()
            )
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                "Cannot read fixed MM-selection provenance"
            ) from error
        validate_fixed_embedding_selection(
            fixed_selection,
            topology_descriptor_sha256=topology_descriptor(pdb.topology),
            embedding_policy_payload=policy,
            embedding_policy_sha256=policy_hash,
            amber_cr2_prmtop=active_prmtop,
            amber_cr2_prmtop_sha256=active_prmtop_hash,
            radius_A=args.radius,
            boundary_exclusion_A=args.boundary_exclusion,
            atom_count=len(atoms),
        )
    elif args.fixed_embedding_selection is not None and args.frame_index != 0:
        raise FileNotFoundError(
            "Fixed MM selection must be created from frame 0 before preparing later frames"
        )
    cr2_names, cr2_symbols, reference, tyr_names, tyr_symbols = reference_mapping()
    reference_cr2 = reference[:29]
    reference_heavy = [i for i, symbol in enumerate(cr2_symbols) if symbol != "H"]
    all_cr2_residues = [
        residue for residue in pdb.topology.residues() if residue.name == "CR2"
    ]
    if args.retain_partner_cr2_charges and len(all_cr2_residues) != 2:
        raise RuntimeError(
            "Partner-CR2 embedding requires exactly two CR2 residues; found "
            f"{len(all_cr2_residues)}"
        )

    site_summaries = {}
    for site, keys in SITE_KEYS.items():
        site_dir = args.output_dir / f"site_{site}"
        site_dir.mkdir()
        cr2_residue = residue_by_key(pdb.topology, keys["cr2"], "CR2")
        tyr_residue = residue_by_key(pdb.topology, keys["tyr"], "TYR")
        cr2_by_name = {atom.name: atom for atom in cr2_residue.atoms()}
        if set(cr2_by_name) != set(cr2_names):
            raise RuntimeError(f"Exact CR2 atom inventory mismatch at site {site}")
        cr2_coords = unwrap(
            np.asarray([positions[cr2_by_name[name].index] for name in cr2_names]), box
        )
        rotation, _, fit_rmsd = kabsch(reference_cr2[reference_heavy], cr2_coords[reference_heavy])

        tyr_by_name = {atom.name: atom for atom in tyr_residue.atoms()}
        if any(name not in tyr_by_name for name in tyr_names + ["CB"]):
            raise RuntimeError(f"Stacked Tyr inventory mismatch at site {site}")
        tyr_all_atoms = list(tyr_residue.atoms())
        tyr_all_coords = image_near(
            np.asarray([positions[atom.index] for atom in tyr_all_atoms]), cr2_coords.mean(axis=0), box
        )
        tyr_coord_by_name = dict(zip((atom.name for atom in tyr_all_atoms), tyr_all_coords))
        tyr_coords = np.asarray([tyr_coord_by_name[name] for name in tyr_names])

        cr2_qm_indices = {atom.index for atom in cr2_residue.atoms()}
        phenol_qm_indices = {tyr_by_name[name].index for name in tyr_names}
        qm_indices = cr2_qm_indices | phenol_qm_indices
        n1_cut = severed_bond(pdb.topology, qm_indices, "N1", cr2_residue)
        c3_cut = severed_bond(pdb.topology, qm_indices, "C3", cr2_residue)
        tyr_cut = severed_bond(pdb.topology, qm_indices, "CG", tyr_residue)
        link_coords = []
        link_diagnostics = []
        cap_specs = (
            (n1_cut, 1.01, cr2_coords[cr2_names.index("N1")]),
            (c3_cut, 1.09, cr2_coords[cr2_names.index("C3")]),
            (tyr_cut, 1.09, tyr_coord_by_name["CG"]),
        )
        for cut, length, parent_coord in cap_specs:
            cap, diagnostic = cap_along_severed_bond(
                *cut, positions, box, length, parent_coord
            )
            link_coords.append(cap)
            link_diagnostics.append(diagnostic)
        physical_coords = np.vstack([cr2_coords, tyr_coords])
        qm_coords = np.vstack([physical_coords, np.asarray(link_coords)])
        qm_symbols = cr2_symbols + tyr_symbols + ["H", "H", "H"]
        geometry_path = site_dir / "geometry.xyz"
        with geometry_path.open("w") as handle:
            handle.write("44\n")
            handle.write(f"Production frame {args.frame_index} site {site}; link atoms 42-44\n")
            for symbol, coord in zip(qm_symbols, qm_coords):
                handle.write(f"{symbol:<2s} {coord[0]: .10f} {coord[1]: .10f} {coord[2]: .10f}\n")

        tree = cKDTree(physical_coords)
        field_indices = []
        field_coords = []
        field_charges = []
        component_counts = {"protein": 0, "water": 0, "ions": 0}
        boundary_redistribution = []
        partner_cr2_indices = set()
        retained_indices = (
            None
            if fixed_selection is None
            else set(int(index) for index in fixed_selection["sites"][site]["atom_indices"])
        )
        for residue in pdb.topology.residues():
            non_qm_atoms = [atom for atom in residue.atoms() if atom.index not in qm_indices]
            if not non_qm_atoms:
                continue
            residue_atoms = list(residue.atoms())
            residue_coords = image_near(
                np.asarray([positions[atom.index] for atom in residue_atoms]),
                physical_coords.mean(axis=0),
                box,
            )
            coord_by_index = dict(zip((atom.index for atom in residue_atoms), residue_coords))
            candidate_coords = np.asarray([coord_by_index[atom.index] for atom in non_qm_atoms])
            is_partner_cr2 = bool(
                args.retain_partner_cr2_charges
                and residue.name == "CR2"
                and residue is not cr2_residue
            )
            if retained_indices is None:
                if is_partner_cr2:
                    keep = np.ones(len(non_qm_atoms), dtype=bool)
                else:
                    distances = tree.query(candidate_coords)[0]
                    if float(distances.min()) >= args.radius:
                        continue
                    keep = distances > args.boundary_exclusion
            else:
                keep = np.asarray(
                    [atom.index in retained_indices for atom in non_qm_atoms],
                    dtype=bool,
                )
                if is_partner_cr2 and not bool(keep.all()):
                    raise RuntimeError(
                        f"Fixed MM selection for site {site} omits atoms from "
                        f"partner CR2 {residue.chain.id}:{residue.id}"
                    )
                if not bool(keep.any()):
                    continue
            if (
                args.conserve_boundary_residue_charge
                and not bool(keep.any())
            ):
                raise RuntimeError(
                    "Cannot conserve a boundary residue with no retained MM "
                    f"atoms: {residue.chain.id}:{residue.id}:{residue.name}"
                )
            selected_atoms = [atom for atom, take in zip(non_qm_atoms, keep) if take]
            selected_coords = candidate_coords[keep]
            selected_charges, charge_diagnostic = (
                redistribute_residue_boundary_charge(
                    charges[
                        np.asarray(
                            [atom.index for atom in non_qm_atoms], dtype=int
                        )
                    ],
                    keep,
                    enabled=args.conserve_boundary_residue_charge,
                )
            )
            field_indices.extend(atom.index for atom in selected_atoms)
            field_coords.extend(selected_coords)
            field_charges.extend(selected_charges)
            if is_partner_cr2:
                partner_cr2_indices.update(atom.index for atom in selected_atoms)
            if 0 < int(keep.sum()) < len(non_qm_atoms):
                boundary_redistribution.append(
                    {
                        "residue": (
                            f"{residue.chain.id}:{residue.id}:{residue.name}"
                        ),
                        **charge_diagnostic,
                    }
                )
            if residue.name in {"HOH", "WAT", "TIP3", "TIP3P"}:
                component = "water"
            elif residue.name in {"NA", "CL", "K", "MG", "CA", "ZN"}:
                component = "ions"
            else:
                component = "protein"
            component_counts[component] += len(selected_atoms)
        if retained_indices is not None and set(field_indices) != retained_indices:
            missing = sorted(retained_indices - set(field_indices))
            extra = sorted(set(field_indices) - retained_indices)
            raise RuntimeError(
                f"Fixed MM atom selection mismatch at site {site}: "
                f"missing={missing[:10]} extra={extra[:10]}"
            )
        field_coords_array = np.asarray(field_coords, dtype=float)
        field_charges_array = np.asarray(field_charges, dtype=float)
        if args.retain_partner_cr2_charges:
            expected_partner_indices = {
                atom.index
                for residue in all_cr2_residues
                if residue is not cr2_residue
                for atom in residue.atoms()
            }
            if partner_cr2_indices != expected_partner_indices:
                raise RuntimeError(
                    f"Site {site} did not retain the complete partner CR2"
                )
        mm_charges_path = site_dir / "mm_charges.dat"
        write_field(
            mm_charges_path,
            field_charges_array,
            field_coords_array,
            policy_description=(
                f"{args.radius:g} A AMBER14/TIP3P embedding; complete "
                "partner CR2 RESP; partial-residue boundary charge conserved"
                if (
                    args.retain_partner_cr2_charges
                    and args.conserve_boundary_residue_charge
                )
                else (
                    f"{args.radius:g} A AMBER14/TIP3P electrostatic embedding"
                )
            ),
        )
        shutil.copy2(
            Path(__file__).resolve().parents[1] / "input_templates/tddft_energy.in",
            site_dir / "tddft.in",
        )
        site_summaries[site] = {
            "cr2_key": list(keys["cr2"]),
            "stacked_tyr_key": list(keys["tyr"]),
            "cr2_reference_heavy_fit_rmsd_A": fit_rmsd,
            "qm_atoms": 44,
            "physical_qm_atoms": 41,
            "excluded_link_atom_indices_one_based": [42, 43, 44],
            "link_caps": link_diagnostics,
            "mm_point_charges": len(field_charges_array),
            "mm_net_charge_e": float(field_charges_array.sum()),
            "mm_component_atom_counts": component_counts,
            "partner_cr2_retained_atom_count": len(partner_cr2_indices),
            "partner_cr2_net_charge_e": float(
                charges[np.asarray(sorted(partner_cr2_indices), dtype=int)].sum()
            ),
            "boundary_charge_redistribution": boundary_redistribution,
            "maximum_boundary_charge_conservation_error_e": float(
                max(
                    (
                        abs(item["charge_conservation_error_e"])
                        for item in boundary_redistribution
                    ),
                    default=0.0,
                )
            ),
            "mm_atom_indices": field_indices,
            "geometry_xyz_sha256": file_sha256(geometry_path),
            "mm_charges_dat_sha256": file_sha256(mm_charges_path),
        }

    if args.fixed_embedding_selection is not None and fixed_selection is None:
        selection_payload = {
            "status": "fixed_embedding_selection",
            "source_frame_index": args.frame_index,
            "topology": str(args.topology.resolve()),
            "topology_descriptor_sha256": topology_descriptor(pdb.topology),
            "embedding_policy": policy,
            "embedding_policy_sha256": policy_hash,
            "amber_cr2_prmtop": (
                str(active_prmtop.resolve())
                if active_prmtop is not None
                else None
            ),
            "amber_cr2_prmtop_sha256": active_prmtop_hash,
            "radius_A_at_selection": args.radius,
            "boundary_exclusion_A_at_selection": args.boundary_exclusion,
            "policy": (
                "The frame-0 validated local QM/MM atom selection is held fixed "
                "for all later frames; coordinates are still instantaneously imaged."
            ),
            "sites": {
                site: {
                    "atom_indices": site_summaries[site]["mm_atom_indices"],
                    "atom_count": site_summaries[site]["mm_point_charges"],
                    "net_charge_e": site_summaries[site]["mm_net_charge_e"],
                }
                for site in ("A", "B")
            },
        }
        atomic_write_json(args.fixed_embedding_selection, selection_payload)

    for site in ("A", "B"):
        # The atom list is retained in the shared provenance file, not repeated
        # in every frame-level summary.
        del site_summaries[site]["mm_atom_indices"]

    canonical = np.asarray(positions, dtype="<f4").tobytes(order="C")
    summary = {
        "status": "prepared_production_frame",
        "topology": str(args.topology.resolve()),
        "topology_file_sha256": hashlib.sha256(
            args.topology.read_bytes()
        ).hexdigest(),
        "trajectory": str(args.trajectory.resolve()),
        "frame_index": args.frame_index,
        "coordinate_sha256": hashlib.sha256(canonical).hexdigest(),
        "topology_descriptor_sha256": topology_descriptor(pdb.topology),
        "dcd": dcd_metadata,
        "box_A": box.tolist(),
        "embedding_radius_A": args.radius,
        "embedding_boundary_exclusion_A": args.boundary_exclusion,
        "embedding_policy": policy,
        "embedding_policy_sha256": policy_hash,
        "retain_partner_cr2_charges": bool(
            args.retain_partner_cr2_charges
        ),
        "conserve_boundary_residue_charge": bool(
            args.conserve_boundary_residue_charge
        ),
        "amber_cr2_prmtop": (
            str(active_prmtop.resolve())
            if active_prmtop is not None
            else None
        ),
        "amber_cr2_prmtop_sha256": active_prmtop_hash,
        "embedding_charges": str(args.embedding_cache.resolve()),
        "embedding_charges_sha256": file_sha256(args.embedding_cache),
        "embedding_selection_mode": (
            "fixed_frame0_atom_indices"
            if args.fixed_embedding_selection is not None
            else "instantaneous_radius_selection"
        ),
        "fixed_embedding_selection": (
            str(args.fixed_embedding_selection.resolve())
            if args.fixed_embedding_selection is not None
            else None
        ),
        "fixed_embedding_selection_sha256": (
            file_sha256(args.fixed_embedding_selection)
            if args.fixed_embedding_selection is not None
            else None
        ),
        "preparation_script_sha256": file_sha256(Path(__file__).resolve()),
        "forcefield": (
            "amber14-all.xml + amber14/tip3p.xml MM charges; partner CR2 "
            "RESP charges from configured AMBER prmtop"
            if args.retain_partner_cr2_charges
            else (
                "amber14-all.xml + amber14/tip3p.xml MM charges; "
                "CR2 wholly QM"
            )
        ),
        "sites": site_summaries,
    }
    atomic_write_json(args.output_dir / "preparation.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Prepare both instantaneous 44-atom QM/MM site jobs from one production DCD frame."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct

import numpy as np
from scipy.spatial import cKDTree
from openmm import NonbondedForce, unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile

from prepare_surrogate_pilot import (
    PHENOL,
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


def load_or_build_charges(pdb: PDBFile, cache: Path) -> np.ndarray:
    descriptor = topology_descriptor(pdb.topology)
    if cache.exists():
        saved = np.load(cache)
        if str(saved["topology_descriptor_sha256"]) != descriptor:
            raise RuntimeError("Embedding-cache topology descriptor mismatch")
        charges = np.asarray(saved["charges_e"], dtype=float)
        if len(charges) != pdb.topology.getNumAtoms():
            raise RuntimeError("Embedding-cache atom-count mismatch")
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
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        charges_e=charges,
        topology_descriptor_sha256=np.asarray(descriptor),
        forcefield=np.asarray("amber14-all.xml + amber14/tip3p.xml; CR2 wholly QM"),
    )
    return charges


def write_field(path: Path, charges: np.ndarray, coords: np.ndarray) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(charges)}\n")
        handle.write("Residue-complete 12 A AMBER14/TIP3P electrostatic embedding\n")
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
    parser.add_argument(
        "--site-b-offset",
        type=int,
        default=0,
        help=(
            "Residue-number offset to add to SITE B's CR2/Tyr keys. The corrected "
            "507-residue construct restores the 8-residue LGMDELYK tail at the "
            "barrel-1 C-terminus, which shifts everything downstream: CR2 328->336 "
            "and Tyr 464->472. Site A (66/202) precedes the insertion and is "
            "unchanged. Default 0 reproduces the original 495-aa numbering exactly."
        ),
    )
    parser.add_argument(
        "--union-embedding-radius",
        type=float,
        default=None,
        help=(
            "Select whole residues within this radius of EITHER site's QM region "
            "and give BOTH sites that identical atom set. For a C2 dimer the union "
            "is itself C2-symmetric, so the two sites see equivalent environments "
            "and the A/B asymmetry that broke the per-site spheres cannot arise. "
            "Much cheaper than --full-system-embedding. Overrides --radius."
        ),
    )
    parser.add_argument(
        "--full-system-embedding",
        action="store_true",
        help=(
            "Retain EVERY non-QM atom in the MM field instead of a radius shell. "
            "Both sites then see the identical whole system, so their MM "
            "complements carry the same net charge by construction. Fixes the "
            "pathology diagnosed 2026-08-28: the residue-complete 12 A cutoff is "
            "charge-arbitrary -- between two structures of the same protein, "
            "Glu141/Glu403 left site A's sphere while Lys106/Arg429/Glu89 "
            "entered it, a +3 e swing worth +2550 cm^-1 on the site energy. "
            "Overrides --radius."
        ),
    )
    parser.add_argument(
        "--link-only-exclusion",
        action="store_true",
        help=(
            "Drop ONLY the MM atoms covalently bonded across the QM boundary "
            "(the three severed-bond partners), instead of every MM atom within "
            "--boundary-exclusion of the QM region. The geometric cut also "
            "deletes hydrogen-bond partners whose membership flips with thermal "
            "H-bond length, injecting site- and frame-dependent noise."
        ),
    )
    parser.add_argument("--radius", type=float, default=12.0)
    parser.add_argument("--boundary-exclusion", type=float, default=1.8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
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
    charges = load_or_build_charges(pdb, args.embedding_cache)
    fixed_selection = None
    if args.fixed_embedding_selection is not None and args.fixed_embedding_selection.exists():
        fixed_selection = json.loads(args.fixed_embedding_selection.read_text())
        if fixed_selection["topology_descriptor_sha256"] != topology_descriptor(pdb.topology):
            raise RuntimeError("Fixed MM selection topology does not match this frame")
    elif args.fixed_embedding_selection is not None and args.frame_index != 0:
        raise FileNotFoundError(
            "Fixed MM selection must be created from frame 0 before preparing later frames"
        )
    cr2_names, cr2_symbols, reference, tyr_names, tyr_symbols = reference_mapping()
    reference_cr2 = reference[:29]
    reference_heavy = [i for i, symbol in enumerate(cr2_symbols) if symbol != "H"]

    site_keys = {
        site: {
            role: (chain, str(int(num) + (args.site_b_offset if site == "B" else 0)))
            for role, (chain, num) in keys.items()
        }
        for site, keys in SITE_KEYS.items()
    }
    if args.site_b_offset:
        print(f"[site-keys] site B offset {args.site_b_offset:+d}: {site_keys['B']}")

    union_indices = None
    if args.union_embedding_radius is not None:
        all_qm = set()
        for site, keys in site_keys.items():
            cr2_r = residue_by_key(pdb.topology, keys["cr2"], "CR2")
            tyr_r = residue_by_key(pdb.topology, keys["tyr"], "TYR")
            all_qm |= {a.index for a in cr2_r.atoms()}
            all_qm |= {a.index for a in tyr_r.atoms()}
        anchor = np.asarray([positions[i] for i in sorted(all_qm)])
        atree = cKDTree(anchor)
        union_indices = set()
        for residue in pdb.topology.residues():
            ridx = [a.index for a in residue.atoms()]
            rc = image_near(np.asarray([positions[i] for i in ridx]),
                            anchor.mean(axis=0), box)
            if float(atree.query(rc)[0].min()) < args.union_embedding_radius:
                union_indices.update(ridx)
        print(f"[union] radius {args.union_embedding_radius} A -> "
              f"{len(union_indices)} atoms shared by both sites")

    site_summaries = {}
    for site, keys in site_keys.items():
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
        with (site_dir / "geometry.xyz").open("w") as handle:
            handle.write("44\n")
            handle.write(f"Production frame {args.frame_index} site {site}; link atoms 42-44\n")
            for symbol, coord in zip(qm_symbols, qm_coords):
                handle.write(f"{symbol:<2s} {coord[0]: .10f} {coord[1]: .10f} {coord[2]: .10f}\n")

        link_partner_indices = {n1_cut[1].index, c3_cut[1].index, tyr_cut[1].index}
        tree = cKDTree(physical_coords)
        field_indices = []
        field_coords = []
        component_counts = {"protein": 0, "water": 0, "ions": 0}
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
            if retained_indices is None and union_indices is not None:
                keep = np.asarray(
                    [(atom.index in union_indices)
                     and (atom.index not in link_partner_indices)
                     for atom in non_qm_atoms],
                    dtype=bool,
                )
            elif retained_indices is None:
                if not args.full_system_embedding:
                    if float(tree.query(candidate_coords)[0].min()) >= args.radius:
                        continue
                if args.link_only_exclusion:
                    # topological boundary: only the severed-bond partners
                    keep = np.asarray(
                        [atom.index not in link_partner_indices for atom in non_qm_atoms],
                        dtype=bool,
                    )
                else:
                    keep = tree.query(candidate_coords)[0] > args.boundary_exclusion
            else:
                keep = np.asarray(
                    [atom.index in retained_indices for atom in non_qm_atoms],
                    dtype=bool,
                )
            selected_atoms = [atom for atom, take in zip(non_qm_atoms, keep) if take]
            selected_coords = candidate_coords[keep]
            field_indices.extend(atom.index for atom in selected_atoms)
            field_coords.extend(selected_coords)
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
        field_charges = charges[np.asarray(field_indices, dtype=int)]
        write_field(site_dir / "mm_charges.dat", field_charges, field_coords_array)
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
            "mm_point_charges": len(field_charges),
            "mm_net_charge_e": float(field_charges.sum()),
            "mm_component_atom_counts": component_counts,
            "mm_atom_indices": field_indices,
        }

    if args.fixed_embedding_selection is not None and fixed_selection is None:
        selection_payload = {
            "status": "fixed_embedding_selection",
            "source_frame_index": args.frame_index,
            "topology": str(args.topology.resolve()),
            "topology_descriptor_sha256": topology_descriptor(pdb.topology),
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
        args.fixed_embedding_selection.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.fixed_embedding_selection.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(selection_payload, indent=2) + "\n")
        temporary.replace(args.fixed_embedding_selection)

    for site in ("A", "B"):
        # The atom list is retained in the shared provenance file, not repeated
        # in every frame-level summary.
        del site_summaries[site]["mm_atom_indices"]

    canonical = np.asarray(positions, dtype="<f4").tobytes(order="C")
    summary = {
        "status": "prepared_production_frame",
        "topology": str(args.topology.resolve()),
        "trajectory": str(args.trajectory.resolve()),
        "frame_index": args.frame_index,
        "coordinate_sha256": hashlib.sha256(canonical).hexdigest(),
        "topology_descriptor_sha256": topology_descriptor(pdb.topology),
        "dcd": dcd_metadata,
        "box_A": box.tolist(),
        "embedding_radius_A": args.radius,
        "embedding_boundary_exclusion_A": args.boundary_exclusion,
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
        "forcefield": "amber14-all.xml + amber14/tip3p.xml MM charges; CR2 wholly QM",
        "sites": site_summaries,
    }
    (args.output_dir / "preparation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

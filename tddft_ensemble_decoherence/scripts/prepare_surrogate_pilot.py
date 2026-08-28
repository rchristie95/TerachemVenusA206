#!/usr/bin/env python3
"""Prepare a labelled three-frame surrogate TDDFT pilot from the local old NVT PDB.

This does not convert the old trajectory into the missing production ensemble.
It preserves instantaneous CR2 heavy atoms and the stacked Tyr phenol, places
the validated production hydrogens/link atoms, and rebuilds the original 12 A
residue-complete AMBER/TIP3P electrostatic field for a feasibility pilot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy.spatial import cKDTree
from openmm import NonbondedForce, unit
from openmm.app import AmberPrmtopFile, ForceField, Modeller, NoCutoff, PDBFile

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
import run_nvt  # noqa: E402

PHENOL = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HD1", "HD2", "HE1", "HE2", "HH"}
SITE_KEYS = {
    "A": {"cr2": ("B", "1"), "tyr": ("C", "136")},
    "B": {"cr2": ("D", "1"), "tyr": ("E", "25")},
}


def parse_residue_key(value: str) -> tuple[str, str]:
    parts = value.split(":", 1)
    if len(parts) != 2 or not all(parts):
        raise argparse.ArgumentTypeError("Residue keys must be CHAIN:RESID")
    return parts[0], parts[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target_center - source_center @ rotation.T
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def unwrap(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    return coords[0] + minimum_image(coords - coords[0], box)


def image_near(coords: np.ndarray, center: np.ndarray, box: np.ndarray) -> np.ndarray:
    coords = unwrap(coords, box)
    shift = minimum_image(coords.mean(axis=0) - center, box) - (coords.mean(axis=0) - center)
    return coords + shift


def reference_mapping() -> tuple[list[str], list[str], np.ndarray, list[str], list[str]]:
    prm = AmberPrmtopFile(str(REPO / "anionic_build/monomer_solv.prmtop"))
    cr2 = next(residue for residue in prm.topology.residues() if residue.name == "CR2")
    tyr = next(residue for residue in prm.topology.residues() if residue.name == "TYR" and residue.id == "202")
    cr2_atoms = list(cr2.atoms())
    tyr_atoms = [atom for atom in tyr.atoms() if atom.name in PHENOL]
    lines = (REPO / "tc_tddft_44/geometry.xyz").read_text().splitlines()[2:46]
    symbols = [line.split()[0] for line in lines]
    coords = np.asarray([[float(value) for value in line.split()[1:4]] for line in lines])
    if len(cr2_atoms) != 29 or len(tyr_atoms) != 12 or len(coords) != 44:
        raise RuntimeError("Validated reference mapping is not 29 CR2 + 12 Tyr + 3 links")
    return (
        [atom.name for atom in cr2_atoms],
        [atom.element.symbol for atom in cr2_atoms],
        coords,
        [atom.name for atom in tyr_atoms],
        [atom.element.symbol for atom in tyr_atoms],
    )


def residue_by_key(topology, key: tuple[str, str], name: str):
    matches = [
        residue for residue in topology.residues()
        if residue.chain.id == key[0] and residue.id == key[1] and residue.name == name
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {name} residue at {key}; found {len(matches)}")
    return matches[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frame_pdb", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--radius", type=float, default=12.0)
    parser.add_argument("--site-a-cr2", type=parse_residue_key, default=SITE_KEYS["A"]["cr2"])
    parser.add_argument("--site-a-tyr", type=parse_residue_key, default=SITE_KEYS["A"]["tyr"])
    parser.add_argument("--site-b-cr2", type=parse_residue_key, default=SITE_KEYS["B"]["cr2"])
    parser.add_argument("--site-b-tyr", type=parse_residue_key, default=SITE_KEYS["B"]["tyr"])
    parser.add_argument(
        "--rigid-reference-cr2",
        action="store_true",
        help="Use the validated rigid CR2 after recording a failed old-trajectory geometry fit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    pdb = PDBFile(str(args.frame_pdb))
    if pdb.topology.getNumAtoms() != 69418:
        raise RuntimeError(f"Surrogate topology must have 69418 atoms; got {pdb.topology.getNumAtoms()}")
    cr2_bonds_before = sum(
        1 for bond in pdb.topology.bonds()
        if bond.atom1.residue.name == "CR2" or bond.atom2.residue.name == "CR2"
    )
    inferred = run_nvt.infer_nonstandard_residue_bonds(pdb.topology, pdb.positions, ["CR2"])
    cr2_bonds_after = sum(
        1 for bond in pdb.topology.bonds()
        if bond.atom1.residue.name == "CR2" or bond.atom2.residue.name == "CR2"
    )
    if not args.rigid_reference_cr2 and cr2_bonds_after < 60:
        raise RuntimeError(
            f"Only {cr2_bonds_after} CR2 bonds remain across two sites; "
            "refusing instantaneous geometry from a disconnected topology"
        )
    # CR2 is wholly QM and its obsolete placeholder template is intentionally
    # excluded from MM charge construction.  Deleting both CR2 residues also
    # prevents a disconnected or newly repaired CR2 bond graph from being
    # forced through that incompatible template.
    original_atoms = list(pdb.topology.atoms())
    retained_atoms = [atom for atom in original_atoms if atom.residue.name != "CR2"]
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.delete([atom for atom in original_atoms if atom.residue.name == "CR2"])
    modelled_atoms = list(modeller.topology.atoms())
    if len(modelled_atoms) != len(retained_atoms):
        raise RuntimeError("CR2 deletion changed the retained MM atom inventory")
    for old, new in zip(retained_atoms, modelled_atoms):
        if (old.name, old.residue.name, old.residue.id, old.residue.chain.id) != (
            new.name, new.residue.name, new.residue.id, new.residue.chain.id
        ):
            raise RuntimeError("MM atom order changed during CR2 deletion")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=NoCutoff, ignoreExternalBonds=True
    )
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    retained_charges = np.asarray([
        nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
        for index in range(system.getNumParticles())
    ])
    charges = np.zeros(len(original_atoms), dtype=float)
    charges[[atom.index for atom in retained_atoms]] = retained_charges
    positions = np.asarray(pdb.getPositions().value_in_unit(unit.angstrom))
    vectors = pdb.topology.getPeriodicBoxVectors()
    if vectors is None:
        raise RuntimeError("Surrogate frame has no periodic box")
    box_matrix = np.asarray([[component.value_in_unit(unit.angstrom) for component in vector] for vector in vectors])
    if not np.allclose(box_matrix, np.diag(np.diag(box_matrix)), atol=1.0e-6):
        raise RuntimeError("This pilot only supports the recorded orthorhombic surrogate box")
    box = np.diag(box_matrix)
    atoms = list(pdb.topology.atoms())
    cr2_names, cr2_symbols, reference, tyr_names, tyr_symbols = reference_mapping()
    reference_cr2 = reference[:29]
    reference_links = reference[41:44]
    reference_heavy = [index for index, symbol in enumerate(cr2_symbols) if symbol != "H"]

    site_keys = {
        "A": {"cr2": args.site_a_cr2, "tyr": args.site_a_tyr},
        "B": {"cr2": args.site_b_cr2, "tyr": args.site_b_tyr},
    }
    all_site_summaries = {}
    for site, keys in site_keys.items():
        site_dir = args.output_dir / f"site_{site}"
        site_dir.mkdir()
        cr2_residue = residue_by_key(pdb.topology, keys["cr2"], "CR2")
        tyr_residue = residue_by_key(pdb.topology, keys["tyr"], "TYR")
        old_cr2 = {atom.name: atom for atom in cr2_residue.atoms()}
        target_heavy_names = [cr2_names[index] for index in reference_heavy]
        if any(name not in old_cr2 for name in target_heavy_names):
            raise RuntimeError(f"Surrogate CR2 lacks shared production heavy atoms at site {site}")
        target_heavy = unwrap(np.asarray([positions[old_cr2[name].index] for name in target_heavy_names]), box)
        rotation, translation, fit_rmsd = kabsch(reference_cr2[reference_heavy], target_heavy)
        if fit_rmsd > 0.8 and not args.rigid_reference_cr2:
            raise RuntimeError(
                f"Site {site} CR2 fit RMSD {fit_rmsd:.3f} A exceeds the 0.8 A "
                "surrogate safety limit; refusing to use malformed instantaneous geometry"
            )
        cr2_coords = reference_cr2 @ rotation.T + translation
        if not args.rigid_reference_cr2:
            for name, coord in zip(target_heavy_names, target_heavy):
                cr2_coords[cr2_names.index(name)] = coord

        tyr_by_name = {atom.name: atom for atom in tyr_residue.atoms()}
        if any(name not in tyr_by_name for name in tyr_names + ["CB"]):
            raise RuntimeError(f"Surrogate Tyr lacks validated phenol/CB atoms at site {site}")
        tyr_all_atoms = list(tyr_residue.atoms())
        tyr_all_coords = image_near(
            np.asarray([positions[atom.index] for atom in tyr_all_atoms]), cr2_coords.mean(axis=0), box
        )
        tyr_coord_by_name = {atom.name: coord for atom, coord in zip(tyr_all_atoms, tyr_all_coords)}
        tyr_coords = np.asarray([tyr_coord_by_name[name] for name in tyr_names])
        tyr_placement_rmsd = None
        if args.rigid_reference_cr2:
            reference_tyr = reference[29:41]
            tyr_heavy = [index for index, symbol in enumerate(tyr_symbols) if symbol != "H"]
            rotation, translation, tyr_placement_rmsd = kabsch(
                reference_tyr[tyr_heavy], tyr_coords[tyr_heavy]
            )
            cr2_coords = reference_cr2 @ rotation.T + translation

        link_coords = []
        link_symbols = ["H", "H", "H"]
        for link_index, parent_index in enumerate((0, 9)):
            direction = reference_links[link_index] - reference[parent_index]
            direction = direction @ rotation.T
            length = 1.01 if parent_index == 0 else 1.09
            link_coords.append(cr2_coords[parent_index] + direction / np.linalg.norm(direction) * length)
        cg = tyr_coord_by_name["CG"]
        cb = tyr_coord_by_name["CB"]
        direction = cg - cb
        link_coords.append(cg + direction / np.linalg.norm(direction) * 1.09)
        physical_coords = np.vstack([cr2_coords, tyr_coords])
        qm_coords = np.vstack([physical_coords, np.asarray(link_coords)])
        qm_symbols = cr2_symbols + tyr_symbols + link_symbols

        with (site_dir / "geometry.xyz").open("w") as handle:
            handle.write("44\n")
            handle.write(f"SURROGATE old-NVT pilot site {site}; not production ensemble\n")
            for symbol, coord in zip(qm_symbols, qm_coords):
                handle.write(f"{symbol:<2s} {coord[0]: .10f} {coord[1]: .10f} {coord[2]: .10f}\n")

        old_qm_indices = {atom.index for atom in cr2_residue.atoms()}
        old_qm_indices.update(tyr_by_name[name].index for name in tyr_names)
        tree = cKDTree(physical_coords)
        field_indices: list[int] = []
        field_coords: list[np.ndarray] = []
        field_components: list[str] = []
        for residue in pdb.topology.residues():
            non_qm_atoms = [atom for atom in residue.atoms() if atom.index not in old_qm_indices]
            if not non_qm_atoms:
                continue
            residue_atoms = list(residue.atoms())
            residue_coords = image_near(
                np.asarray([positions[atom.index] for atom in residue_atoms]), physical_coords.mean(axis=0), box
            )
            coord_by_index = {atom.index: coord for atom, coord in zip(residue_atoms, residue_coords)}
            non_qm_coords = np.asarray([coord_by_index[atom.index] for atom in non_qm_atoms])
            if float(tree.query(non_qm_coords)[0].min()) < args.radius:
                field_indices.extend(atom.index for atom in non_qm_atoms)
                field_coords.extend(non_qm_coords)
                if residue.name in run_nvt.WATER_RESIDUE_NAMES:
                    component = "water"
                elif residue.name in run_nvt.COMMON_ION_RESIDUES:
                    component = "ions"
                else:
                    component = "protein"
                field_components.extend([component] * len(non_qm_atoms))
        field_coords_array = np.asarray(field_coords)
        field_charges = charges[field_indices]
        field_components_array = np.asarray(field_components)
        keep = tree.query(field_coords_array)[0] > 1.8
        field_coords_array = field_coords_array[keep]
        field_charges = field_charges[keep]
        field_components_array = field_components_array[keep]

        def write_field(path: Path, selected: np.ndarray, label: str) -> None:
            with path.open("w") as handle:
                handle.write(f"{int(selected.sum())}\n")
                handle.write(f"SURROGATE {label}; residue-complete 12 A AMBER14/TIP3P field; QM CR2 excluded\n")
                for charge, coord in zip(field_charges[selected], field_coords_array[selected]):
                    handle.write(f"{charge: .8f} {coord[0]: .8f} {coord[1]: .8f} {coord[2]: .8f}\n")

        write_field(site_dir / "mm_charges.dat", np.ones(len(field_charges), dtype=bool), "full MM")
        for component in ("protein", "water", "ions"):
            write_field(
                site_dir / f"mm_charges_{component}.dat",
                field_components_array == component,
                f"{component} component",
            )
        shutil.copy2(REPO / "tddft_ensemble_decoherence/input_templates/tddft_energy.in", site_dir / "tddft.in")
        all_site_summaries[site] = {
            "surrogate_cr2_key": list(keys["cr2"]),
            "surrogate_stacked_tyr_key": list(keys["tyr"]),
            "shared_heavy_fit_rmsd_A": fit_rmsd,
            "reference_tyr_placement_fit_rmsd_A": tyr_placement_rmsd,
            "mm_point_charges": int(len(field_charges)),
            "mm_net_charge_e": float(field_charges.sum()),
            "mm_components": {
                component: {
                    "count": int(np.count_nonzero(field_components_array == component)),
                    "net_charge_e": float(field_charges[field_components_array == component].sum()),
                }
                for component in ("protein", "water", "ions")
            },
            "qm_atoms": 44,
            "physical_qm_atoms": 41,
        }

    summary = {
        "status": "prepared_rigid_reference_environment_surrogate" if args.rigid_reference_cr2 else "prepared_valid_initial_surrogate",
        "warning": (
            "The local 69418-atom/31-atom-CR2 trajectory predates the exact 69471-atom/29-atom-CR2 production windows and its CR2 heavy-atom fits fail the production geometry criterion. Results are electrostatic-environment feasibility diagnostics only."
            if args.rigid_reference_cr2 else
            "This geometrically valid pre-minimization 69418-atom/31-atom-CR2 structure is not the requested 69471-atom/29-atom production ensemble. Results are initial-structure diagnostics only."
        ),
        "source": str(args.frame_pdb.resolve()),
        "source_sha256": sha256(args.frame_pdb),
        "box_A": box.tolist(),
        "forcefield": "amber14-all.xml + amber14/tip3p.xml on MM atoms; both wholly-QM CR2 residues deleted before charge construction",
        "cr2_handling": (
            "validated reference CR2+Tyr stacking geometry placed by the instantaneous stacked-Tyr heavy atoms because old CR2 geometry failed preflight"
            if args.rigid_reference_cr2 else
            "19 instantaneous heavy atoms retained; validated production H/link positions placed from the same reference mapping"
        ),
        "embedding_radius_A": args.radius,
        "embedding_boundary_exclusion_A": 1.8,
        "inferred_cr2_bonds": inferred,
        "cr2_bonds_before_inference": cr2_bonds_before,
        "cr2_bonds_after_inference": cr2_bonds_after,
        "sites": all_site_summaries,
    }
    (args.output_dir / "preparation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

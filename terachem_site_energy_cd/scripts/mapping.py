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


__all__ = ["image_near", "kabsch", "minimum_image", "reference_mapping", "residue_by_key", "unwrap"]

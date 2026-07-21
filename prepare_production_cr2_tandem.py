#!/usr/bin/env python3
"""Replace both tandem CR2 residues with the retained production AMBER CR2.

The tandem crystal model contains only the 19 heavy CR2 atoms.  The retained
AMBER monomer topology instead uses the 29-atom, anionic CR2 residue.  This
script rigidly fits those 19 shared heavy atoms at each tandem site and writes
all 29 production atoms, including the production hydrogen names and geometry.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def atom_name(line: str) -> str:
    return line[12:16].strip()


def residue_key(line: str) -> tuple[str, str, str]:
    return line[21:22], line[22:26], line[26:27]


def xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


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


def format_atom(
    serial: int,
    name: str,
    chain: str,
    resseq: str,
    icode: str,
    position: np.ndarray,
    element: str,
) -> str:
    atom_field = f" {name:<3}" if len(name) < 4 else f"{name:>4}"
    return (
        f"HETATM{serial:5d} {atom_field} CR2 {chain:1}{int(resseq):4d}{icode:1}   "
        f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}"
        f"  1.00  0.00          {element:>2}  "
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tandem", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tandem_lines = args.tandem.read_text(errors="replace").splitlines()
    reference_lines = [
        line
        for line in args.reference.read_text(errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CR2"
    ]
    if len(reference_lines) != 29:
        raise RuntimeError(f"Expected 29 production CR2 atoms, found {len(reference_lines)}")

    reference_by_name = {atom_name(line): line for line in reference_lines}
    if len(reference_by_name) != len(reference_lines):
        raise RuntimeError("Production CR2 atom names are not unique")

    target_groups: dict[tuple[str, str, str], list[str]] = {}
    for line in tandem_lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CR2":
            target_groups.setdefault(residue_key(line), []).append(line)
    if len(target_groups) != 2:
        raise RuntimeError(f"Expected two tandem CR2 residues, found {len(target_groups)}")

    replacements: dict[tuple[str, str, str], list[tuple[str, np.ndarray, str]]] = {}
    heavy_rmsds = []
    reference_heavy_names = [
        name
        for name, line in reference_by_name.items()
        if line[76:78].strip().upper() != "H"
    ]
    for key, target_group in target_groups.items():
        target_by_name = {atom_name(line): line for line in target_group}
        shared = sorted(set(reference_heavy_names) & set(target_by_name))
        if len(shared) != 19:
            raise RuntimeError(f"Expected 19 shared CR2 heavy atoms at {key}, found {len(shared)}")
        reference_xyz = np.array([xyz(reference_by_name[name]) for name in shared])
        target_xyz = np.array([xyz(target_by_name[name]) for name in shared])
        rotation, translation, rmsd = kabsch(reference_xyz, target_xyz)
        heavy_rmsds.append(rmsd)

        transformed = []
        for line in reference_lines:
            name = atom_name(line)
            position = xyz(line) @ rotation.T + translation
            element = line[76:78].strip() or name[0]
            transformed.append((name, position, element))
        replacements[key] = transformed

    output_lines = []
    emitted = set()
    serial = 1
    for line in tandem_lines:
        is_atom = line.startswith(("ATOM  ", "HETATM"))
        is_cr2 = is_atom and line[17:20].strip() == "CR2"
        if is_cr2:
            key = residue_key(line)
            if key in emitted:
                continue
            emitted.add(key)
            chain, resseq, icode = key
            for name, position, element in replacements[key]:
                output_lines.append(format_atom(serial, name, chain, resseq, icode, position, element))
                serial += 1
            continue
        if is_atom:
            output_lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        elif not line.startswith("CONECT"):
            output_lines.append(line)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(output_lines) + "\n")
    print(f"[CR2] wrote {args.output}")
    print(f"[CR2] sites={len(replacements)}, atoms/site=29, shared-heavy/site=19")
    print("[CR2] rigid-fit RMSD (A): " + ", ".join(f"{value:.6f}" for value in heavy_rmsds))


if __name__ == "__main__":
    main()

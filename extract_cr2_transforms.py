#!/usr/bin/env python3
"""Extract per-frame rigid CR2 fits for a binary tandem trajectory.

This is run in the PyMOL environment because PyMOL can read the DCD directly.
The resulting small NPZ cache is then consumed by coupling_dcd_steom.py.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def cr2_atoms(pdb_path: Path) -> dict[str, np.ndarray]:
    atoms = {}
    with open(pdb_path) as handle:
        for line in handle:
            if line[:6].strip() in ("ATOM", "HETATM") and line[17:20].strip() == "CR2":
                element = line[76:78].strip().upper()
                if element == "H" or (not element and line[12:16].strip().upper().startswith("H")):
                    continue
                atoms[line[12:16].strip()] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float
                )
    return atoms


def kabsch(source: np.ndarray, target: np.ndarray):
    source_c = source - source.mean(axis=0)
    target_c = target - target.mean(axis=0)
    u, _, vt = np.linalg.svd(source_c.T @ target_c)
    handedness = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, handedness]) @ u.T
    translation = target.mean(axis=0) - rotation @ source.mean(axis=0)
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traj", type=Path, required=True)
    ap.add_argument("--topology", type=Path, required=True)
    ap.add_argument("--monomer", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    try:
        from pymol import cmd
        cmd.get_version()
    except Exception:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd

    cmd.reinitialize()
    cmd.load(str(args.topology), "traj")
    cmd.load_traj(str(args.traj), "traj", state=1)
    n_states = cmd.count_states("traj")

    first = cmd.get_model("traj and resn CR2", state=1)
    residue_keys = []
    for atom in first.atom:
        key = ((atom.chain or "").strip(), str(atom.resi).strip())
        if key not in residue_keys:
            residue_keys.append(key)
    if len(residue_keys) != 2:
        raise RuntimeError(f"Expected two CR2 residues, found {residue_keys}")

    reference = cr2_atoms(args.monomer)
    rotations = np.empty((n_states, 2, 3, 3), dtype=np.float64)
    translations = np.empty((n_states, 2, 3), dtype=np.float64)
    rmsd = np.empty((n_states, 2), dtype=np.float64)
    common_names = None

    for state in range(1, n_states + 1):
        model = cmd.get_model("traj and resn CR2", state=state)
        sites = {key: {} for key in residue_keys}
        for atom in model.atom:
            key = ((atom.chain or "").strip(), str(atom.resi).strip())
            if key in sites:
                sites[key][atom.name.strip()] = np.asarray(atom.coord, dtype=float)

        for site_index, key in enumerate(residue_keys):
            common = sorted(set(reference) & set(sites[key]))
            if common_names is None:
                common_names = common
            elif common != common_names:
                raise RuntimeError(f"CR2 atom-name set changed in state {state}, site {key}")
            if len(common) < 4:
                raise RuntimeError(f"Only {len(common)} shared CR2 atoms in state {state}, site {key}")
            source = np.array([reference[name] for name in common])
            target = np.array([sites[key][name] for name in common])
            rotations[state - 1, site_index], translations[state - 1, site_index], rmsd[
                state - 1, site_index
            ] = kabsch(source, target)

        if state == 1 or state % 100 == 0 or state == n_states:
            print(f"[transforms] {state}/{n_states}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        rotation=rotations,
        translation=translations,
        rmsd_A=rmsd,
        common_names=np.asarray(common_names),
        residue_keys=np.asarray(residue_keys),
    )
    print(
        f"[transforms] wrote {args.out}; {n_states} frames, {len(common_names)} shared atoms, "
        f"fit RMSD mean {rmsd.mean():.4f} A (max {rmsd.max():.4f} A)"
    )


# PyMOL executes scripts with a non-standard module name, so run unconditionally.
main()
try:
    from pymol import cmd as _cmd
    _cmd.quit(0)
except Exception:
    pass

#!/usr/bin/env python3
"""Build candidate alternative tandem registers from the 1MYW crystal lattice.

The overnight steered-MD run refuted linker tension as the cause of the angle
discrepancy, and the superradiance sign analysis then sharpened the target:
the emitting state is BRIGHTER than the monomer, which requires cos(alpha) < 0,
so the anisotropy's |cos alpha| = 0.660 is alpha = 131.3 deg (obtuse), not the
48.7 deg Nguyen et al. quote. Both experiments therefore agree on an obtuse
angle -- the same side as our computed 96-101 deg, but ~30 deg further over.

notes/lattice_dimer_scan.py finds that the P3(1)12 lattice offers exactly one
packing contact with an obtuse angle near that target AND a C-term-to-N-term
span the 33-residue linker can bridge comfortably:

    op 5, lattice (1, 0, -1): alpha = 125.2 deg, |cos| = 0.576,
    CR2 separation 35.9 A, linker span 26.5 A, 55 heavy-atom contacts

against the biological dimer's alpha = 106.6 deg, |cos| = 0.285, separation
25.4 A, span 54.2 A, 229 contacts.

This script materialises those registers as dimer PDBs so they can be scored
with the same rigid STEOM density placement used for the production ensemble,
and then solvated for MD if they survive.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

PDB = Path("/home/robson/PetaChem/1MYW.pdb")


def cell_matrix(a, b, c, al, be, ga):
    al, be, ga = np.radians([al, be, ga])
    v = np.sqrt(1 - np.cos(al) ** 2 - np.cos(be) ** 2 - np.cos(ga) ** 2
                + 2 * np.cos(al) * np.cos(be) * np.cos(ga))
    return np.array([
        [a, b * np.cos(ga), c * np.cos(be)],
        [0, b * np.sin(ga), c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)],
        [0, 0, c * v / np.sin(ga)],
    ])


def symmetry_ops(path):
    cur = {}
    for line in open(path):
        if line.startswith("REMARK 290   SMTRY"):
            k = int(line[18]) - 1
            idx = int(line[19:23])
            vals = [float(x) for x in line[23:].split()[:4]]
            cur.setdefault(idx, np.zeros((3, 4)))[k] = vals
    return [(cur[i][:, :3], cur[i][:, 3]) for i in sorted(cur)]


def read_atom_lines(path):
    keep = []
    for line in open(path):
        if line[:6] in ("ATOM  ", "HETATM") and line[17:20].strip() != "HOH":
            keep.append(line.rstrip("\n"))
    return keep


def transform_line(line, rot, trans, chain):
    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    new = rot @ xyz + trans
    return (line[:21] + chain + line[22:30]
            + f"{new[0]:8.3f}{new[1]:8.3f}{new[2]:8.3f}" + line[54:])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--op", type=int, required=True, help="1-based SMTRY operator index")
    ap.add_argument("--lattice", type=int, nargs=3, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    lines = read_atom_lines(PDB)
    for raw in open(PDB):
        if raw.startswith("CRYST1"):
            cryst1 = raw.rstrip("\n")
            cell = cell_matrix(*[float(raw[i:i + 9]) for i in (6, 15, 24)],
                               *[float(raw[i:i + 7]) for i in (33, 40, 47)])
            break

    rot, trans = symmetry_ops(PDB)[args.op - 1]
    offset = cell @ np.array(args.lattice, dtype=float)
    total_trans = trans + offset

    out_lines = [cryst1]
    out_lines += [line[:21] + "A" + line[22:] for line in lines]
    out_lines.append("TER")
    out_lines += [transform_line(line, rot, total_trans, "B") for line in lines]
    out_lines.append("TER")
    out_lines.append("END")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out_lines) + "\n")

    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    xyz_b = xyz @ rot.T + total_trans
    print(f"op {args.op}, lattice {tuple(args.lattice)} -> {args.out}")
    print(f"  chain A centroid {np.round(xyz.mean(0), 2)}")
    print(f"  chain B centroid {np.round(xyz_b.mean(0), 2)}")
    print(f"  centroid separation {np.linalg.norm(xyz.mean(0) - xyz_b.mean(0)):.2f} A")
    print(f"  {len(lines)} atoms per chain")


if __name__ == "__main__":
    main()

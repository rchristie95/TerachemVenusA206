#!/usr/bin/env python3
"""Enumerate every crystal-lattice dimer of 1MYW and score it.

The tandem construct was built in the register of the author-determined
biological dimer (REMARK 350 BIOMT2). That register puts the 33-residue linker
under tension and gives an inter-chromophore transition-dipole angle of ~107 deg
against ~131 deg required by the limiting anisotropy. This asks whether the
lattice offers any other packing contact that is a better candidate: a real
interface (buried contact area), a chromophore separation near the ~25 A that
both papers assume, an angle nearer the experimental requirement, and a
C-term-to-N-term span the linker can bridge without strain.

Applies all six P3(1)12 symmetry operators over a 3x3x3 block of lattice
translations, keeps images that actually touch chain A, and reports one row per
distinct interface.
"""

from __future__ import annotations

import itertools
import numpy as np
from scipy.spatial import cKDTree

PDB = "/home/robson/PetaChem/1MYW.pdb"
CONTACT_A = 4.5          # heavy-atom contact cutoff
MIN_CONTACTS = 40        # below this it is a lattice brush, not an interface
RING = ["CA2", "C2", "N2", "C1", "N3"]


def read_pdb(path):
    rows = []
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if line[17:20].strip() == "HOH":
            continue
        element = line[76:78].strip().upper()
        if element == "H":
            continue
        rows.append((int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                     float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return rows


def cell_matrix(a, b, c, al, be, ga):
    al, be, ga = np.radians([al, be, ga])
    v = np.sqrt(1 - np.cos(al)**2 - np.cos(be)**2 - np.cos(ga)**2
                + 2*np.cos(al)*np.cos(be)*np.cos(ga))
    return np.array([
        [a, b*np.cos(ga), c*np.cos(be)],
        [0, b*np.sin(ga), c*(np.cos(al) - np.cos(be)*np.cos(ga))/np.sin(ga)],
        [0, 0, c*v/np.sin(ga)],
    ])


def symmetry_ops(path):
    rot, trans, cur = [], [], {}
    for line in open(path):
        if line.startswith("REMARK 290   SMTRY"):
            k = int(line[18]) - 1
            idx = int(line[19:23])
            vals = [float(x) for x in line[23:].split()[:4]]
            cur.setdefault(idx, np.zeros((3, 4)))[k] = vals
    for idx in sorted(cur):
        rot.append(cur[idx][:, :3])
        trans.append(cur[idx][:, 3])
    return list(zip(rot, trans))


def chromophore(rows):
    d = {n: np.array([x, y, z]) for r, rn, n, x, y, z in rows if rn == "CR2"}
    ring = [d[r] for r in RING if r in d]
    axis = d["OH"] - np.mean(ring, axis=0)
    return np.mean(list(d.values()), axis=0), axis / np.linalg.norm(axis)


def main():
    rows = read_pdb(PDB)
    X = np.array([r[3:] for r in rows])
    for line in open(PDB):
        if line.startswith("CRYST1"):
            cell = cell_matrix(*[float(line[i:i+9]) for i in (6, 15, 24)],
                               *[float(line[i:i+7]) for i in (33, 40, 47)])
            break
    ops = symmetry_ops(PDB)
    cenA, axA = chromophore(rows)
    resid = np.array([r[0] for r in rows])
    nterm, cterm = resid.min(), resid.max()
    ca = {r[0]: np.array(r[3:]) for r in rows if r[2] == "CA"}
    treeA = cKDTree(X)

    print(f"{len(ops)} symmetry operators, {len(rows)} heavy atoms, "
          f"residues {nterm}-{cterm}")
    print(f"reference chromophore centroid {np.round(cenA,2)}\n")
    print(f"{'op':>3} {'lattice':>10} {'contacts':>9} {'CR2 sep':>8} "
          f"{'angle':>7} {'|cos a|':>8} {'Cterm-Nterm':>12}  note")

    seen = []
    inv_cell = np.linalg.inv(cell)
    for k, (R, t) in enumerate(ops, 1):
        # Centre the translation search on the shift that brings this operator's
        # image closest to chain A; the deposited molecule sits far from the
        # origin, so a fixed block around (0,0,0) misses real contacts entirely
        # -- the biological dimer itself needs (2,2,-1).
        base = np.round(inv_cell @ (X.mean(0) - (X @ R.T + t).mean(0))).astype(int)
        for delta in itertools.product((-1, 0, 1), repeat=3):
            shift = tuple(base + np.array(delta))
            off = cell @ np.array(shift, float)
            Y = X @ R.T + t + off
            if np.linalg.norm(Y.mean(0) - X.mean(0)) < 1e-6:
                continue
            pairs = treeA.query_ball_point(Y, CONTACT_A)
            n_contact = sum(len(p) for p in pairs)
            if n_contact < MIN_CONTACTS:
                continue
            cenB = cenA @ R.T + t + off
            axB = axA @ R.T
            cos = float(np.dot(axA, axB))
            sep = float(np.linalg.norm(cenA - cenB))
            span = float(np.linalg.norm(ca[cterm] - (ca[nterm] @ R.T + t + off)))
            key = (round(sep, 1), round(cos, 2), n_contact)
            if key in seen:
                continue
            seen.append(key)
            interface = {rows[i][0] for p in pairs for i in p}
            note = "A206 interface" if 206 in interface else ""
            print(f"{k:>3} {str(shift):>10} {n_contact:>9} {sep:>8.2f} "
                  f"{np.degrees(np.arccos(cos)):>7.2f} {abs(cos):>8.3f} "
                  f"{span:>12.1f}  {note}")


if __name__ == "__main__":
    main()

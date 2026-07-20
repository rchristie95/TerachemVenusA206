#!/usr/bin/env python3
r"""
tandem_to_coupling.py  --  Split a single-molecule tandem-dimer NVT trajectory
into the two-chain (barrel1 -> chain A, barrel2 -> chain B) PBC-whole convention
that coupling_ensemble.py expects.

The tandem is one covalent chain, so run_nvt.py's PDBFixer fragments it into
several chains at the CR2 chromophores, with one chain spanning BOTH barrels
(FP1 tail + linker + FP2 head). Residue/chain ranges are therefore unreliable.
Instead we split SPATIALLY: each protein atom is assigned to the nearer of the
two CR2 chromophore centroids (with periodic imaging), keeping only atoms within
a barrel radius so the flexible linker is dropped.

    python tandem_to_coupling.py --in tc_tandem_nvt/tandem_nvt.pdb \
        --out tandem_nvt_clean.pdb --rcut 27
"""
import argparse
from pathlib import Path

import numpy as np

WATER_IONS = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "CA", "ZN", "CL-", "NA+"}


def parse_box(lines):
    for l in lines:
        if l.startswith("CRYST1"):
            return np.array([float(l[6:15]), float(l[15:24]), float(l[24:33])])
    return None


def split_models(lines):
    if not any(l.startswith("MODEL") for l in lines):
        yield [l for l in lines if l[:6].strip() in ("ATOM", "HETATM")]
        return
    cur = None
    for l in lines:
        if l.startswith("MODEL"):
            cur = []
        elif l.startswith("ENDMDL"):
            if cur is not None:
                yield cur
            cur = None
        elif cur is not None and l[:6].strip() in ("ATOM", "HETATM"):
            cur.append(l)


def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def set_xyz(l, p, chain):
    return f"{l[:21]}{chain}{l[22:30]}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{l[54:]}"


def two_cr2_centroids(prot, box):
    """The two CR2 chromophore centroids (min-imaged relative to each other)."""
    cr2 = {}
    for l in prot:
        if l[17:20].strip() == "CR2":
            cr2.setdefault(l[21], []).append(xyz(l))
    keys = sorted(cr2, key=lambda k: -len(cr2[k]))[:2]
    c1 = np.mean(cr2[keys[0]], axis=0)
    c2 = np.mean(cr2[keys[1]], axis=0)
    if box is not None:
        c2 = c2 - box * np.round((c2 - c1) / box)
    return c1, c2


def process_frame(atom_lines, box, n_fp1, n_link):
    """Split by residue FILE ORDER (connectivity is preserved through PDBFixer):
    the first n_fp1 protein residues are barrel A, the next n_link are the linker
    (dropped), the rest are barrel B. Then make each barrel PBC-whole around its
    CR2 and image B next to A."""
    prot = [l for l in atom_lines if l[17:20].strip() not in WATER_IONS]
    # residues in order of first appearance
    order = []
    groups = {}
    for l in prot:
        key = (l[21], l[22:27])
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(l)

    segA = [groups[k] for k in order[:n_fp1]]
    segB = [groups[k] for k in order[n_fp1 + n_link:]]

    def assemble(seg, chain):
        atoms = [l for g in seg for l in g]
        cr2 = [xyz(l) for l in atoms if l[17:20].strip() == "CR2"]
        ref = np.mean(cr2, axis=0) if cr2 else np.mean([xyz(l) for l in atoms], axis=0)
        coords = np.array([xyz(l) for l in atoms])
        if box is not None:                       # make barrel whole around its CR2
            coords = coords - box * np.round((coords - ref) / box)
        return atoms, coords, ref

    aA, cA, refA = assemble(segA, "A")
    aB, cB, refB = assemble(segB, "B")
    if box is not None:                            # image barrel B next to barrel A
        cr2B = np.mean([cB[i] for i, l in enumerate(aB) if l[17:20].strip() == "CR2"], axis=0)
        cB = cB - box * np.round((cr2B - refA) / box)

    out = []
    for l, p in zip(aA, cA):
        out.append(set_xyz(l, p, "A"))
    for l, p in zip(aB, cB):
        out.append(set_xyz(l, p, "B"))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-fp1", type=int, default=229, help="residues in barrel 1 (chain A)")
    ap.add_argument("--n-link", type=int, default=33, help="linker residues (dropped)")
    args = ap.parse_args(argv)

    lines = Path(args.inp).read_text().splitlines(keepends=True)
    box = parse_box(lines)
    cryst = next((l for l in lines if l.startswith("CRYST1")), None)

    n = 0
    a0 = b0 = 0
    with open(args.out, "w") as f:
        if cryst:
            f.write(cryst)
        for atoms in split_models(lines):
            kept = process_frame(atoms, box, args.n_fp1, args.n_link)
            n += 1
            if n == 1:
                a0 = sum(1 for l in kept if l[21] == "A")
                b0 = sum(1 for l in kept if l[21] == "B")
            f.write(f"MODEL     {n:4d}\n")
            f.writelines(k if k.endswith("\n") else k + "\n" for k in kept)
            f.write("ENDMDL\n")
        f.write("END\n")
    print(f"[tandem-split] wrote {args.out}: {n} frames, chain A={a0} atoms, "
          f"chain B={b0} atoms (linker+water/ions dropped, barrels PBC-whole).")


if __name__ == "__main__":
    main()

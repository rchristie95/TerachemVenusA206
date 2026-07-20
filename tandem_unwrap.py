#!/usr/bin/env python3
r"""
tandem_unwrap.py  --  Make a single-molecule tandem NVT trajectory PBC-whole for
visualisation, keeping the whole construct (FP1 + linker + FP2).

The tandem is one covalent chain but FP2 swings out to ~70 A (> half the ~90 A
box), so simple imaging around one point wraps it. We instead unwrap SEQUENTIALLY
along the covalent chain (residue file order is the covalent order through
PDBFixer): each residue is imaged next to the previous one, so consecutive
peptide-bonded residues (~3.8 A apart) stay continuous and the fully-extended
linker unwraps correctly. Waters/ions are dropped; the three segments are
relabelled FP1 -> chain A, linker -> chain B, FP2 -> chain C for a clean
camera fit.

    python tandem_unwrap.py --in tc_tandem_nvt/tandem_nvt.pdb \
        --out tandem_nvt_whole.pdb --n-fp1 229 --n-link 33
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


def set_line(l, p, chain, resnum=None):
    resfield = f"{resnum:4d} " if resnum is not None else l[22:27]
    return f"{l[:21]}{chain}{resfield}{l[27:30]}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{l[54:]}"


def renumber_chains(lines):
    """Give each output chain consecutive residue numbers so PyMOL cartoon
    connects the backbone (the relabelled fragments otherwise reset per chain)."""
    out = []
    counters = {}
    last = {}
    for l in lines:
        ch = l[21]
        key = l[22:27]
        if ch not in counters:
            counters[ch] = 0
            last[ch] = None
        if last[ch] != key:
            counters[ch] += 1
            last[ch] = key
        out.append(f"{l[:22]}{counters[ch]:4d} {l[27:]}")
    return out


def cr2_ref(keys, groups):
    """First CR2 atom of a segment (a single, unambiguous reference point)."""
    for k in keys:
        for l in groups[k]:
            if l[17:20].strip() == "CR2":
                return xyz(l)
    return xyz(groups[keys[0]][0])


def process_frame(atom_lines, box, n_fp1, n_link):
    """Unwrap along the covalent chain (file order = covalent order through
    PDBFixer). For each residue: (1) make it internally whole around its own first
    atom (the CR2 chromophore straddles the box in most frames, so this is
    essential), then (2) image it next to the PREVIOUS residue. Following the
    FP1 -> linker -> FP2 path this way gives the true, linker-connected geometry
    without the min-image errors that corrupted the near-half-box frames.
    FP1 -> chain A, linker -> B, FP2 -> C."""
    prot = [l for l in atom_lines if l[17:20].strip() not in WATER_IONS]
    order, groups = [], {}
    for l in prot:
        key = (l[21], l[22:27])
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(l)

    out = []
    prev = None
    for i, key in enumerate(order):
        atoms = groups[key]
        coords = np.array([xyz(l) for l in atoms])
        if box is not None:
            coords = coords - box * np.round((coords - coords[0]) / box)   # (1) residue whole
            cen = coords.mean(0)
            if prev is not None:                                            # (2) next to previous
                shift = -box * np.round((cen - prev) / box)
                coords = coords + shift
                cen = cen + shift
            prev = cen
        chain = "A" if i < n_fp1 else ("B" if i < n_fp1 + n_link else "C")
        for l, p in zip(atoms, coords):
            out.append(set_line(l, p, chain))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-fp1", type=int, default=229)
    ap.add_argument("--n-link", type=int, default=33)
    args = ap.parse_args(argv)

    lines = Path(args.inp).read_text().splitlines(keepends=True)
    box = parse_box(lines)
    n = 0
    with open(args.out, "w") as f:
        for atoms in split_models(lines):
            kept = renumber_chains(process_frame(atoms, box, args.n_fp1, args.n_link))
            n += 1
            f.write(f"MODEL     {n:4d}\n")
            f.writelines(k if k.endswith("\n") else k + "\n" for k in kept)
            f.write("ENDMDL\n")
        f.write("END\n")
    print(f"[tandem-unwrap] wrote {args.out}: {n} frames, protein made PBC-whole "
          f"(FP1->A, linker->B, FP2->C).")


if __name__ == "__main__":
    main()

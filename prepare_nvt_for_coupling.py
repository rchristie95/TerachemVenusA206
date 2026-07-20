#!/usr/bin/env python3
"""
prepare_nvt_for_coupling.py

Convert a solvated OpenMM NVT trajectory (from run_dimer_nvt.py) into the clean,
PBC-whole, two-chain protein-only convention that sample_coupling_md.py / the
Stage-3 TDC alignment expect (matching venus_dimer.pdb: monomer 1 -> chain A,
monomer 2 -> chain B).

Two corrections are applied per frame:

  1. Periodic imaging. OpenMM's PDBReporter wraps molecules into the primary
     box, which splits each ~42 A barrel across the periodic boundary. A
     wrapped monomer superposes with good global RMSD yet maps the (central)
     chromophore transition density tens of Angstrom away. We therefore make
     each monomer whole by minimum-imaging every atom relative to that
     monomer's CR2 chromophore centroid (the barrel radius < L/2, so this is
     safe), then image monomer 2 next to monomer 1.

  2. Chain/solvent cleanup. PDBFixer splits each monomer across several chain
     IDs at the non-standard CR2 residue; water/ions get their own chains. We
     drop water/ions and collapse the protein chains of each monomer into a
     single chain (A, B). (PyMOL `super` is sequence-independent, so residue
     numbering need not be unique.)

Usage:
    python prepare_nvt_for_coupling.py --in tc_dimer_nvt_restrained/dimer_nvt_restrained.pdb \
        --out dimer_nvt_clean.pdb
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
    have = any(l.startswith("MODEL") for l in lines)
    if not have:
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


def protein_chains_in_order(atom_lines):
    order = []
    for l in atom_lines:
        if l[17:20].strip() in WATER_IONS:
            continue
        ch = l[21]
        if ch not in order:
            order.append(ch)
    return order


def build_chain_map(atom_lines, chains_a, chains_b):
    if chains_a and chains_b:
        amap = {c: "A" for c in chains_a}
        amap.update({c: "B" for c in chains_b})
        return amap
    chains = protein_chains_in_order(atom_lines)
    if len(chains) % 2 != 0:
        raise SystemExit(f"Cannot auto-split protein chains {chains}; pass --chains-a/--chains-b.")
    half = len(chains) // 2
    return {**{c: "A" for c in chains[:half]}, **{c: "B" for c in chains[half:]}}


def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def set_xyz(l, p):
    return f"{l[:30]}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{l[54:]}"


def process_frame(atom_lines, chain_map, box):
    """Return cleaned, PBC-whole protein atom lines relabelled to chains A/B."""
    # Partition protein atoms by target monomer.
    mono = {"A": [], "B": []}
    for l in atom_lines:
        if l[17:20].strip() in WATER_IONS:
            continue
        tgt = chain_map.get(l[21])
        if tgt:
            mono[tgt].append(l)

    out = []
    centroids = {}
    for tgt in ("A", "B"):
        lines = mono[tgt]
        coords = np.array([xyz(l) for l in lines])
        cr2 = np.array([xyz(l) for l in lines if l[17:20].strip() == "CR2"])
        ref = cr2.mean(axis=0) if len(cr2) else coords.mean(axis=0)
        if box is not None:
            # Make this monomer whole: image every atom near its chromophore.
            coords = coords - box * np.round((coords - ref) / box)
        centroids[tgt] = (np.array([xyz(l) for l in lines if l[17:20].strip() == "CR2"]) if len(cr2) else coords)
        mono[tgt] = (lines, coords)

    # Image monomer B next to monomer A using CR2 centroids (recomputed post-unwrap).
    def cr2_centroid(tgt):
        lines, coords = mono[tgt]
        idx = [i for i, l in enumerate(lines) if l[17:20].strip() == "CR2"]
        return coords[idx].mean(axis=0) if idx else coords.mean(axis=0)

    if box is not None:
        cA = cr2_centroid("A")
        linesB, coordsB = mono["B"]
        cB = coordsB[[i for i, l in enumerate(linesB) if l[17:20].strip() == "CR2"]].mean(axis=0)
        shift = box * np.round((cB - cA) / box)
        mono["B"] = (linesB, coordsB - shift)

    for tgt in ("A", "B"):
        lines, coords = mono[tgt]
        for l, p in zip(lines, coords):
            out.append(set_xyz(l[:21] + tgt + l[22:], p))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--chains-a", default="")
    ap.add_argument("--chains-b", default="")
    args = ap.parse_args()

    lines = Path(args.inp).read_text().splitlines(keepends=True)
    box = parse_box(lines)
    chains_a = [c.strip() for c in args.chains_a.split(",") if c.strip()]
    chains_b = [c.strip() for c in args.chains_b.split(",") if c.strip()]

    models = list(split_models(lines))
    if not models:
        raise SystemExit("No models/atoms found.")
    chain_map = build_chain_map(models[0], chains_a, chains_b)
    print(f"[*] {len(models)} model(s); box={None if box is None else box.round(1)}; "
          f"chain map {chain_map}")

    cryst = next((l for l in lines if l.startswith("CRYST1")), None)
    with open(args.out, "w") as f:
        if cryst:
            f.write(cryst)
        for i, atoms in enumerate(models, start=1):
            kept = process_frame(atoms, chain_map, box)
            f.write(f"MODEL     {i:4d}\n")
            f.writelines(kept if kept and kept[0].endswith("\n") else [k + "\n" for k in kept])
            f.write("ENDMDL\n")
        f.write("END\n")

    k0 = process_frame(models[0], chain_map, box)
    a = sum(1 for l in k0 if l[21] == "A")
    b = sum(1 for l in k0 if l[21] == "B")
    print(f"[*] Wrote {args.out}: per frame chain A={a} atoms, chain B={b} atoms "
          f"(water/ions dropped, monomers made PBC-whole).")


if __name__ == "__main__":
    main()

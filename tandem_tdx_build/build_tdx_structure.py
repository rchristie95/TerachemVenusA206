#!/usr/bin/env python3
"""Build the dVenus-TDX starting structure from the dVenus-TD v4 structure.

TDX is TD with one barrel's chromophore extinguished (Amber, Venus-local Y67C).
The two structures differ in one place only: the dark barrel's mature CR2
becomes the unreacted Gly-Cys-Gly triad it would have been had maturation never
happened, which is one residue in the sequence and three in the PDB.  Everything else -- both barrel folds,
the linker, and above all the A206 inter-barrel register the coupling depends
on -- is carried over atom-for-atom from

    tandem_v4_correct_seq/tandem_dimer_v4_af3tails.pdb

That register is deliberately NOT taken from AlphaFold.  alphafold_tandem/
ALPHAFOLD_REGISTER.md establishes that AF3 folds each barrel to 0.70 A and then
reports 18.7-21.8 A PAE between them; its register sits 11.7 A from crystal.
Only LOCAL conformation is ever imported from AF, which is what --from-af3 does
here and what tandem_v4_correct_seq/graft_af3_tails.py did for the C-terminal
tails.

Two sources for the dead triad's geometry:

  --from-af3   Take the triad from an AlphaFold3 model of the TDX sequence
               (alphafold_tandem/af3_tdx_job_request_dark2.json), superposed on
               the 15 flanking barrel residues either side, backbone atoms only.
               AF folds the UNREACTED triad natively -- open ring, sp3 CA, real
               phi/psi -- so nothing has to be invented.  This is the path.

  default      Un-mature the existing CR2 in place, by atom correspondence:

                   Gly  N  CA  C        <- N1  CA1  C1
                   Cys  N  CA  CB  C O  <- N2  CA2  CB2 C2 O2
                   Gly  N  CA  C  O     <- N3  CA3  C3  O3

               with the phenol dropped, SG placed on the tyrosine's CB->CG
               vector so chi1 is kept, and the first glycine's carbonyl O -- the
               one lost as water when the ring closed -- rebuilt on C1's free
               sp2 lobe.  It does not work, and the script refuses to write the
               result unless you pass --force-strained.  The reason is worth
               stating because it is not obvious: in the mature chromophore C1
               already HAS three substituents, CA1, N2 and the ring N3, so the
               lobe the carbonyl O has to occupy is the one N3 is sitting in.
               The rebuilt O lands 0.21 A from N3.  Opening the ring is a
               conformational change of the triad inside a fixed barrel, not an
               atom substitution, and no amount of minimisation should be asked
               to start from a 0.21 A overlap.  Hence AlphaFold.

Usage:
    python3 tandem_tdx_build/build_tdx_structure.py --from-af3 fold_dvenus_tdx_model_0.cif
    python3 tandem_tdx_build/build_tdx_structure.py --from-af3 <cif> --dead-barrel 1
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SRC = REPO / "tandem_v4_correct_seq" / "tandem_dimer_v4_af3tails.pdb"

# Mature -> immature heavy-atom correspondence, per replacement residue.
TRIAD = (
    ("GLY", (("N", "N1"), ("CA", "CA1"), ("C", "C1"))),
    ("CYS", (("N", "N2"), ("CA", "CA2"), ("CB", "CB2"), ("C", "C2"), ("O", "O2"))),
    ("GLY", (("N", "N3"), ("CA", "CA3"), ("C", "C3"), ("O", "O3"))),
)
DROPPED = ("CD1", "CD2", "CE1", "CE2", "CG2", "CZ", "OH")  # the phenol, sorted

C_O_BOND = 1.231      # A, peptide carbonyl
CA_CB_BOND = 1.526    # A, sp3 (the mature CA2=CB2 is a short 1.35 A double bond)
CB_SG_BOND = 1.808    # A, cysteine
BACKBONE = ("N", "CA", "C", "O")


def unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def read_pdb(path: Path) -> list[dict]:
    """Group ATOM/HETATM lines into residues, in file order."""
    residues: list[dict] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        number, name = int(line[22:26]), line[17:20].strip()
        if not residues or residues[-1]["number"] != number or residues[-1]["name"] != name:
            residues.append({"number": number, "name": name, "chain": line[21], "atoms": []})
        residues[-1]["atoms"].append({
            "name": line[12:16].strip(),
            "element": line[76:78].strip() or line[12:16].strip()[0],
            "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        })
    return residues


def atom_map(residue: dict) -> dict[str, np.ndarray]:
    return {atom["name"]: atom["xyz"] for atom in residue["atoms"]}


def kabsch(mobile: np.ndarray, target: np.ndarray):
    mc, tc = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rotation, tc - rotation @ mc


def unmature_in_place(cr2: dict) -> list[tuple[str, list[tuple[str, np.ndarray]]]]:
    """Split one mature CR2 into the Gly-Cys-Gly triad it formed from."""
    xyz = atom_map(cr2)
    missing = {source for _, pairs in TRIAD for _, source in pairs} - set(xyz)
    if missing:
        raise SystemExit(f"CR2 {cr2['number']} is missing atoms {sorted(missing)}")

    built = []
    for name, pairs in TRIAD:
        built.append((name, [(target, xyz[source].copy()) for target, source in pairs]))

    # First glycine's carbonyl oxygen: sp2, in the CA-C-N plane, away from both.
    ca1, c1, n2 = xyz["CA1"], xyz["C1"], xyz["N2"]
    built[0][1].append(("O", c1 - C_O_BOND * unit(unit(ca1 - c1) + unit(n2 - c1))))

    # Cysteine: relax the double bond to sp3 length and put SG on the old
    # CB->CG vector, which keeps the tyrosine's chi1.
    ca2, cb2, cg2 = xyz["CA2"], xyz["CB2"], xyz["CG2"]
    cb = ca2 + CA_CB_BOND * unit(cb2 - ca2)
    sg = cb + CB_SG_BOND * unit(cg2 - cb2)
    cys = dict(built[1][1])
    cys["CB"] = cb
    built[1] = ("CYS", [("N", cys["N"]), ("CA", cys["CA"]), ("CB", cb),
                        ("SG", sg), ("C", cys["C"]), ("O", cys["O"])])
    return built


def graft_from_af3(cif: Path, residues: list[dict], dead_index: int,
                   anchor: int) -> list[tuple[str, list[tuple[str, np.ndarray]]]]:
    """Take the triad from an AF3 TDX model, superposed on flanking barrel residues."""
    af: dict[int, dict[str, np.ndarray]] = {}
    names: dict[int, str] = {}
    inloop = False
    for line in cif.read_text().splitlines():
        if line.startswith("_atom_site."):
            inloop = True
            continue
        if not inloop:
            continue
        parts = line.split()
        if len(parts) < 18:
            if af:
                break
            continue
        if parts[0] != "ATOM":
            continue
        number = int(parts[15])
        names[number] = parts[5]
        af.setdefault(number, {})[parts[3]] = np.array(
            [float(parts[10]), float(parts[11]), float(parts[12])])
    if not af:
        raise SystemExit(f"no ATOM records parsed from {cif}")

    # v4 numbering -> 511-aa numbering: each CR2 stands for three residues.
    def to_af(number: int) -> int | None:
        if number <= 65:
            return number
        if number == 66:
            return None
        if number <= 335:
            return number + 2
        if number == 336:
            return None
        return number + 4

    dead_number = residues[dead_index]["number"]
    triad_start = to_af(dead_number - 1) + 1
    triad = [triad_start, triad_start + 1, triad_start + 2]
    observed = [names.get(n) for n in triad]
    if observed != ["GLY", "CYS", "GLY"]:
        raise SystemExit(f"AF3 residues {triad} are {observed}, expected GLY CYS GLY -- "
                         f"is {cif.name} a model of the TDX sequence?")

    mobile, target = [], []
    window = [residues[i] for i in range(max(0, dead_index - anchor),
                                        min(len(residues), dead_index + anchor + 1))
              if residues[i]["name"] != "CR2"]
    for residue in window:
        number = to_af(residue["number"])
        if number is None or number not in af:
            continue
        local = atom_map(residue)
        for name in BACKBONE:
            if name in local and name in af[number]:
                mobile.append(af[number][name])
                target.append(local[name])
    mobile, target = np.array(mobile), np.array(target)
    rotation, shift = kabsch(mobile, target)
    rmsd = float(np.sqrt((((mobile @ rotation.T) + shift - target) ** 2).sum(1).mean()))
    print(f"  AF3 graft: anchored on {len(mobile)} backbone atoms of {len(window)} residues, "
          f"superposition RMSD {rmsd:.2f} A")
    if rmsd > 1.5:
        raise SystemExit(f"anchor RMSD {rmsd:.2f} A is too large to trust the graft")

    built = []
    for name, number in zip(("GLY", "CYS", "GLY"), triad):
        atoms = [(atom, rotation @ position + shift)
                 for atom, position in af[number].items()
                 if not atom.startswith("H")]
        built.append((name, atoms))
    return built


def triad_contacts(residues: list[dict], first: int, cutoff: float):
    """Non-bonded heavy-atom contacts made by the three new residues.

    Pairs inside one residue, and pairs between neighbouring residues, are
    bonds and their 1-3 partners; everything else within `cutoff` is a clash.
    """
    triad = range(first, first + 3)
    positions = [(index, residue["name"], atom["name"], atom["xyz"])
                 for index, residue in enumerate(residues) for atom in residue["atoms"]]
    coordinates = np.array([entry[3] for entry in positions])
    found = []
    for index in triad:
        for _, name, atom, xyz in [e for e in positions if e[0] == index]:
            distances = np.linalg.norm(coordinates - xyz, axis=1)
            for other in np.where(distances < cutoff)[0]:
                partner = positions[other]
                if abs(partner[0] - index) <= 1:
                    continue
                if partner[0] in triad and partner[0] < index:
                    continue        # already reported from the other side
                found.append((float(distances[other]),
                              f"{name}{residues[index]['number']} {atom}",
                              f"{partner[1]}{residues[partner[0]]['number']} {partner[2]}"))
    return sorted(found)


def write_pdb(path: Path, residues: list[dict]) -> None:
    serial = 0
    with path.open("w") as handle:
        for index, residue in enumerate(residues, start=1):
            record = "HETATM" if residue["name"] == "CR2" else "ATOM  "
            for atom in residue["atoms"]:
                serial += 1
                x, y, z = atom["xyz"]
                # Columns 13-16 hold the atom name with the element right-justified
                # in 13-14, so a one-letter element takes a leading space.
                name = atom["name"]
                field = name if len(name) == 4 else f" {name:<3s}"
                handle.write(
                    f"{record}{serial:5d} {field} {residue['name']:>3s} "
                    f"{residue['chain']}{index:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                    f"{atom['element']:>2s}\n")
        handle.write("END\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--dead-barrel", type=int, choices=(1, 2), default=2,
                    help="Barrel whose chromophore is extinguished (default 2)")
    ap.add_argument("--from-af3", type=Path,
                    help="AlphaFold3 .cif of the TDX sequence; the triad is taken from it")
    ap.add_argument("--anchor", type=int, default=15,
                    help="Residues either side used to superpose an AF3 graft")
    ap.add_argument("--output", type=Path)
    ap.add_argument("--clash-a", type=float, default=2.2,
                    help="Report non-bonded contacts to the new triad below this (A)")
    ap.add_argument("--force-strained", action="store_true",
                    help="Write the structure even if the new triad clashes")
    args = ap.parse_args()

    residues = read_pdb(args.source)
    cr2 = [i for i, residue in enumerate(residues) if residue["name"] == "CR2"]
    if len(cr2) != 2:
        raise SystemExit(f"{args.source} has {len(cr2)} CR2 residues, expected 2")
    dead_index = cr2[args.dead_barrel - 1]
    live_index = cr2[1 - (args.dead_barrel - 1)]
    dead = residues[dead_index]
    print(f"source {args.source}: {len(residues)} residues, CR2 at "
          f"{[residues[i]['number'] for i in cr2]}")
    print(f"extinguishing barrel {args.dead_barrel}: CR2 {dead['number']} -> GLY-CYS-GLY")

    if args.from_af3:
        built = graft_from_af3(args.from_af3, residues, dead_index, args.anchor)
        source_label = args.from_af3.name
    else:
        built = unmature_in_place(dead)
        source_label = "in-place un-maturation"

    live_before = atom_map(residues[live_index])
    dead_before = atom_map(dead)
    replacement = [{"number": dead["number"] + offset, "name": name, "chain": dead["chain"],
                    "atoms": [{"name": atom, "element": atom[0], "xyz": xyz}
                              for atom, xyz in atoms]}
                   for offset, (name, atoms) in enumerate(built)]
    residues[dead_index:dead_index + 1] = replacement

    output = args.output or (Path(__file__).resolve().parent /
                             f"tandem_dimer_tdx_dark{args.dead_barrel}.pdb")

    # --- checks -------------------------------------------------------------
    live_after = [r for r in residues if r["name"] == "CR2"]
    cys = next(r for r in replacement if r["name"] == "CYS")
    cys_xyz = atom_map(cys)
    mapped = {source for _, pairs in TRIAD for _, source in pairs}
    dropped = sorted(set(dead_before) - mapped)
    if args.from_af3 is None and tuple(dropped) != DROPPED:
        raise SystemExit(f"expected to drop the phenol {DROPPED}, dropped {tuple(dropped)}")
    kept = np.array([live_before[n] for n in sorted(live_before)])
    after = np.array([atom_map(live_after[0])[n] for n in sorted(live_before)])
    live_shift = float(np.abs(kept - after).max())

    print(f"  triad geometry from: {source_label}")
    print(f"  atoms dropped from the mature chromophore: {dropped}")
    print(f"  surviving CR2: {len(live_after)} residue, max coordinate change {live_shift:.3f} A")
    print(f"  chromophore separation CR2 CA2 <-> dead-triad Cys CA: "
          f"{np.linalg.norm(live_before['CA2'] - cys_xyz['CA']):.3f} A "
          f"(TD CR2-CR2 CA2-CA2 was "
          f"{np.linalg.norm(live_before['CA2'] - dead_before['CA2']):.3f} A)")
    gly_c = atom_map(replacement[0])["C"]
    gly_n = atom_map(replacement[2])["N"]
    opening = float(np.linalg.norm(gly_c - gly_n))
    print(f"  former ring bond Gly(1) C ... Gly(3) N: {opening:.3f} A "
          f"({'still closed' if opening < 2.5 else 'open'})")

    contacts = triad_contacts(residues, dead_index, args.clash_a)
    for distance, left, right in contacts[:5]:
        print(f"  CLASH {distance:.2f} A: {left} ... {right}")
    if contacts and not args.force_strained:
        raise SystemExit(
            f"refusing to write: the new triad makes {len(contacts)} contact(s) below "
            f"{args.clash_a} A, the worst at {contacts[0][0]:.2f} A. Build the triad "
            f"from an AlphaFold model of the TDX sequence (--from-af3) instead of "
            f"un-maturing the chromophore in place; see this script's docstring for "
            f"why the in-place route cannot place the carbonyl oxygen. "
            f"--force-strained overrides.")

    write_pdb(output, residues)
    print(f"wrote {output}: {len(residues)} residues, "
          f"{sum(len(r['atoms']) for r in residues)} atoms")


if __name__ == "__main__":
    main()

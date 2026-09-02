#!/usr/bin/env python3
"""Build the dVenus-TDX starting structure from the dVenus-TD v4 structure.

TDX is TD with one barrel's chromophore extinguished.  Addgene #250939 does
that by DELETING the chromophore tyrosine (Y339, barrel-local 67), not by
mutating it -- checked against the deposit on 2026-09-02, see
alphafold_tandem/TDX_CONSTRUCT.md.  So the two structures differ in one place
only: the dark barrel's mature CR2 becomes the unreacted Gly-Gly it would have
been had the tyrosine ever been there to cyclise, which is one residue in the
PDB becoming two and the chain becoming one residue shorter than TD.
Everything else -- both barrel folds, the linker, and above all the A206
inter-barrel register the coupling depends on -- is carried over atom-for-atom
from

    tandem_v4_correct_seq/tandem_dimer_v4_af3tails.pdb

That register is deliberately NOT taken from AlphaFold.  alphafold_tandem/
ALPHAFOLD_REGISTER.md establishes that AF3 folds each barrel to 0.70 A and then
reports 18.7-21.8 A PAE between them; its register sits 11.7 A from crystal.
Only LOCAL conformation is ever imported from AF, which is what --from-af3 does
here and what tandem_v4_correct_seq/graft_af3_tails.py did for the C-terminal
tails.

Two sources for the dead pair's geometry:

  --from-af3   Take Gly-Gly from an AlphaFold3 model of the 510-aa TDX sequence
               (alphafold_tandem/af3_tdx_job_request_dark2.json), superposed on
               the 15 flanking barrel residues either side, backbone atoms only.
               AF folds the shortened strand natively, which matters more here
               than it did for the substitution reading of the construct:
               deleting a residue from a barrel strand is a register change
               inside the dark barrel, not a side-chain swap, and nothing local
               to the crystal structure knows where the shortened strand sits.
               This is the path.

  default      Collapse the existing CR2 in place, by atom correspondence:

                   Gly  N  CA  C  O     <- N1  CA1  C1  (O rebuilt)
                   Gly  N  CA  C  O     <- N3  CA3  C3  O3

               dropping the whole tyrosine-derived middle of the chromophore.
               Unlike the Y->C reading of the construct this is geometrically
               well posed, and it is worth knowing why.  In the mature
               chromophore the imidazolinone ring is C1-N2-CA2-C2-N3, so C1 is
               already bonded to N3: deleting the tyrosine turns that 1.355 A
               ring bond straight into the Gly-Gly peptide bond, and C1, having
               lost N2, has a free sp2 lobe for the carbonyl oxygen that the
               Y->C route had nowhere to put.  What it cannot do is relax the
               strand: the two glycines stay exactly where the ring held them,
               with ring-sized bond angles.  Useful as a cross-check on the
               graft and as a fallback, not as the production structure.

The script refuses to write a structure whose new residues clash; that refusal
is the point, and --force-strained is for inspecting a refused build, not for
simulating one.

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

# Mature CR2 -> the Gly-Gly the deleted-tyrosine triad leaves behind.
# The tyrosine-derived atoms (N2 CA2 CB2 CG2 CD1 CD2 CE1 CE2 CZ OH C2 O2) go.
DEAD_PAIR = (
    ("GLY", (("N", "N1"), ("CA", "CA1"), ("C", "C1"))),
    ("GLY", (("N", "N3"), ("CA", "CA3"), ("C", "C3"), ("O", "O3"))),
)
DROPPED = ("C2", "CA2", "CB2", "CD1", "CD2", "CE1", "CE2", "CG2", "CZ",
           "N2", "O2", "OH")            # the whole tyrosine, sorted
LIVE_TRIAD = 3                          # residues a surviving CR2 stands for
DEAD_PAIR_SIZE = len(DEAD_PAIR)         # residues the extinguished one stands for

C_O_BOND = 1.231      # A, peptide carbonyl
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


def sequence_numbering(residues: list[dict], dead_index: int) -> dict[int, int | None]:
    """Map each v4 residue's list index to its first residue number in the
    510-aa TDX chain.  A surviving CR2 stands for three sequence residues and
    the extinguished one for two, so the map is built by walking the chain
    rather than by hard-coded offsets -- which is what keeps it correct for
    either choice of dead barrel.  CR2 entries map to the number of the FIRST
    residue they stand for.
    """
    mapping: dict[int, int] = {}
    number = 1
    for index, residue in enumerate(residues):
        mapping[index] = number
        if residue["name"] == "CR2":
            number += DEAD_PAIR_SIZE if index == dead_index else LIVE_TRIAD
        else:
            number += 1
    return mapping


def collapse_in_place(cr2: dict) -> list[tuple[str, list[tuple[str, np.ndarray]]]]:
    """Split one mature CR2 into the Gly-Gly a deleted chromophore tyrosine leaves."""
    xyz = atom_map(cr2)
    missing = {source for _, pairs in DEAD_PAIR for _, source in pairs} - set(xyz)
    if missing:
        raise SystemExit(f"CR2 {cr2['number']} is missing atoms {sorted(missing)}")

    built = [(name, [(target, xyz[source].copy()) for target, source in pairs])
             for name, pairs in DEAD_PAIR]

    # The first glycine's carbonyl oxygen, lost as water when the ring closed.
    # C1 keeps CA1 and N3 once the tyrosine goes, so the third sp2 lobe is free:
    # this is exactly the placement the Y->C reading of the construct could not
    # make, because there N2 was still occupying it.
    ca1, c1, n3 = xyz["CA1"], xyz["C1"], xyz["N3"]
    built[0][1].append(("O", c1 - C_O_BOND * unit(unit(ca1 - c1) + unit(n3 - c1))))
    return built


def graft_from_af3(cif: Path, residues: list[dict], dead_index: int,
                   anchor: int) -> list[tuple[str, list[tuple[str, np.ndarray]]]]:
    """Take the dead pair from an AF3 TDX model, superposed on flanking barrel residues."""
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

    to_af = sequence_numbering(residues, dead_index)
    pair = [to_af[dead_index] + offset for offset in range(DEAD_PAIR_SIZE)]
    observed = [names.get(n) for n in pair]
    if observed != ["GLY"] * DEAD_PAIR_SIZE:
        raise SystemExit(
            f"AF3 residues {pair} are {observed}, expected "
            f"{['GLY'] * DEAD_PAIR_SIZE} -- is {cif.name} a model of the 510-aa TDX "
            f"sequence? A model of TD, or of the 511-aa Y->C reading of the "
            f"construct, has a residue between the two glycines and will land here.")
    if len(af) != len(to_af) - 2 + LIVE_TRIAD + DEAD_PAIR_SIZE:
        expected = len(to_af) - 2 + LIVE_TRIAD + DEAD_PAIR_SIZE
        raise SystemExit(f"{cif.name} has {len(af)} residues, expected {expected} "
                         f"for the TDX chain")

    mobile, target = [], []
    window = [residues[i] for i in range(max(0, dead_index - anchor),
                                        min(len(residues), dead_index + anchor + 1))
              if residues[i]["name"] != "CR2"]
    for residue in window:
        index = residues.index(residue)
        number = to_af[index]
        if number not in af:
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

    return [(names[number],
             [(atom, rotation @ position + shift)
              for atom, position in af[number].items() if not atom.startswith("H")])
            for number in pair]


def new_residue_contacts(residues: list[dict], first: int, count: int, cutoff: float):
    """Non-bonded heavy-atom contacts made by the newly built residues.

    Pairs inside one residue, and pairs between neighbouring residues, are
    bonds and their 1-3 partners; everything else within `cutoff` is a clash.
    """
    new = range(first, first + count)
    positions = [(index, residue["name"], atom["name"], atom["xyz"])
                 for index, residue in enumerate(residues) for atom in residue["atoms"]]
    coordinates = np.array([entry[3] for entry in positions])
    found = []
    for index in new:
        for _, name, atom, xyz in [e for e in positions if e[0] == index]:
            distances = np.linalg.norm(coordinates - xyz, axis=1)
            for other in np.where(distances < cutoff)[0]:
                partner = positions[other]
                if abs(partner[0] - index) <= 1:
                    continue
                if partner[0] in new and partner[0] < index:
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
                    help="Barrel whose chromophore is extinguished (default 2, "
                         "which is the barrel Addgene #250939 darkens)")
    ap.add_argument("--from-af3", type=Path,
                    help="AlphaFold3 .cif of the 510-aa TDX sequence; the dead "
                         "pair is taken from it")
    ap.add_argument("--anchor", type=int, default=15,
                    help="Residues either side used to superpose an AF3 graft")
    ap.add_argument("--output", type=Path)
    ap.add_argument("--clash-a", type=float, default=2.2,
                    help="Report non-bonded contacts to the new residues below this (A)")
    ap.add_argument("--force-strained", action="store_true",
                    help="Write the structure even if the new residues clash")
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
    print(f"extinguishing barrel {args.dead_barrel}: CR2 {dead['number']} -> GLY-GLY "
          f"(chromophore tyrosine deleted, Addgene #250939)")

    if args.from_af3:
        built = graft_from_af3(args.from_af3, residues, dead_index, args.anchor)
        source_label = args.from_af3.name
    else:
        built = collapse_in_place(dead)
        source_label = "in-place collapse of the mature chromophore"

    before = {residue["number"]: atom_map(residue) for residue in residues}
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
    mapped = {source for _, pairs in DEAD_PAIR for _, source in pairs}
    dropped = tuple(sorted(set(dead_before) - mapped))
    if dropped != DROPPED:
        raise SystemExit(f"expected to drop the tyrosine {DROPPED}, dropped {dropped}")
    kept = np.array([live_before[n] for n in sorted(live_before)])
    after = np.array([atom_map(live_after[0])[n] for n in sorted(live_before)])
    live_shift = float(np.abs(kept - after).max())

    # Everything that is not one of the two new residues must be untouched;
    # this, not any single distance, is what "the crystal register is intact"
    # means, and it is checked over every atom rather than sampled.
    moved = 0.0
    for index, residue in enumerate(residues):
        if dead_index <= index < dead_index + DEAD_PAIR_SIZE:
            continue
        reference = before[residue["number"]]
        for atom in residue["atoms"]:
            moved = max(moved, float(np.abs(atom["xyz"] - reference[atom["name"]]).max()))

    print(f"  dead pair geometry from: {source_label}")
    print(f"  atoms dropped from the mature chromophore: {list(dropped)}")
    print(f"  surviving CR2: {len(live_after)} residue, max coordinate change {live_shift:.3f} A")
    print(f"  all other residues: max coordinate change {moved:.3f} A "
          f"(the inter-barrel register is carried over, not rebuilt)")
    gly1, gly2 = (atom_map(r) for r in replacement)
    print(f"  chromophore separation CR2 CA2 <-> dead Gly(1) CA: "
          f"{np.linalg.norm(live_before['CA2'] - gly1['CA']):.3f} A "
          f"(TD CR2-CR2 CA2-CA2 reference "
          f"{np.linalg.norm(live_before['CA2'] - dead_before['CA2']):.3f} A)")
    bond = float(np.linalg.norm(gly1["C"] - gly2["N"]))
    print(f"  Gly(1) C - Gly(2) N: {bond:.3f} A "
          f"({'peptide bond' if 1.20 <= bond <= 1.50 else 'NOT a peptide bond'})")
    if not 1.20 <= bond <= 1.50:
        raise SystemExit(
            f"the two glycines are {bond:.3f} A apart, which is not a peptide bond. "
            f"With the tyrosine deleted they are consecutive residues; a gap here "
            f"means the model is of a different construct or the graft is misaligned.")

    # Where the graft meets the crystal.  The dead pair spans a slot the crystal
    # cut for a three-residue triad, so a shortened strand does not fit it
    # exactly and the mismatch lands entirely in these two bonds.  It is small --
    # the mature triad is ring-compressed, N1..C3 only 4.3 A, so Gly-Gly is
    # nearly the same length -- but it is real, and it is strain the
    # minimisation before MD has to take out, not something to leave unsaid.
    flanks = []
    for offset, (left, right, label) in enumerate((
            (dead_index - 1, dead_index, "N-terminal"),
            (dead_index + DEAD_PAIR_SIZE - 1, dead_index + DEAD_PAIR_SIZE, "C-terminal"))):
        if not 0 <= left or right >= len(residues):
            continue
        before_atoms, after_atoms = atom_map(residues[left]), atom_map(residues[right])
        if "C" not in before_atoms or "N" not in after_atoms:
            continue
        flanks.append((label, float(np.linalg.norm(before_atoms["C"] - after_atoms["N"]))))
    for label, distance in flanks:
        note = "relaxed" if distance <= 1.45 else "stretched, minimisation will take it out"
        print(f"  {label} junction C-N: {distance:.3f} A ({note})")
    worst = max((d for _, d in flanks), default=0.0)
    if worst > 1.8:
        raise SystemExit(
            f"a junction C-N bond is {worst:.3f} A, which is a break rather than a "
            f"stretched bond. The dead pair does not fit the slot the crystal "
            f"structure leaves; check that the model is of the 510-aa TDX sequence "
            f"and that --dead-barrel matches it.")

    contacts = new_residue_contacts(residues, dead_index, DEAD_PAIR_SIZE, args.clash_a)
    for distance, left, right in contacts[:5]:
        print(f"  CLASH {distance:.2f} A: {left} ... {right}")
    if contacts and not args.force_strained:
        raise SystemExit(
            f"refusing to write: the new residues make {len(contacts)} contact(s) below "
            f"{args.clash_a} A, the worst at {contacts[0][0]:.2f} A. Build the pair "
            f"from an AlphaFold model of the 510-aa TDX sequence (--from-af3) rather "
            f"than collapsing the chromophore in place; deleting a residue shortens a "
            f"barrel strand and the crystal structure does not know where the "
            f"shortened strand sits. --force-strained overrides.")

    write_pdb(output, residues)
    print(f"wrote {output}: {len(residues)} residues, "
          f"{sum(len(r['atoms']) for r in residues)} atoms")


if __name__ == "__main__":
    main()

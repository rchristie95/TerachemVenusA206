#!/usr/bin/env python3
"""Derive the dVenus-TDX sequence from the verified dVenus-TD deposit sequence.

dVenus-TDX is the single-chromophore control of Nguyen et al. (Biophys. J. 124,
4293-4309, 2025): the same tandem dimer, with the chromophore of ONE barrel
extinguished so the construct carries one absorber instead of two.  The
measured extinction coefficients are in exact 2:1 ratio (eps_516 = 184,400 for
TD against 92,200 for TDX), which is the whole content of "one chromophore
instead of two" -- the fold, the linker and the barrel-barrel geometry are
meant to be unchanged.

HOW the chromophore is extinguished was checked against the deposit on
2026-09-02 and is a DELETION, not a substitution.  Addgene #250939 carries a
510-aa insert that differs from the 511-aa #250938 parent by exactly one
opcode: the barrel-2 chromophore tyrosine Y339 is absent.  At the codon level
barrel 2 reads CTG GGC GGG CTG CAG (L-G-G-L-Q) against barrel 1's intact
CTC GGC TAC GGC CTG CAG, and Addgene's own metadata agrees -- insert 1549 bp
and vector 4450 bp against 1552 bp and 4453 bp, exactly three base pairs
shorter, which a point substitution cannot produce.  See TDX_CONSTRUCT.md.

So the Gly-Tyr-Gly chromophore triad becomes Gly-Gly.  Removing the tyrosine
removes the phenol that becomes the chromophore's anionic arm, so the triad can
never cyclise: the barrel still folds and still docks, but it holds no
chromophore and no -1 e.

    python3 alphafold_tandem/make_tdx_sequence.py                 # dark barrel 2
    python3 alphafold_tandem/make_tdx_sequence.py --dead-barrel 1 # dark barrel 1

Everything is therefore ONE character removed from the 511-aa Addgene #250938
sequence.  Nothing is transcribed by hand: the parent FASTA is read, the two
GYG triads are located, and the chosen tyrosine is dropped.  With
--dead-barrel 2 the result must equal the deposited insert byte for byte, and
this script checks that against pdVenus-TDX_Addgene250939_deposit_510aa.fasta
before writing.

`--mode substitute` reproduces the earlier Y->C reading of the construct, which
is Amber, the Vogel-lab non-fluorescent Venus (Koushik, Chen, Thaler, Puhl &
Vogel, Biophys. J. 91, L99, 2006; Addgene #27798).  Amber is a real and
different protein; it is simply not what #250939 contains, so it is not the
default and its output is named for what it is.

Writes the FASTA and an AlphaFold Server job request next to the parent FASTA.
The job request is the exact JSON the server's "upload JSON" button takes, and
reuses the model seed of the TD run in af3_server_run/ so the two are
comparable.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT_FASTA = HERE / "pdVenus-TD_Addgene250938_511aa.fasta"
DEPOSIT_FASTA = HERE / "pdVenus-TDX_Addgene250939_deposit_510aa.fasta"
TD_SEED = "425000176"  # seed of af3_server_run/fold_dvenus_td_tandem_job_request.json


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if line and not line.startswith(">"))


def chromophore_tyrosines(sequence: str) -> list[int]:
    """1-based positions of the tyrosine in each GYG chromophore-forming triad."""
    return [match.start() + 2 for match in re.finditer("GYG", sequence)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dead-barrel", type=int, choices=(1, 2), default=2,
                    help="Which barrel loses its chromophore (default: 2, the "
                         "C-terminal barrel, which is the one Addgene #250939 "
                         "darkens and which keeps the surviving chromophore at the "
                         "residue numbering it has in every TD artefact)")
    ap.add_argument("--mode", choices=("delete", "substitute"), default="delete",
                    help="delete: drop the chromophore tyrosine, which is what "
                         "Addgene #250939 does (510 aa, the default). substitute: "
                         "mutate it to cysteine, which is Amber and a different "
                         "protein (511 aa); provided so the earlier reading of the "
                         "construct stays reproducible, not because it is #250939.")
    ap.add_argument("--parent-fasta", type=Path, default=PARENT_FASTA)
    ap.add_argument("--deposit-fasta", type=Path, default=DEPOSIT_FASTA)
    ap.add_argument("--seed", default=TD_SEED)
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    td = read_fasta(args.parent_fasta)
    tyrosines = chromophore_tyrosines(td)
    if len(tyrosines) != 2:
        raise SystemExit(f"expected two GYG triads in {args.parent_fasta}, found {len(tyrosines)}")

    position = tyrosines[args.dead_barrel - 1]          # 1-based, in the 511-aa numbering
    if td[position - 1] != "Y":
        raise SystemExit(f"residue {position} is {td[position - 1]}, not Y")

    if args.mode == "delete":
        tdx = td[:position - 1] + td[position:]
        lesion = f"Y{position}del"
        detail = "chromophore tyrosine deleted, Gly-Tyr-Gly -> Gly-Gly"
    else:
        tdx = td[:position - 1] + "C" + td[position:]
        lesion = f"Y{position}C"
        detail = "chromophore tyrosine -> cysteine (Amber, Venus-local Y67C)"

    # Venus-local numbering of the same residue, so the lesion can also be named
    # the way the Vogel lab names its chromophore mutants: position 67.
    starts = [m.start() + 1 for m in re.finditer("MVSKGEELFTG", td)]
    if len(starts) != 2:
        raise SystemExit(f"expected two Venus barrels, found {len(starts)} start motifs")
    local = position - starts[args.dead_barrel - 1] + 1
    if local != 67:
        raise SystemExit(f"chromophore tyrosine is barrel-local {local}, expected 67")

    # The surviving triad must still be intact, and it must be the only one left.
    remaining = chromophore_tyrosines(tdx)
    expected_remaining = [tyrosines[0]] if args.dead_barrel == 2 else \
        [tyrosines[1] - (1 if args.mode == "delete" else 0)]
    if remaining != expected_remaining:
        raise SystemExit(f"expected one surviving GYG triad at {expected_remaining}, "
                         f"found {remaining}")

    # With the deposited lesion in the deposited barrel, this must BE the deposit.
    if args.mode == "delete" and args.dead_barrel == 2 and args.deposit_fasta.exists():
        deposit = read_fasta(args.deposit_fasta)
        if tdx != deposit:
            raise SystemExit(
                f"derived sequence does not match {args.deposit_fasta.name}; "
                f"lengths {len(tdx)} vs {len(deposit)}")
        print(f"checked against {args.deposit_fasta.name}: identical")

    length = len(tdx)
    name = f"dVenus-TDX_dark{args.dead_barrel}"
    if args.mode == "substitute":
        name += "_amberY67C"
    header = (f">{name}_derived_{length}aa  {lesion} ({detail}) "
              f"in barrel {args.dead_barrel} of Addgene #250938")

    stem = f"pdVenus-TDX_derived_{length}aa_dark{args.dead_barrel}"
    if args.mode == "substitute":
        stem += "_amberY67C"
    fasta_path = args.outdir / f"{stem}.fasta"
    lines = [header] + [tdx[i:i + 60] for i in range(0, length, 60)]
    fasta_path.write_text("\n".join(lines) + "\n")

    job = [{
        "name": name,
        "modelSeeds": [args.seed],
        "sequences": [{
            "proteinChain": {
                "sequence": tdx,
                "count": 1,
                "useStructureTemplate": True,
            }
        }],
        "dialect": "alphafoldserver",
        "version": 3,
    }]
    suffix = f"dark{args.dead_barrel}" + ("_amberY67C" if args.mode == "substitute" else "")
    json_path = args.outdir / f"af3_tdx_job_request_{suffix}.json"
    json_path.write_text(json.dumps(job, indent=1) + "\n")

    live = remaining[0]
    print(f"parent   {args.parent_fasta.name}: {len(td)} aa, GYG tyrosines at {tyrosines}")
    print(f"lesion   {lesion}  (barrel {args.dead_barrel}, Venus-local position 67)")
    print(f"         {detail}")
    print(f"live chromophore tyrosine stays at {live}")
    print(f"TDX      {length} aa ({length - len(td):+d} from TD)")
    print(f"wrote {fasta_path}")
    print(f"wrote {json_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Derive the dVenus-TDX sequence from the verified dVenus-TD deposit sequence.

dVenus-TDX is the single-chromophore control of Nguyen et al. (Biophys. J. 124,
4293-4309, 2025): the same tandem dimer, with the chromophore of ONE barrel
extinguished so the construct carries one absorber instead of two.  The
measured extinction coefficients are in exact 2:1 ratio (eps_516 = 184,400 for
TD against 92,200 for TDX), which is the whole content of "one chromophore
instead of two" -- the fold, the linker and the barrel-barrel geometry are
meant to be unchanged.

The dark barrel is Amber, the Vogel-lab non-fluorescent Venus, which is Venus
carrying the single point mutation Y67C in the Gly66-Tyr67-Gly68 chromophore
triad (Koushik, Chen, Thaler, Puhl & Vogel, Biophys. J. 91, L99, 2006;
Addgene #27798).  Removing the tyrosine removes the phenol that becomes the
chromophore's anionic arm, so the triad never cyclises: the barrel still folds
and still docks, but it holds no chromophore and no -1 e.

Everything here is therefore ONE character of the 511-aa Addgene #250938
sequence.  Nothing is transcribed by hand: the parent FASTA is read, the two
GYG triads are located, and the chosen tyrosine is replaced.

    python3 alphafold_tandem/make_tdx_sequence.py                 # dark barrel 2
    python3 alphafold_tandem/make_tdx_sequence.py --dead-barrel 1 # dark barrel 1

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
                    help="Which barrel carries the Amber mutation (default: 2, the "
                         "C-terminal barrel, so the surviving chromophore keeps the "
                         "residue numbering it has in every TD artefact)")
    ap.add_argument("--parent-fasta", type=Path, default=PARENT_FASTA)
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
    tdx = td[:position - 1] + "C" + td[position:]

    # Venus-local numbering of the same residue, so the mutation can also be
    # named the way the Vogel lab names it: Y67C.
    starts = [m.start() + 1 for m in re.finditer("MVSKGEELFTG", td)]
    if len(starts) != 2:
        raise SystemExit(f"expected two Venus barrels, found {len(starts)} start motifs")
    local = position - starts[args.dead_barrel - 1] + 1
    if local != 67:
        raise SystemExit(f"chromophore tyrosine is barrel-local {local}, expected 67")

    name = f"dVenus-TDX_dark{args.dead_barrel}"
    mutation = f"Y{position}C"
    header = (f">{name}_derived_511aa  {mutation} (Venus-local Y67C, Amber) "
              f"in barrel {args.dead_barrel} of Addgene #250938")

    fasta_path = args.outdir / f"pdVenus-TDX_derived_511aa_dark{args.dead_barrel}.fasta"
    lines = [header] + [tdx[i:i + 60] for i in range(0, len(tdx), 60)]
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
    json_path = args.outdir / f"af3_tdx_job_request_dark{args.dead_barrel}.json"
    json_path.write_text(json.dumps(job, indent=1) + "\n")

    live = tyrosines[0] if args.dead_barrel == 2 else tyrosines[1]
    print(f"parent   {args.parent_fasta.name}: {len(td)} aa, GYG tyrosines at {tyrosines}")
    print(f"mutation {mutation}  (barrel {args.dead_barrel}, Venus-local Y67C = Amber)")
    print(f"live chromophore tyrosine stays at {live}")
    print(f"TDX      {len(tdx)} aa, {sum(a != b for a, b in zip(td, tdx))} residue differs from TD")
    print(f"wrote {fasta_path}")
    print(f"wrote {json_path}")


if __name__ == "__main__":
    main()

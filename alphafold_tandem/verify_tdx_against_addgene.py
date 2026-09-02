#!/usr/bin/env python3
"""Compare the derived TDX sequences with the deposited Addgene #250939 insert.

The deposit was read on 2026-09-02 from the Addgene-verified NGS sequence of the
whole 4450 bp plasmid at https://www.addgene.org/250939/sequences/, its single
559-aa ORF translated, and the 49-aa pRSETb 6His leader split off.  Both halves
are stored beside this script so the check reruns without network access:

    pdVenus-TDX_Addgene250939_deposit_510aa.fasta
    pdVenus-TDX_Addgene250939_prsetb_leader_49aa.fasta

Run:  python3 alphafold_tandem/verify_tdx_against_addgene.py

Exit status 0 means the derived dark-barrel-2 sequence matches the deposit.
It does not, as of 2026-09-02: the deposit deletes the barrel-2 chromophore
tyrosine (510 aa) where the derivation substitutes it with cysteine (511 aa).
See TDX_CONSTRUCT.md, "What is derived and what is verified".
"""

import difflib
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


def read_fasta(name):
    return "".join(
        line.strip()
        for line in (HERE / name).read_text().splitlines()
        if not line.startswith(">")
    )


def report(label, ref, dep):
    ops = [
        op
        for op in difflib.SequenceMatcher(None, ref, dep, autojunk=False).get_opcodes()
        if op[0] != "equal"
    ]
    print("\n== Addgene #250939 insert vs %s (%d aa) ==" % (label, len(ref)))
    if not ops:
        print("   identical")
        return True
    for tag, i1, i2, j1, j2 in ops:
        print(
            "   %-7s ref %d-%d %r  ->  deposit %r"
            % (tag, i1 + 1, i2, ref[i1:i2], dep[j1:j2])
        )
    return False


def main():
    deposit = read_fasta("pdVenus-TDX_Addgene250939_deposit_510aa.fasta")
    td = read_fasta("pdVenus-TD_Addgene250938_511aa.fasta")
    dark2 = read_fasta("pdVenus-TDX_derived_511aa_dark2.fasta")
    dark1 = read_fasta("pdVenus-TDX_derived_511aa_dark1.fasta")

    print("deposit insert  %d aa" % len(deposit))
    print("TD  #250938     %d aa" % len(td))

    print("\nchromophore triad regions in the deposit:")
    print("   barrel 1, residues 64-73 : %s" % deposit[63:73])
    print("   barrel 2, residues 334-343: %s" % deposit[333:343])

    report("TD #250938", td, deposit)
    ok1 = report("derived TDX dark barrel 1", dark1, deposit)
    ok2 = report("derived TDX dark barrel 2", dark2, deposit)

    if ok1 or ok2:
        print("\nOK: a derived sequence matches the deposit.")
        return 0
    print(
        "\nMISMATCH: neither derived sequence matches Addgene #250939.\n"
        "See alphafold_tandem/TDX_CONSTRUCT.md before running the AlphaFold job\n"
        "or tandem_tdx_build/build_tdx_structure.py."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())

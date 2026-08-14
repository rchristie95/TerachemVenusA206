#!/usr/bin/env python3
"""Collect whatever the weekend queue finished. Safe to run at any time.

Reports the bright-state wavelength, oscillator strength and TRANSITION DIPOLE
for each basis. The dipole is the point of the exercise: the manuscript uses
|mu| = 9.6 D from the spectroscopically normalised STEOM density, while
extinction and Strickler-Berg give 7.5-7.9 D, and J scales as |mu|^2.

Unfinished jobs are reported as pending or failed, never guessed at.
"""
from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
AU_TO_DEBYE = 2.541746473

JOBS = [
    ("tddft_def2_svp",      "TDDFT wB97X-D3", "def2-SVP"),
    ("tddft_def2_svpd",     "TDDFT wB97X-D3", "def2-SVPD"),
    ("tddft_def2_tzvp",     "TDDFT wB97X-D3", "def2-TZVP"),
    ("tddft_def2_tzvpd",    "TDDFT wB97X-D3", "def2-TZVPD"),
    ("steom_def2_svp",      "DLPNO-STEOM",    "def2-SVP"),
    ("steom_def2_svpd",     "DLPNO-STEOM",    "def2-SVPD *production*"),
    ("steom_def2_tzvp",     "DLPNO-STEOM",    "def2-TZVP"),
    ("eomccsd_neutral_svp", "EOM-CCSD (neutral ladder)", "def2-SVP"),
]

ROW = re.compile(r"\s*\d+-\d+\w*\s+->\s+\d+-\d+\w*\s+([\d.]+)\s+([\d.]+)\s+"
                 r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)")


def bright(stem: str):
    p = HERE / f"{stem}.out"
    if not p.exists():
        return None, "pending"
    text = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        for marker, why in (("aborting the run", "aborted"),
                            ("Error", "error")):
            if marker in text:
                return None, why
        return None, "incomplete"
    # Last electric-dipole block; in a STEOM run that is the geometric-mean
    # transition moment, which is the correct non-Hermitian convention.
    blocks = text.split("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE")
    if len(blocks) < 2:
        return None, "no absorption block"
    tail = blocks[-1]
    cut = tail.find("SPECTRUM", 200)
    if cut > 0:
        tail = tail[:cut]
    best = None
    for line in tail.splitlines():
        m = ROW.match(line)
        if m:
            ev, cm, nm, f, d2, dx, dy, dz = (float(m.group(i)) for i in range(1, 9))
            mu_debye = (d2 ** 0.5) * AU_TO_DEBYE
            if best is None or f > best[3]:
                best = (ev, cm, nm, f, mu_debye)
    return (best, "ok") if best else (None, "no roots parsed")


def main():
    print("\nTransition-dipole basis convergence, 44-atom anion in the MM field")
    print("(the EOM-CCSD row is the 32-atom neutral ladder model)\n")
    print(f"{'method':<28}{'basis':<24}{'nm':>9}{'f':>8}{'|mu| D':>9}   status")
    print("-" * 88)
    for stem, method, basis in JOBS:
        b, why = bright(stem)
        if b:
            ev, cm, nm, f, mu = b
            print(f"{method:<28}{basis:<24}{nm:>9.2f}{f:>8.4f}{mu:>9.2f}   {why}")
        else:
            print(f"{method:<28}{basis:<24}{'--':>9}{'--':>8}{'--':>9}   {why}")
    print("\nreference: manuscript uses |mu| = 9.6 D (spec-normalised STEOM density);")
    print("           extinction + Strickler-Berg give 7.5-7.9 D.")
    print("           J scales as |mu|^2, so 9.6 -> 7.9 would move J by x0.68.")


if __name__ == "__main__":
    main()

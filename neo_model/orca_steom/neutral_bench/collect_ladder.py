#!/usr/bin/env python3
"""Collect every finished ladder rung into one table.

Safe to run at any time: rungs that have not finished are reported as pending
rather than guessed at. Nothing here is written into the manuscript
automatically -- the point is that every number can be traced to the output
file it came from.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
EV_NM = 1239.841984

ORCA_ROWS = [
    ("tddft_camb3lyp",      "TDDFT, CAM-B3LYP"),
    ("tddft_wb97xd3",       "TDDFT, wB97X-D3"),
    ("steom_neutral_check", "DLPNO-STEOM-CCSD"),
    ("eomccsd_neutral",     "EOM-CCSD (canonical)"),
]


def orca_bright(stem: str):
    """Brightest root from the last electric-dipole absorption block."""
    path = HERE / f"{stem}.out"
    if not path.exists():
        return None, "no output"
    text = path.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        return None, "did not terminate normally"
    blocks = text.split("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE")
    if len(blocks) < 2:
        return None, "no absorption block"
    # Stop at the next section header. ORCA prints VELOCITY DIPOLE and CD
    # SPECTRUM tables after this one with IDENTICAL row shape but different
    # column meanings -- the CD table's 4th column is a rotatory strength, which
    # can be tens or hundreds and is often negative. Scanning to end-of-file
    # therefore picks a rotatory strength and calls it an oscillator strength.
    tail = blocks[-1]
    cut = tail.find("SPECTRUM", 200)
    if cut > 0:
        tail = tail[:cut]
    best = None
    for line in tail.splitlines():
        m = re.match(r"\s*\d+-\d+\w*\s+->\s+\d+-\d+\w*\s+([\d.]+)\s+([\d.]+)\s+"
                     r"([\d.]+)\s+([\d.]+)", line)
        if m:
            ev, cm, nm, f = (float(m.group(i)) for i in range(1, 5))
            if best is None or f > best[3]:
                best = (ev, cm, nm, f)
    return (best, "ok") if best else (None, "no roots parsed")


def main():
    rows = []
    for stem, label in ORCA_ROWS:
        best, why = orca_bright(stem)
        rows.append((label, best, why))

    adc = HERE / "adc2_neutral.json"
    if adc.exists():
        d = json.loads(adc.read_text())["bright_state"]
        rows.append(("ADC(2)", (d["eV"], d["cm-1"], d["nm"], d["fosc"]), "ok"))
    else:
        rows.append(("ADC(2)", None, "pending"))

    print("\nNeutral gas-phase chromophore, def2-SVP, identical geometry")
    print(f"{'method':<24}{'lambda (nm)':>13}{'cm^-1':>11}{'f':>9}   status")
    print("-" * 70)
    for label, best, why in rows:
        if best:
            ev, cm, nm, f = best
            print(f"{label:<24}{nm:>13.2f}{cm:>11.1f}{f:>9.4f}   {why}")
        else:
            print(f"{label:<24}{'--':>13}{'--':>11}{'--':>9}   {why}")
    print("\nEOM-CCSD(fT) is absent by necessity: the Q-Chem licence expired")
    print("2026-07-25, ORCA has no excited-state triples correction, and CC3/CCSDT")
    print("are N^7/N^8 -- unreachable at 331 basis functions.")


if __name__ == "__main__":
    main()

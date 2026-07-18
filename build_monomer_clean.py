#!/usr/bin/env python3
"""
Stage A (offline): prepare the clean monomer PDB for the anionic tleap build.

From 1MYW chain A: apply the classic anionic-GFP pocket protonation
(His148 -> HID, Glu222 -> GLH), drop crystal waters (solvate fresh), convert the
CR2 chromophore HETATM -> ATOM, and write a single TER at the end so tleap keeps
the chromophore covalently bonded into the chain (no spurious mid-chain break).
"""
import os
from pathlib import Path

SRC = "1MYW.pdb"
OUT = "anionic_build/monomer_clean2.pdb"
Path("anionic_build").mkdir(exist_ok=True)

# Glu222 protonation is the variable under test: "GLH" (neutral, classic
# anionic-GFP B-state) or "GLU" (charged, as in the experiment-matching
# original run). Set via env var GLU222_STATE.
glu222 = os.environ.get("GLU222_STATE", "GLH").strip().upper()
assert glu222 in ("GLH", "GLU"), glu222

out = []
for l in open(SRC):
    if l[:6].strip() in ("ATOM", "HETATM") and l[21] == "A":
        rn = l[17:20].strip()
        if rn == "HOH":
            continue
        r = int(l[22:26])
        if rn == "HIS" and r == 148:
            l = l[:17] + "HID" + l[20:]
        if rn == "GLU" and r == 222 and glu222 == "GLH":
            l = l[:17] + "GLH" + l[20:]
        if l.startswith("HETATM"):
            l = "ATOM  " + l[6:]
        out.append(l)
out.append("TER\nEND\n")
Path(OUT).write_text("".join(out))
print(f"[A] wrote {OUT}: {sum(1 for x in out if x.startswith('ATOM'))} atoms "
      f"(HIS148->HID, GLU222->{glu222}, waters dropped, chromophore chained)")

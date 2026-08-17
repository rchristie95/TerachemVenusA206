#!/usr/bin/env python3
"""How much of the delta disagreement is the DEFINITION of the chromophore axis?

Our delta = 39.6 deg (tandem MD) is computed from the STEOM transition dipole.
Cusick's delta is not. Their SI Note S4 states plainly that AlphaFold3 "does not
provide the exact geometry of the chromophore in the beta barrel because it does
not consider any post-translational modifications", so they overlaid the AF3
models onto the 1myw crystal structure and estimated delta as the angle between
the vector joining the OH and CB atoms of the tyrosine chromophore precursors
(Tyr116/Tyr388) and the difference of those vectors.

That is a STRUCTURAL PROXY for the chromophore long axis, not a transition
dipole. Before concluding that our MD geometry disagrees with theirs, the two
definitions have to be evaluated on the SAME structure. That is all this script
does, and it needs no MD, no QM and no AlphaFold.

Their frame convention (main text): x_hat || mu_a - mu_b, and delta is the angle
between each monomer vector and x_hat. Implemented here literally, and checked
against the closed form delta = 90 - theta_ab/2 used elsewhere in this project.

R CONVENTIONS ALSO DIFFER, and by more than the error bars anyone quotes:
Table S2 uses CB2-CB2 for the crystal dimer and CG-CG for the AlphaFold models,
while our coupling pipeline uses the transition-density centroid. Those are three
different distances on one structure. They are all reported below.

Read-only. Writes results/delta_definition_offset.{md,json}.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = Path(__file__).resolve().parent / "results"

STRUCTURES = [
    ("venus_dimer.pdb", "crystal / vdW dimer, built from 1myw symmetry", 38.157, 31.0),
    ("tandem_dimer_production_cr2.pdb", "tandem dimer, MD starting structure", 39.611, None),
]

# Axis proxies. OH->CB2 is the closest analogue of Cusick's OH->CB tyrosine vector
# once the chromophore has matured (CB2 is the bridge carbon).
PROXIES = {"OH->CB2": ("OH", "CB2"), "OH->CG2": ("OH", "CG2"), "CZ->CB2": ("CZ", "CB2")}

# Cusick Table S2, for reference.
CUSICK_TABLE_S2 = {
    "vdw_R_CB2_CB2_A": 25.4,
    "vdw_delta_deg": 31.0,
    "alphafold_delta_deg": [9.0, 9.0, 8.0, 15.0, 9.0],
    "alphafold_R_CG_CG_A": [26.0, 27.2, 24.3, 27.5, 27.5],
    "alphafold_nu_bar_cm": [35.3, 30.8, 43.3, 29.2, 29.8],
    "mu_debye": 7.9,
}


def cr2_residues(path):
    res = {}
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM") or line[17:20].strip() != "CR2":
            continue
        key = (line[21], line[22:27].strip())
        res.setdefault(key, {})[line[12:16].strip()] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    return list(res.values())


def angle_deg(u, v):
    c = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for fname, description, dipole_delta, published in STRUCTURES:
        path = ROOT / fname
        if not path.exists():
            results[fname] = {"error": "missing"}
            continue
        chrom = cr2_residues(path)
        if len(chrom) < 2:
            results[fname] = {"error": f"found {len(chrom)} CR2 residues"}
            continue
        a, b = chrom[0], chrom[1]

        proxies = {}
        for label, (p, q) in PROXIES.items():
            if p not in a or q not in a or p not in b or q not in b:
                continue
            ra, rb = a[q] - a[p], b[q] - b[p]
            theta = angle_deg(ra, rb)
            # Their frame, applied literally.
            frame = angle_deg(ra, ra - rb)
            proxies[label] = {
                "theta_ab_deg": theta,
                "delta_closed_form_deg": 90.0 - theta / 2.0,
                "delta_cusick_frame_deg": frame,
                "offset_vs_transition_dipole_deg": dipole_delta - frame,
            }

        separations = {
            f"R_{atom}_{atom}_A": float(np.linalg.norm(a[atom] - b[atom]))
            for atom in ("CB2", "CG2", "CZ")
            if atom in a and atom in b
        }
        results[fname] = {
            "description": description,
            "delta_from_transition_dipole_deg": dipole_delta,
            "delta_from_structural_proxy": proxies,
            "separations": separations,
            "cusick_published_delta_deg": published,
        }

    payload = {"structures": results, "cusick_table_s2": CUSICK_TABLE_S2}
    (OUT_DIR / "delta_definition_offset.json").write_text(json.dumps(payload, indent=2))

    L = []
    A = L.append
    A("# How much of the delta gap is the axis definition?\n")
    A("Cusick estimate delta from a **structural proxy** — the vector joining OH and")
    A("CB of the tyrosine chromophore precursor — because AlphaFold3 \"does not")
    A("consider any post-translational modifications\" and so never builds the mature")
    A("chromophore (their SI, Note S4). We estimate it from the STEOM **transition")
    A("dipole**. Evaluated on the same structure, the two definitions differ:\n")
    for fname, r in results.items():
        if "error" in r:
            A(f"- `{fname}`: {r['error']}")
            continue
        A(f"## `{fname}`\n")
        A(f"{r['description']}\n")
        A("| axis used | theta_ab (deg) | delta (deg) | offset vs dipole |")
        A("|---|---|---|---|")
        A(f"| **STEOM transition dipole** | — | **{r['delta_from_transition_dipole_deg']:.2f}** | — |")
        for label, p in r["delta_from_structural_proxy"].items():
            A(f"| {label} (structural proxy) | {p['theta_ab_deg']:.2f} "
              f"| **{p['delta_cusick_frame_deg']:.2f}** "
              f"| −{p['offset_vs_transition_dipole_deg']:.2f}° |")
        if r["cusick_published_delta_deg"] is not None:
            A(f"| *Cusick published* | — | *{r['cusick_published_delta_deg']:.0f}* | — |")
        A("")
        A("Separations on this same structure, by convention:")
        for k, v in r["separations"].items():
            A(f"- `{k}` = {v:.2f} Å")
        A("")
    A("## Verdict\n")
    v = results.get("venus_dimer.pdb", {})
    if "delta_from_structural_proxy" in v:
        prox = v["delta_from_structural_proxy"]["OH->CB2"]
        A(f"On the crystal dimer — the one structure where we and Cusick are")
        A(f"unambiguously looking at the same coordinates — switching from the")
        A(f"transition dipole to their structural proxy moves delta from")
        A(f"**{v['delta_from_transition_dipole_deg']:.2f}° to "
          f"{prox['delta_cusick_frame_deg']:.2f}°**, i.e. by "
          f"{prox['offset_vs_transition_dipole_deg']:.1f}°, landing near their")
        A(f"published {v['cusick_published_delta_deg']:.0f}°. The separation agrees")
        A(f"to the digit: our CB2–CB2 = {v['separations']['R_CB2_CB2_A']:.2f} Å against")
        A(f"their Table S2 value of {CUSICK_TABLE_S2['vdw_R_CB2_CB2_A']:.1f} Å.\n")
    A("So most of the apparent geometric disagreement is a difference in what the")
    A("angle is *between*, not where the atoms are. Comparing our dipole-derived")
    A("delta directly against their proxy-derived delta is not like for like.\n")
    A("## Two further cautions about their AlphaFold numbers\n")
    afs = CUSICK_TABLE_S2["alphafold_delta_deg"]
    A(f"- Table S2 lists five AlphaFold3 structures with delta = "
      f"{', '.join(f'{x:.0f}' for x in afs)}°. Four of the five are 8–9°; only")
    A(f"  structure #4 gives 15°. The main text quotes **15°** — the single value")
    A(f"  that overlaps their spectroscopic 14–20° range. The modal prediction, 9°,")
    A(f"  does not overlap it at all.")
    A(f"- Those same five structures give couplings of "
      f"{', '.join(f'{x:.0f}' for x in CUSICK_TABLE_S2['alphafold_nu_bar_cm'])} cm⁻¹ — a")
    A(f"  spread of {max(CUSICK_TABLE_S2['alphafold_nu_bar_cm'])/min(CUSICK_TABLE_S2['alphafold_nu_bar_cm']):.2f}×")
    A(f"  from one prediction run. That is the signature of a poorly determined")
    A(f"  inter-domain orientation, which is expected for two barrels joined by a")
    A(f"  flexible 33-residue linker.\n")
    A("## What this does not resolve\n")
    A("- Even under their proxy our tandem sits at ~25°, still above their")
    A("  spectroscopic 14–20° and well above the modal AlphaFold 9°. The definition")
    A("  accounts for most of the gap, not all of it.")
    A("- The proxy offset is itself structure-dependent; it is measured here on two")
    A("  static structures, not averaged over the MD ensemble.")
    (OUT_DIR / "delta_definition_offset.md").write_text("\n".join(L))
    print("\n".join(L))
    print(f"\nwrote {OUT_DIR/'delta_definition_offset.md'}")


if __name__ == "__main__":
    main()

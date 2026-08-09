#!/usr/bin/env python3
"""Summarize the two-site smoke test without promoting it to production data."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from terachem_site_energy_cd.nondegenerate_spectra import frame_observables

EV_TO_CM = 8065.544005


def candidate(summary: dict) -> dict:
    candidates = []
    for root in summary["roots"]:
        ci = root["largest_excitation"]
        if ci["occupied"] == 93 and ci["virtual"] == 95:
            candidates.append(root)
    if not candidates:
        candidates = summary["roots"]
    return max(candidates, key=lambda item: item["oscillator_strength"])


def cr2_centroid(xyz: Path) -> np.ndarray:
    lines = xyz.read_text().splitlines()[2:31]
    return np.mean([[float(x) for x in line.split()[1:4]] for line in lines], axis=0)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("frame_dir", type=Path)
    p.add_argument("--coupling-csv", type=Path, required=True)
    p.add_argument("--frame", type=int, default=499)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    selected = {}
    for site in ("A", "B"):
        summary = json.loads((args.frame_dir / f"site_{site}" / "energy_summary.json").read_text())
        selected[site] = candidate(summary)
    with args.coupling_csv.open() as handle:
        rows = {int(row["frame"]): row for row in csv.DictReader(handle)}
    j_cm = float(rows[args.frame]["J_cm"])
    mu_a = np.asarray(selected["A"]["transition_dipole_au"], float)
    mu_b = np.asarray(selected["B"]["transition_dipole_au"], float)
    r_a = cr2_centroid(args.frame_dir / "site_A" / "geometry.xyz")
    r_b = cr2_centroid(args.frame_dir / "site_B" / "geometry.xyz")
    obs = frame_observables(
        selected["A"]["energy_cm-1"], selected["B"]["energy_cm-1"], j_cm,
        mu_a, mu_b, r_a, r_b,
    )
    payload = {
        "status": "smoke_converged_state_identity_ambiguous_production_join_forbidden",
        "frame": args.frame,
        "state_selection": {
            site: {
                "root": root["root"], "energy_eV": root["energy_eV"],
                "energy_cm-1": root["energy_cm-1"], "oscillator_strength": root["oscillator_strength"],
                "transition_dipole_au": root["transition_dipole_au"],
                "rule": "brightest root with dominant 93->95 configuration",
                "ambiguity": True,
                "reason": "No NTO/transition-density overlap was generated for the corrected-charge smoke input.",
            } for site, root in selected.items()
        },
        "detuning_cm-1": selected["A"]["energy_cm-1"] - selected["B"]["energy_cm-1"],
        "archived_J_cm_algorithmic_only": j_cm,
        "archived_J_join_valid": False,
        "join_warning": "The smoke DCD SHA-256 differs from the archived coupling trajectory; this J is used only to exercise invariants.",
        "hamiltonian": {
            "energies_cm": obs["energies_cm"].tolist(), "omega_cm": obs["omega_cm"],
            "mixing": obs["mixing"], "localization_weights": obs["localization_weights"].tolist(),
            "dipole_strengths": obs["dipole_strengths"].tolist(),
            "rotational_strengths_relative": obs["rotational_strengths_relative"].tolist(),
            "cd_units": "relative interaction-induced exciton chirality; not absolute molar CD",
        },
        "scientific_conclusion_allowed": False,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

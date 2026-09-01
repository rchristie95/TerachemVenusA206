#!/usr/bin/env python3
"""Do X-H constraints change the chromophore coordinate that carries the gap?

`constraints=HBonds` freezes the X-H stretches at ~3000 cm^-1, inside the band a
~12 fs T2 is sensitive to. This compares two otherwise identical 500 fs runs from
the same seed frame, one constrained and one not, on the METHINE BRIDGE LENGTH
and the bond-length alternation:

    bridge = d(CA2,CB2) + d(CB2,CG2)
    BLA    = d(CA2,CB2) - d(CB2,CG2)

using the same definitions as bridge_detuning_analysis.py:40-41. The v5 analysis
put 25.9% of the detuning variance on the bridge and 4.1% on electrostatics, so
this is the channel that matters. An electrostatic Stark proxy would be blind to
it -- that decomposition holds the chromophore geometry rigid and varies only the
MM field, so it cannot see a change in the chromophore's own internal motion.

Reports the sd, the ACF and its positive-sequence correlation time for each. If
the constrained and unconstrained runs agree, freezing X-H does not distort the
gap dynamics and the production T2 stands as computed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from openmm.app import PDBFile

from register_nvt_trajectories import dcd_frames
from summarize_100fs50 import acf_unbiased, positive_correlation_time

REPO = Path("/home/robson/PetaChem")
TOPOLOGY = REPO / "tc_tandem_v5_campaign/solvated_protonated.pdb"
NEEDED = ("CA2", "CB2", "CG2")


def cr2_atom_indices(topology) -> list[dict[str, int]]:
    sites = []
    for residue in topology.residues():
        if residue.name != "CR2":
            continue
        idx = {a.name: a.index for a in residue.atoms() if a.name in NEEDED}
        if set(idx) != set(NEEDED):
            raise RuntimeError(f"CR2 {residue.id} missing {set(NEEDED) - set(idx)}")
        sites.append(idx)
    if len(sites) != 2:
        raise RuntimeError(f"Expected 2 CR2 residues, found {len(sites)}")
    return sites


def series(dcd: Path, sites) -> dict:
    bridge, bla = [[], []], [[], []]
    for xyz in dcd_frames(dcd):
        for s, idx in enumerate(sites):
            a, b, c = (xyz[idx["CA2"]], xyz[idx["CB2"]], xyz[idx["CG2"]])
            d1 = float(np.linalg.norm(a - b))
            d2 = float(np.linalg.norm(b - c))
            bridge[s].append(d1 + d2)
            bla[s].append(d1 - d2)
    return {"bridge": np.asarray(bridge), "BLA": np.asarray(bla)}


def describe(name: str, values: np.ndarray, dt_fs: float) -> dict:
    """values: (2, nframes) -- one row per site."""
    out = {}
    for label, row in (("A", values[0]), ("B", values[1])):
        c = acf_unbiased(row)
        out[label] = {
            "mean_A": float(row.mean()),
            "sd_mA": float(row.std(ddof=1) * 1000.0),
            "tau_c_fs": positive_correlation_time(c, dt_fs),
        }
    # The gap-relevant quantity is the MISMATCH between the two sites, which is
    # what bridge_detuning_analysis correlates against the detuning.
    diff = values[0] - values[1]
    c = acf_unbiased(diff)
    out["mismatch"] = {
        "mean_mA": float(diff.mean() * 1000.0),
        "sd_mA": float(diff.std(ddof=1) * 1000.0),
        "tau_c_fs": positive_correlation_time(c, dt_fs),
        "acf_normalized": (c / c[0]).tolist()[:60],
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dt-fs", type=float, default=1.0)
    parser.add_argument("--output", type=Path,
                        default=REPO / "tandem_v6_t2_campaign/constraint_control.json")
    args = parser.parse_args()

    sites = cr2_atom_indices(PDBFile(str(TOPOLOGY)).topology)
    runs = {
        "constrained": REPO / "tc_tandem_v6_ctrl_constrained/ctrl_constrained.dcd",
        "unconstrained": REPO / "tc_tandem_v6_ctrl_unconstrained/ctrl_unconstrained.dcd",
    }
    result = {}
    for tag, dcd in runs.items():
        if not dcd.exists():
            raise FileNotFoundError(dcd)
        s = series(dcd, sites)
        result[tag] = {k: describe(k, v, args.dt_fs) for k, v in s.items()}
        result[tag]["n_frames"] = int(s["bridge"].shape[1])

    print("=== X-H constraint control: does freezing the stretches change the bridge? ===")
    print(f"    {result['constrained']['n_frames']} frames at {args.dt_fs:g} fs, "
          f"same seed frame (v5 900) and velocity seed\n")
    for coord in ("bridge", "BLA"):
        print(f"  {coord}:")
        print(f"    {'':16s} {'site A sd':>10s} {'site B sd':>10s} "
              f"{'mismatch sd':>12s} {'mismatch tau_c':>15s}")
        for tag in runs:
            r = result[tag][coord]
            print(f"    {tag:16s} {r['A']['sd_mA']:9.2f}m {r['B']['sd_mA']:9.2f}m "
                  f"{r['mismatch']['sd_mA']:11.2f}m {r['mismatch']['tau_c_fs']:14.2f}")
        cs = result["constrained"][coord]["mismatch"]["sd_mA"]
        us = result["unconstrained"][coord]["mismatch"]["sd_mA"]
        print(f"    -> mismatch sd ratio unconstrained/constrained = {us / cs:.3f}")
        print()

    args.output.write_text(json.dumps(result, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

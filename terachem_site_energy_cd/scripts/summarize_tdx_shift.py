#!/usr/bin/env python3
"""Paired TD-minus-TDX site-energy shift.

For every frame with both a completed TD calculation and its TDX control, this
reports E(TDX) - E(TD) per site.  A positive value means the presence of the
anionic partner chromophore lowers the site energy, i.e. red-shifts absorption,
which is the sign Nguyen and Kim measure (35.3 and 33.8 cm^-1 respectively for
the whole dVenus band).

The comparison is paired at fixed geometry: the two calculations differ only in
the 29 partner-CR2 point charges, so the per-frame difference is far better
determined than either absolute site energy, and the ensemble error is the
standard error of the differences rather than of the energies.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np


def bright_state(site_dir: Path, python: str) -> dict | None:
    """Highest-oscillator-strength root, with the runner-up margin."""
    summary = site_dir / "energy_summary.json"
    if not summary.is_file():
        parser = Path(__file__).with_name("parse_terachem_tddft.py")
        out = site_dir / "tddft.out"
        if not out.is_file() or "Job finished:" not in out.read_text(errors="replace"):
            return None
        subprocess.run(
            [python, str(parser), str(out), "--json", str(summary)],
            check=True, capture_output=True,
        )
    data = json.loads(summary.read_text())
    roots = [r for r in data.get("roots", []) if r.get("oscillator_strength") is not None]
    if len(roots) < 2:
        return None
    ranked = sorted(roots, key=lambda r: -r["oscillator_strength"])
    best, runner = ranked[0], ranked[1]
    return {
        "energy_cm": best["energy_cm-1"],
        "root": best["root"],
        "f": best["oscillator_strength"],
        "margin": best["oscillator_strength"] / max(runner["oscillator_strength"], 1e-12),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path,
                    default=Path("terachem_site_energy_cd/results"))
    ap.add_argument("--td-prefix", default="v2_camb3lyp_frame_")
    ap.add_argument("--tdx-prefix", default="tdx_neutral_frame_")
    ap.add_argument("--python", default="/home/robson/anaconda3/envs/TeraChem/bin/python")
    ap.add_argument("--min-margin", type=float, default=3.0)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    rows, dropped = [], []
    for tdx_dir in sorted(args.results_dir.glob(f"{args.tdx_prefix}*")):
        frame = int(re.search(r"(\d+)$", tdx_dir.name).group(1))
        td_dir = args.results_dir / f"{args.td_prefix}{frame:04d}"
        if not td_dir.is_dir():
            continue
        entry = {"frame": frame}
        usable = True
        for site in ("A", "B"):
            td = bright_state(td_dir / f"site_{site}", args.python)
            tdx = bright_state(tdx_dir / f"site_{site}", args.python)
            if td is None or tdx is None:
                usable = False
                break
            if min(td["margin"], tdx["margin"]) < args.min_margin:
                dropped.append((frame, site, round(min(td["margin"], tdx["margin"]), 2)))
                usable = False
                break
            entry[f"td_{site}"] = td["energy_cm"]
            entry[f"tdx_{site}"] = tdx["energy_cm"]
            entry[f"shift_{site}"] = tdx["energy_cm"] - td["energy_cm"]
        if usable:
            entry["shift_mean"] = 0.5 * (entry["shift_A"] + entry["shift_B"])
            rows.append(entry)

    if not rows:
        raise SystemExit("no paired TD/TDX frames found")

    shifts = np.array([r["shift_mean"] for r in rows])
    a = np.array([r["shift_A"] for r in rows])
    b = np.array([r["shift_B"] for r in rows])
    n = len(shifts)
    sem = shifts.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan")

    print(f"paired TD/TDX frames: {n}"
          + (f"  ({len(dropped)} dropped low-margin: {dropped})" if dropped else ""))
    print(f"  site A   E(TDX)-E(TD) = {a.mean():8.1f} +/- {a.std(ddof=1)/np.sqrt(n):.1f} cm^-1"
          if n > 1 else f"  site A   {a.mean():8.1f} cm^-1")
    print(f"  site B   E(TDX)-E(TD) = {b.mean():8.1f} +/- {b.std(ddof=1)/np.sqrt(n):.1f} cm^-1"
          if n > 1 else f"  site B   {b.mean():8.1f} cm^-1")
    print(f"  BAND     E(TDX)-E(TD) = {shifts.mean():8.1f} +/- {sem:.1f} cm^-1"
          f"   (sd {shifts.std(ddof=1):.1f})" if n > 1 else
          f"  BAND     {shifts.mean():8.1f} cm^-1")
    print(f"  measured: 35.3 (Nguyen dVenus-TD vs TDX), 33.8 (Kim free dimer)")

    if args.output:
        args.output.write_text(json.dumps(
            {"frames": rows, "dropped_low_margin": dropped,
             "mean_shift_cm": float(shifts.mean()),
             "sem_cm": float(sem) if n > 1 else None,
             "n": n}, indent=2) + "\n")


if __name__ == "__main__":
    main()

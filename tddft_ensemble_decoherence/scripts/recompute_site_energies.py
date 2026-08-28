#!/usr/bin/env python3
"""Recompute site energies from stored roots with the bright-moment estimator.

No new TDDFT: every frame's seven root energies and oscillator strengths are
already on disk. Site B shares its bright pi-pi* intensity with a dark state in
a third of the 4 fs frames, which makes a single-root pick ill-defined; the
oscillator-strength-weighted first moment recovers the diabatic energy exactly.
See audit/SITE_B_DIAGNOSIS.md for the seven causes tested and rejected.

    python3 scripts/recompute_site_energies.py --window short
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_ensemble import _bright_moment                     # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
EV_TO_CM = 8065.544


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", default="short")
    parser.add_argument("--timestep-fs", type=float, default=4.0)
    parser.add_argument("--window-cm", type=float, default=4000.0)
    parser.add_argument("--f-min", type=float, default=0.20)
    parser.add_argument("--out", type=Path, default=None)
    return parser.parse_args()


def collect(window: str, window_cm: float, f_min: float) -> tuple[list[int], dict]:
    base = ROOT / "frames" / window
    frames = sorted(int(path.name.split("_")[-1]) for path in base.glob("frame_*"))
    complete = [
        frame for frame in frames
        if all((base / f"frame_{frame:04d}" / f"site_{site}" / "density_summary.json").exists()
               for site in "AB")
    ]
    series: dict[str, dict[str, list]] = {}
    for site in "AB":
        moment, selected, n_bright = [], [], []
        for frame in complete:
            site_dir = base / f"frame_{frame:04d}" / f"site_{site}"
            roots = json.loads((site_dir / "density_summary.json").read_text())["roots"]
            table = {item["root"]: item for item in roots}
            moment.append(_bright_moment(table, window_cm, f_min) * EV_TO_CM)
            strengths = [item["oscillator_strength"] for item in roots]
            n_bright.append(sum(1 for value in strengths if value > 0.4))
            assignment = site_dir / "state_assignment.json"
            selected.append(
                table[json.loads(assignment.read_text())["selected_root"]]["energy_eV"] * EV_TO_CM
                if assignment.exists() else np.nan
            )
        series[site] = {
            "bright_moment_cm": np.array(moment),
            "selected_root_cm": np.array(selected),
            "n_bright_roots": np.array(n_bright),
        }
    return complete, series


def describe(values: np.ndarray, timestep_fs: float) -> dict:
    n = len(values)
    centred = values - values.mean()
    lag1 = float(np.dot(centred[:n - 1], centred[1:]) / np.dot(centred, centred))
    # residual spread after averaging over hbar/2J = 81 fs, the window on which
    # the exciton sees detuning as static
    block = max(1, int(round(81.0 / timestep_fs)))
    usable = n // block * block
    slow = float(values[:usable].reshape(-1, block).mean(1).std()) if usable >= 4 * block else float("nan")
    return {"mean": float(values.mean()), "sd": float(values.std()), "acf_lag1": lag1, "slow_sd": slow}


def main() -> None:
    args = parse_args()
    frames, series = collect(args.window, args.window_cm, args.f_min)
    print(f"window {args.window}: {len(frames)} complete frames "
          f"({(len(frames) - 1) * args.timestep_fs:.0f} fs)\n")
    print(f"{'site':>5} {'estimator':>16} {'mean':>9} {'sd':>8} {'ACF1':>6} {'slow(81fs)':>11}")
    stats = {}
    for site in "AB":
        for key, label in (("selected_root_cm", "selected root"), ("bright_moment_cm", "bright moment")):
            values = series[site][key]
            if np.isnan(values).any():
                continue
            summary = describe(values, args.timestep_fs)
            stats[(site, key)] = summary
            print(f"{site:>5} {label:>16} {summary['mean']:>9.1f} {summary['sd']:>8.1f} "
                  f"{summary['acf_lag1']:>6.2f} {summary['slow_sd']:>11.1f}")
    for key, label in (("selected_root_cm", "selected root"), ("bright_moment_cm", "bright moment")):
        if ("A", key) in stats and ("B", key) in stats:
            ratio = stats[("B", key)]["sd"] / stats[("A", key)]["sd"]
            gap = series["A"][key] - series["B"][key]
            print(f"\n{label}: sd_B/sd_A = {ratio:.2f}   "
                  f"gap sd {gap.std():.1f}   <|D|> {np.abs(gap).mean():.1f}")
    print("\nreference -- 1 ns production: sd 519/516, ratio 0.99, gap sd 686.5, <|D|> 576.4")
    out = args.out or ROOT / "audit" / f"site_energies_{args.window}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out,
        frame=np.array(frames),
        bright_moment_A_cm=series["A"]["bright_moment_cm"],
        bright_moment_B_cm=series["B"]["bright_moment_cm"],
        selected_root_A_cm=series["A"]["selected_root_cm"],
        selected_root_B_cm=series["B"]["selected_root_cm"],
        n_bright_A=series["A"]["n_bright_roots"],
        n_bright_B=series["B"]["n_bright_roots"],
        timestep_fs=args.timestep_fs,
        moment_window_cm=args.window_cm,
        moment_f_min=args.f_min,
    )
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

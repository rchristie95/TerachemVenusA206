#!/usr/bin/env python3
"""Morning readout for the overnight linker-release experiment.

Prints phase-resolved statistics of both arms and writes a four-panel figure.
The decision variable: after the bias is released, does the inter-CR2-axis
angle alpha settle back at the crystal-register value (~95-107 deg,
|cos| ~ 0.1-0.3) or move toward the anisotropy-required 131.3 deg
(|cos| = 0.66)? The control arm shows what unbiased sampling does alone.

Run with the TeraChem env python:
  /home/robson/anaconda3/envs/TeraChem/bin/python analyze_morning.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET_ALPHA = 131.3          # deg, from Nguyen limiting anisotropy
TARGET_ABS_COS = 0.660
CRYSTAL_ALPHA = 106.6         # deg, 1MYW register


def load(arm: str):
    path = HERE / arm / f"{arm}_metrics.csv"
    rows = list(csv.DictReader(open(path)))
    out = {}
    for key in rows[0]:
        if key == "phase":
            out[key] = np.array([r[key] for r in rows])
        else:
            out[key] = np.array([float(r[key]) for r in rows])
    return out


def block_stats(x, label, tail_frac=0.5):
    """Mean +/- SEM from 5 blocks over the last tail_frac of the data."""
    tail = x[int(len(x) * (1.0 - tail_frac)):]
    blocks = np.array_split(tail, 5)
    means = np.array([b.mean() for b in blocks])
    print(f"    {label:<22} {tail.mean():8.2f} +/- {means.std(ddof=1)/np.sqrt(5):5.2f} "
          f"(n={len(tail)}, last {tail_frac*100:.0f}%)")
    return tail.mean(), means.std(ddof=1) / np.sqrt(5)


def main():
    data = {}
    for arm in ("control", "release"):
        try:
            data[arm] = load(arm)
        except FileNotFoundError:
            print(f"[!] {arm}: no metrics file yet")
    if not data:
        return

    for arm, d in data.items():
        ns = d["time_ps"][-1] / 1000.0
        print(f"\n=== {arm}: {ns:.1f} ns sampled ===")
        if arm == "release":
            for phase in ("ramp", "hold", "released"):
                mask = d["phase"] == phase
                if not mask.any():
                    continue
                print(f"  phase {phase} ({mask.sum()} frames, "
                      f"{(d['time_ps'][mask][-1]-d['time_ps'][mask][0])/1000.0:.1f} ns):")
                block_stats(d["alpha_deg"][mask], "alpha (deg)")
                block_stats(np.abs(d["cos_alpha"][mask]), "|cos alpha|")
                block_stats(d["linker_e2e_ang"][mask], "linker e2e (A)")
                block_stats(d["cr2_sep_ang"][mask], "CR2 sep (A)")
                block_stats(d["triple_product_ang"][mask], "triple prod (A)")
        else:
            print(f"  unbiased ({len(d['alpha_deg'])} frames):")
            block_stats(d["alpha_deg"], "alpha (deg)")
            block_stats(np.abs(d["cos_alpha"]), "|cos alpha|")
            block_stats(d["linker_e2e_ang"], "linker e2e (A)")
            block_stats(d["cr2_sep_ang"], "CR2 sep (A)")
            block_stats(d["triple_product_ang"], "triple prod (A)")
        bad_t = np.abs(d["temperature_k"] - 300.0) > 15.0
        if bad_t.any():
            print(f"  [!] {bad_t.sum()} frames with |T-300| > 15 K")

    print(f"\nReference points: crystal alpha = {CRYSTAL_ALPHA} deg, "
          f"anisotropy requires alpha = {TARGET_ALPHA} deg (|cos| = {TARGET_ABS_COS})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
    colors = {"control": "tab:blue", "release": "tab:red"}
    for arm, d in data.items():
        t = d["time_ps"] / 1000.0
        axes[0].plot(t, d["alpha_deg"], lw=0.6, color=colors[arm], label=arm)
        axes[1].plot(t, np.abs(d["cos_alpha"]), lw=0.6, color=colors[arm])
        axes[2].plot(t, d["linker_e2e_ang"], lw=0.6, color=colors[arm])
        axes[3].plot(t, d["triple_product_ang"], lw=0.6, color=colors[arm])
    if "release" in data:
        d = data["release"]
        for phase, style in (("hold", ":"), ("released", "--")):
            mask = d["phase"] == phase
            if mask.any():
                t0 = d["time_ps"][mask][0] / 1000.0
                for ax in axes:
                    ax.axvline(t0, color="k", ls=style, lw=0.8, alpha=0.6)
    axes[0].axhline(TARGET_ALPHA, color="g", ls="--", lw=1, label="anisotropy 131.3")
    axes[0].axhline(CRYSTAL_ALPHA, color="gray", ls="--", lw=1, label="crystal 106.6")
    axes[0].set_ylabel("alpha (deg)")
    axes[0].legend(fontsize=8, ncol=4)
    axes[1].axhline(TARGET_ABS_COS, color="g", ls="--", lw=1)
    axes[1].set_ylabel("|cos alpha|")
    axes[2].axhline(38.0, color="g", ls="--", lw=1)
    axes[2].set_ylabel("linker e2e (A)")
    axes[3].axhline(0.0, color="gray", lw=0.8)
    axes[3].set_ylabel("R.(a x b) (A)")
    axes[3].set_xlabel("time (ns)")
    fig.suptitle("Overnight linker-release experiment (dotted: hold start, dashed: release)")
    fig.tight_layout()
    out = HERE / "overnight_summary.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure: {out}")


if __name__ == "__main__":
    main()

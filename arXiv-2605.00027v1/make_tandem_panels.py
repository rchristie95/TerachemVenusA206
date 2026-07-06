#!/usr/bin/env python3
"""
Generate individual panel PDFs for the tandem-dimer figure (Figure 4):
  (a) Chromophore separation vs frame
  (b) Coupling |J| vs frame
  (c) J histogram (from the spectra coupling panel data)

Reads: coupling_tandem_1000/coupling_samples.csv
       coupling_tandem_1000/coupling_distribution.json
Writes: Fig_Tandem_Separation.pdf, Fig_Tandem_Coupling.pdf, Fig_Tandem_Histogram.pdf
"""

import json
import sys
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).parent.parent / "coupling_tandem_1000"
OUT_DIR = Path(__file__).parent  # arXiv dir
CSV = DATA_DIR / "coupling_samples.csv"
JSON = DATA_DIR / "coupling_distribution.json"


def load_csv(path):
    import csv
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}


def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "legend.fontsize": 9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.formatter.use_mathtext": True,
    })
    return plt


C_ABS = "#1f6aa5"
C_MEAN = "#222222"
FIGSIZE = (4.2, 3.4)


def style(ax):
    ax.grid(alpha=0.25, lw=0.6)
    for s in ax.spines.values():
        s.set_linewidth(0.8)


def leg(ax, **kw):
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False, **kw)


def hist_bins(samples, binwidth):
    lo = np.floor(samples.min() / binwidth) * binwidth
    hi = np.ceil(samples.max() / binwidth) * binwidth + binwidth
    return np.arange(lo, hi, binwidth)


def main():
    data = load_csv(CSV)
    with open(JSON) as f:
        stats = json.load(f)

    plt = mpl()
    frames = data["frame"]
    sep = data["separation_A"]
    J = data["J_cm"]
    mean_J = stats["mean"]
    std_J = stats["std"]
    n = len(frames)

    # Time axis (1 ns over 1000 frames)
    time_ns = frames / frames.max()

    # --- Panel (a): separation vs time ---
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(time_ns, sep, color=C_ABS, lw=0.8, alpha=0.85)
    mean_sep = np.mean(sep)
    ax.axhline(mean_sep, color=C_MEAN, lw=1.8, ls="--",
               label=fr"mean $= {mean_sep:.1f}\,\mathrm{{\AA}}$")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(r"chromophore separation ($\mathrm{\AA}$)")
    ax.set_ylim(22, 28)
    style(ax)
    leg(ax, loc="upper right")
    fig.tight_layout()
    p_a = OUT_DIR / "Fig_Tandem_Separation.pdf"
    fig.savefig(p_a)
    plt.close(fig)
    print(f"  wrote {p_a}")

    # --- Panel (b): coupling vs time ---
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(time_ns, J, color=C_ABS, lw=0.8, alpha=0.85)
    ax.axhline(mean_J, color=C_MEAN, lw=1.8, ls="--",
               label=fr"$\bar J = {mean_J:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.axhspan(mean_J - std_J, mean_J + std_J, color=C_MEAN, alpha=0.10,
               label=fr"$\pm\sigma = {std_J:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(r"coupling $|J|$ (cm$^{-1}$)")
    style(ax)
    leg(ax, loc="upper right")
    fig.tight_layout()
    p_b = OUT_DIR / "Fig_Tandem_Coupling.pdf"
    fig.savefig(p_b)
    plt.close(fig)
    print(f"  wrote {p_b}")

    # --- Panel (c): coupling histogram ---
    fig, ax = plt.subplots(figsize=FIGSIZE)
    bw = 0.4
    ax.hist(J, bins=hist_bins(J, bw), color=C_ABS, alpha=0.80,
            edgecolor="white", lw=0.4)
    ax.axvline(mean_J, color=C_MEAN, lw=1.8,
               label=fr"$\bar J = {mean_J:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.axvspan(mean_J - std_J, mean_J + std_J, color=C_MEAN, alpha=0.10,
               label=fr"$\pm\sigma = {std_J:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.set_xlabel(r"Davydov coupling $J$ (cm$^{-1}$)")
    ax.set_ylabel("MD snapshots")
    style(ax)
    leg(ax, loc="upper right")
    fig.tight_layout()
    p_c = OUT_DIR / "Fig_Tandem_Histogram.pdf"
    fig.savefig(p_c)
    plt.close(fig)
    print(f"  wrote {p_c}")

    print(f"\n  Panel (a): {p_a.name}  separation vs time")
    print(f"  Panel (b): {p_b.name}  coupling vs time")
    print(f"  Panel (c): {p_c.name}  coupling histogram")


if __name__ == "__main__":
    main()

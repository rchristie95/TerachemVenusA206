#!/usr/bin/env python3
"""Generate the three production-NVT panels used in manuscript Figure 5.

The coupling and geometry statistics are calculated from ``coupling_samples.csv``
at run time; no manuscript result is hard-coded.  By default the script reads the
1000-frame, 1 ns production calculation and writes PDF and PNG versions of:

* ``Fig_Tandem_Separation``: chromophore separation versus time;
* ``Fig_Tandem_Coupling``: mutual transition-density coupling versus time;
* ``Fig_Tandem_Histogram``: the sampled coupling distribution.

Use ``--data-dir`` to regenerate the panels from another audited coupling run.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = HERE.parent / "coupling_nvt_production_cr2_1000_20260721"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="directory containing coupling_samples.csv (default: production 1000-frame run)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=HERE,
        help="output directory for panel PDFs/PNGs (default: manuscript directory)",
    )
    parser.add_argument(
        "--duration-ns",
        type=float,
        default=1.0,
        help="trajectory duration represented by all saved frames (default: 1.0 ns)",
    )
    return parser.parse_args()


def load_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No samples found in {path}")
    required = {"frame", "J_cm", "separation_A"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def validate_samples(data: dict[str, np.ndarray], csv_path: Path) -> None:
    lengths = {len(values) for values in data.values()}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent column lengths in {csv_path}")
    for key in ("frame", "J_cm", "separation_A"):
        if not np.all(np.isfinite(data[key])):
            raise ValueError(f"Non-finite {key} values in {csv_path}")
    if np.any(np.diff(data["frame"]) <= 0):
        raise ValueError(f"Frame numbers must be strictly increasing in {csv_path}")


def mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
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
        }
    )
    return plt


C_ABS = "#1f6aa5"
C_MEAN = "#222222"
FIGSIZE = (4.2, 3.4)


def style(ax) -> None:
    ax.grid(alpha=0.25, lw=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def leg(ax, **kwargs) -> None:
    ax.legend(
        frameon=True,
        framealpha=0.95,
        edgecolor="0.75",
        fancybox=False,
        **kwargs,
    )


def histogram_bins(samples: np.ndarray) -> np.ndarray:
    """Return stable Freedman-Diaconis bins with sensible small-sample fallbacks."""
    lo = float(samples.min())
    hi = float(samples.max())
    if np.isclose(lo, hi):
        return np.asarray([lo - 0.5, hi + 0.5])
    q25, q75 = np.percentile(samples, [25.0, 75.0])
    width = 2.0 * (q75 - q25) / np.cbrt(samples.size)
    if not np.isfinite(width) or width <= 0.0:
        width = (hi - lo) / max(10, int(np.sqrt(samples.size)))
    nbins = max(5, int(np.ceil((hi - lo) / width)))
    return np.linspace(lo, hi, nbins + 1)


def save_panel(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")


def check_summary_json(data_dir: Path, mean_j: float, std_j: float, n: int) -> None:
    """Warn if a stale summary JSON disagrees with the samples used for the panels."""
    path = data_dir / "coupling_distribution.json"
    if not path.exists():
        return
    with path.open(encoding="utf-8") as handle:
        summary = json.load(handle)
    mismatches = []
    if int(summary.get("n", n)) != n:
        mismatches.append("n")
    if not np.isclose(float(summary.get("mean", mean_j)), mean_j, atol=1e-8):
        mismatches.append("mean")
    if not np.isclose(float(summary.get("std", std_j)), std_j, atol=1e-8):
        mismatches.append("std")
    if mismatches:
        print(f"[warning] {path.name} disagrees with CSV fields: {', '.join(mismatches)}")


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    csv_path = data_dir / "coupling_samples.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Production coupling samples are not available yet: {csv_path}\n"
            "Run again after the 1000-frame GPU coupling calculation has completed, "
            "or supply --data-dir."
        )
    if args.duration_ns <= 0.0:
        raise ValueError("--duration-ns must be positive")

    data = load_csv(csv_path)
    validate_samples(data, csv_path)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    frames = data["frame"]
    separation = data["separation_A"]
    coupling = data["J_cm"]
    n = len(frames)
    mean_j = float(np.mean(coupling))
    # The coupling audit JSON reports the usual sample standard deviation.
    std_j = float(np.std(coupling, ddof=1)) if n > 1 else 0.0
    mean_sep = float(np.mean(separation))
    std_sep = float(np.std(separation, ddof=1)) if n > 1 else 0.0
    check_summary_json(data_dir, mean_j, std_j, n)

    # Frame zero is the first saved post-integration snapshot, not t = 0.
    time_ns = np.arange(1, n + 1, dtype=float) * (args.duration_ns / n)
    plt = mpl()

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(time_ns, separation, color=C_ABS, lw=0.8, alpha=0.85)
    ax.axhline(
        mean_sep,
        color=C_MEAN,
        lw=1.8,
        ls="--",
        label=fr"mean $= {mean_sep:.2f}\,\mathrm{{\AA}}$",
    )
    ax.axhspan(
        mean_sep - std_sep,
        mean_sep + std_sep,
        color=C_MEAN,
        alpha=0.08,
        label=fr"$\pm\sigma = {std_sep:.2f}\,\mathrm{{\AA}}$",
    )
    ax.set_xlim(0.0, args.duration_ns)
    margin = max(0.25, 0.08 * np.ptp(separation))
    ax.set_ylim(float(separation.min() - margin), float(separation.max() + margin))
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(r"chromophore separation ($\mathrm{\AA}$)")
    style(ax)
    leg(ax, loc="best")
    fig.tight_layout()
    separation_stem = out_dir / "Fig_Tandem_Separation"
    save_panel(fig, separation_stem)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(time_ns, coupling, color=C_ABS, lw=0.8, alpha=0.85)
    ax.axhline(
        mean_j,
        color=C_MEAN,
        lw=1.8,
        ls="--",
        label=fr"$\bar J = {mean_j:.1f}\,\mathrm{{cm^{{-1}}}}$",
    )
    ax.axhspan(
        mean_j - std_j,
        mean_j + std_j,
        color=C_MEAN,
        alpha=0.10,
        label=fr"$\pm\sigma = {std_j:.1f}\,\mathrm{{cm^{{-1}}}}$",
    )
    ax.set_xlim(0.0, args.duration_ns)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(r"TDC coupling $J$ (cm$^{-1}$)")
    style(ax)
    leg(ax, loc="best")
    fig.tight_layout()
    coupling_stem = out_dir / "Fig_Tandem_Coupling"
    save_panel(fig, coupling_stem)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.hist(
        coupling,
        bins=histogram_bins(coupling),
        color=C_ABS,
        alpha=0.80,
        edgecolor="white",
        lw=0.4,
    )
    ax.axvline(
        mean_j,
        color=C_MEAN,
        lw=1.8,
        label=fr"$\bar J = {mean_j:.1f}\,\mathrm{{cm^{{-1}}}}$",
    )
    ax.axvspan(
        mean_j - std_j,
        mean_j + std_j,
        color=C_MEAN,
        alpha=0.10,
        label=fr"$\pm\sigma = {std_j:.1f}\,\mathrm{{cm^{{-1}}}}$",
    )
    ax.set_xlabel(r"TDC coupling $J$ (cm$^{-1}$)")
    ax.set_ylabel("MD snapshots")
    style(ax)
    leg(ax, loc="best")
    fig.tight_layout()
    histogram_stem = out_dir / "Fig_Tandem_Histogram"
    save_panel(fig, histogram_stem)
    plt.close(fig)

    print(f"source: {csv_path}")
    print(f"frames: {n}; represented duration: {args.duration_ns:.6g} ns")
    print(f"J: {mean_j:.6f} +/- {std_j:.6f} cm^-1 (sample SD)")
    print(f"separation: {mean_sep:.6f} +/- {std_sep:.6f} A (sample SD)")
    for stem in (separation_stem, coupling_stem, histogram_stem):
        print(f"wrote {stem.with_suffix('.pdf').name} and {stem.with_suffix('.png').name}")


if __name__ == "__main__":
    main()

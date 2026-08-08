#!/usr/bin/env python3
"""Plot the 24 fs and 60 fs numerical CD spectra on one amplitude scale.

The manuscript spectra are normalized independently for shape comparison.  That
normalization hides the loss of raw model amplitude caused by faster homogeneous
dephasing.  This companion plot divides both raw spectra by the 60 fs peak so
that linewidth, extrema motion, and amplitude loss can be inspected together.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANUSCRIPT = ROOT / "manuscript"
sys.path.insert(0, str(MANUSCRIPT))

from make_nguyen_style_spectra import (  # noqa: E402
    homogeneous_hwhm_cm,
    load_ensemble,
    load_site_energy,
    numerical_cd_components,
    sha256,
    wavelength_order,
)

DEFAULT_DATA = ROOT / "coupling_nvt_production_cr2_1000_20260721"
DEFAULT_VALIDATION = ROOT / "reference" / "orca_validation.json"
DEFAULT_PME = ROOT / "solvation_decoherence_test" / "pme_validation_8ps" / "summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--validation-json", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--pme-summary", type=Path, default=DEFAULT_PME)
    parser.add_argument("--reference-t2-fs", type=float, default=60.0)
    parser.add_argument("--npts", type=int, default=8000)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    return parser.parse_args()


def difference_spectrum(
    grid_cm: np.ndarray,
    ensemble: dict[str, np.ndarray],
    e0_cm: float,
    t2_fs: float,
) -> tuple[np.ndarray, float]:
    gamma_cm = homogeneous_hwhm_cm(t2_fs)
    td_high, td_low, uncoupled = numerical_cd_components(
        grid_cm, ensemble, e0_cm, gamma_cm
    )
    return uncoupled - (td_high + td_low), gamma_cm


def extremum(wavelength: np.ndarray, signal: np.ndarray, maximum: bool) -> tuple[float, float]:
    index = int(np.argmax(signal) if maximum else np.argmin(signal))
    return float(wavelength[index]), float(signal[index])


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return f"external/{resolved.name}"


def main() -> None:
    args = parse_args()
    if args.reference_t2_fs <= 0.0 or args.npts < 2000:
        raise ValueError("Require --reference-t2-fs > 0 and --npts >= 2000")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    data_dir = args.data_dir.resolve()
    validation_path = args.validation_json.resolve()
    pme_path = args.pme_summary.resolve()
    ensemble, csv_source, npz_source = load_ensemble(data_dir)
    e0_cm, _, _ = load_site_energy(validation_path)
    with pme_path.open(encoding="utf-8") as handle:
        pme_record = json.load(handle)
    t2_fast = float(pme_record["PME"]["classical_cumulant_T2_1e_fs"])

    grid_cm = np.linspace(1.0e7 / 600.0, 1.0e7 / 400.0, args.npts)
    raw60_grid, hwhm60 = difference_spectrum(
        grid_cm, ensemble, e0_cm, args.reference_t2_fs
    )
    raw24_grid, hwhm24 = difference_spectrum(
        grid_cm, ensemble, e0_cm, t2_fast
    )
    wavelength, order = wavelength_order(grid_cm)
    wl60, raw60 = wavelength, raw60_grid[order]
    wl24, raw24 = wavelength, raw24_grid[order]

    scale = float(np.max(np.abs(raw60)))
    if scale <= 0.0:
        raise ValueError("The 60 fs reference spectrum has zero amplitude")
    y60 = raw60 / scale
    y24 = raw24 / scale
    peak_ratio = float(np.max(np.abs(y24)))
    reduction = 1.0 - peak_ratio

    mean_j = float(np.mean(ensemble["J_cm"]))
    short_nm = 1.0e7 / (e0_cm + mean_j)
    long_nm = 1.0e7 / (e0_cm - mean_j)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9.5,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(
        wl60,
        y60,
        color="#2468a2",
        lw=2.2,
        label=rf"$T_2^*=60$ fs (HWHM {hwhm60:.1f} cm$^{{-1}}$)",
    )
    ax.plot(
        wl24,
        y24,
        color="#c43c39",
        lw=2.2,
        label=rf"$T_2^*={t2_fast:.2f}$ fs (HWHM {hwhm24:.1f} cm$^{{-1}}$)",
    )
    ax.axhline(0.0, color="0.25", lw=0.8)
    for centre in (short_nm, long_nm):
        ax.axvline(centre, color="0.55", lw=0.9, ls=":")

    ax.set_xlim(500.0, 545.0)
    ax.set_ylim(-1.08, 1.08)
    ax.set_xlabel("Photon wavelength (nm)")
    ax.set_ylabel("Interaction-induced CD (common scale)")
    ax.set_title(
        "Homogeneous-dephasing sensitivity\n"
        rf"Corrected production coupling: $J={mean_j:.2f}$ cm$^{{-1}}$",
        fontsize=12,
    )
    ax.grid(alpha=0.18, lw=0.6)
    ax.legend(loc="upper left", frameon=False)
    ax.text(
        0.985,
        0.055,
        f"24 fs peak = {peak_ratio:.3f} of 60 fs\n"
        f"raw peak reduction = {100.0 * reduction:.1f}%",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.94},
    )

    pdf_path = args.out_dir / "Fig_T2_CommonScale.pdf"
    png_path = args.out_dir / "Fig_T2_CommonScale.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    min60 = extremum(wl60, raw60, maximum=False)
    max60 = extremum(wl60, raw60, maximum=True)
    min24 = extremum(wl24, raw24, maximum=False)
    max24 = extremum(wl24, raw24, maximum=True)
    record = {
        "description": "Common-amplitude-scale comparison of the numerical 60 fs and 24 fs tandem CD spectra",
        "normalization": "Both raw difference spectra divided by the absolute 60 fs peak",
        "T2_60_fs": args.reference_t2_fs,
        "T2_24_fs": t2_fast,
        "hwhm_60_cm-1": hwhm60,
        "hwhm_24_cm-1": hwhm24,
        "raw_peak_60": scale,
        "raw_peak_24": float(np.max(np.abs(raw24))),
        "peak_ratio_24_over_60": peak_ratio,
        "peak_reduction_fraction": reduction,
        "coupling_unit_status": "corrected reciprocal-distance conversion",
        "mean_J_cm-1": mean_j,
        "mean_splitting_cm-1": 2.0 * abs(mean_j),
        "extrema_60_nm": {"negative": min60[0], "positive": max60[0]},
        "extrema_24_nm": {"negative": min24[0], "positive": max24[0]},
        "exciton_centres_nm": {"short": short_nm, "long": long_nm},
        "sources": {
            "coupling_csv": portable_path(csv_source),
            "coupling_csv_sha256": sha256(csv_source),
            "geometry_npz": portable_path(npz_source),
            "geometry_npz_sha256": sha256(npz_source),
            "orca_validation": portable_path(validation_path),
            "orca_validation_sha256": sha256(validation_path),
            "pme_summary": portable_path(pme_path),
            "pme_summary_sha256": sha256(pme_path),
        },
        "outputs": {"pdf": portable_path(pdf_path), "png": portable_path(png_path)},
    }
    with (args.out_dir / "decoherence_note_figure_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")

    print(f"[note-figure] wrote {pdf_path}")
    print(f"[note-figure] wrote {png_path}")
    print(f"[note-figure] 24/60 raw peak ratio = {peak_ratio:.6f}")


if __name__ == "__main__":
    main()

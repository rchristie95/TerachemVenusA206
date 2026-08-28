#!/usr/bin/env python3
"""Create compact statistics for the corrected 50-frame 100 ns campaign."""

from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "manifests/long100ns50_frames.csv"
OUTDIR = ROOT / "rapid_100ns50"
HC_EV_NM = 1239.8419843320026


def correlation_block_length(values: np.ndarray) -> tuple[float, int]:
    centered = values - np.mean(values)
    variance = float(np.dot(centered, centered))
    if len(values) < 4 or variance == 0.0:
        return 0.5, 1
    acf = np.correlate(centered, centered, mode="full")[len(values) - 1 :] / variance
    positive = []
    for value in acf[1:]:
        if value <= 0.0:
            break
        positive.append(float(value))
    tau = 0.5 + sum(positive)
    return tau, max(1, min(len(values) // 4, int(np.ceil(2.0 * tau))))


def block_bootstrap_mean_ci(
    values: np.ndarray, block: int, rng: np.random.Generator, draws: int = 5000
) -> list[float]:
    n = len(values)
    starts = np.arange(n)
    means = np.empty(draws)
    blocks_needed = int(np.ceil(n / block))
    offsets = np.arange(block)
    for draw in range(draws):
        chosen = rng.choice(starts, blocks_needed, replace=True)
        indices = (chosen[:, None] + offsets[None, :]).ravel()[:n] % n
        means[draw] = float(np.mean(values[indices]))
    return [float(value) for value in np.quantile(means, [0.025, 0.975])]


def stats(values: np.ndarray, rng: np.random.Generator) -> dict[str, float | int | list[float]]:
    tau, block = correlation_block_length(values)
    return {
        "n": int(len(values)),
        "mean": float(np.mean(values)),
        "sample_sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "median": float(np.median(values)),
        "q025": float(np.quantile(values, 0.025)),
        "q05": float(np.quantile(values, 0.05)),
        "q25": float(np.quantile(values, 0.25)),
        "q75": float(np.quantile(values, 0.75)),
        "q95": float(np.quantile(values, 0.95)),
        "q975": float(np.quantile(values, 0.975)),
        "integrated_autocorrelation_time_samples_positive_sequence": tau,
        "bootstrap_block_length_samples": block,
        "block_bootstrap_mean_ci95": block_bootstrap_mean_ci(values, block, rng),
    }


def main() -> None:
    calibration = json.loads((OUTDIR / "calibration.json").read_text())
    delta = float(calibration["steom_constant_shift_eV"])
    rng = np.random.default_rng(20260723)
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    complete = [row for row in rows if row["convergence_status"] == "complete"]
    failed = [row for row in rows if row["convergence_status"] == "failed"]
    if len(complete) + len(failed) != 50:
        raise RuntimeError("Campaign still has unclassified frames")
    if not complete:
        raise RuntimeError("No completed frames")

    ea = np.asarray([float(row["site_A_energy_eV"]) for row in complete])
    eb = np.asarray([float(row["site_B_energy_eV"]) for row in complete])
    eac = np.asarray([float(row["site_A_corrected_eV"]) for row in complete])
    ebc = np.asarray([float(row["site_B_corrected_eV"]) for row in complete])
    jvac = np.asarray([float(row["J_vacuum_cm-1"]) for row in complete])
    jscr = np.asarray([float(row["J_screened_cm-1"]) for row in complete])
    jpda = np.asarray([float(row["J_pda_screened_cm-1"]) for row in complete])
    sep = np.asarray([float(row["separation_A"]) for row in complete])
    gap = ea - eb
    corrected_gap = eac - ebc
    mean_site = (ea + eb) / 2.0
    abs_j = np.abs(jscr)
    splitting = 2.0 * abs_j
    gap_invariance = float(np.max(np.abs(gap - corrected_gap)))
    if gap_invariance > 2.0e-9:
        raise RuntimeError(f"Constant-shift gap invariance failed: {gap_invariance}")

    result = {
        "status": "complete" if not failed else "complete_with_classified_failures",
        "requested_frames": 50,
        "complete_frames": len(complete),
        "failed_frames": len(failed),
        "failed_frame_indices": [int(row["frame_index"]) for row in failed],
        "failure_classes": dict(Counter(row["failure_class"] for row in failed)),
        "coverage_ns": [2.0, 100.0],
        "sampling_interval_ns": 2.0,
        "method": calibration["method"],
        "basis": calibration["basis"],
        "steom_constant_shift_eV": delta,
        "constant_shift_gap_invariance_max_abs_eV": gap_invariance,
        "distributions": {
            "site_A_raw_eV": stats(ea, rng),
            "site_B_raw_eV": stats(eb, rng),
            "site_A_raw_nm": stats(HC_EV_NM / ea, rng),
            "site_B_raw_nm": stats(HC_EV_NM / eb, rng),
            "site_A_STEOM_shifted_eV": stats(eac, rng),
            "site_B_STEOM_shifted_eV": stats(ebc, rng),
            "site_A_STEOM_shifted_nm": stats(HC_EV_NM / eac, rng),
            "site_B_STEOM_shifted_nm": stats(HC_EV_NM / ebc, rng),
            "mean_site_raw_eV": stats(mean_site, rng),
            "differential_gap_eV": stats(gap, rng),
            "J_vacuum_cm-1": stats(jvac, rng),
            "J_screened_cm-1": stats(jscr, rng),
            "absolute_J_screened_cm-1": stats(abs_j, rng),
            "Davydov_splitting_screened_cm-1": stats(splitting, rng),
            "J_PDA_screened_cm-1": stats(jpda, rng),
            "separation_A": stats(sep, rng),
        },
        "selected_root_counts": {
            "site_A": dict(Counter(row["site_A_selected_root"] for row in complete)),
            "site_B": dict(Counter(row["site_B_selected_root"] for row in complete)),
        },
        "interpretation_limit": (
            "The 2 ns sampling interval characterizes slow inhomogeneous disorder; "
            "it cannot resolve femtosecond decoherence or determine T2*."
        ),
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MANIFEST, OUTDIR / "frame_manifest.csv")
    (OUTDIR / "summary.json").write_text(json.dumps(result, indent=2) + "\n")

    figure, axes = plt.subplots(2, 3, figsize=(13.5, 7.5))
    bins_energy = np.histogram_bin_edges(np.concatenate([ea, eb]), bins="auto")
    axes[0, 0].hist(ea, bins=bins_energy, alpha=0.62, label="site A")
    axes[0, 0].hist(eb, bins=bins_energy, alpha=0.62, label="site B")
    axes[0, 0].set(xlabel="Raw excitation energy (eV)", ylabel="Frames")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].hist(mean_site, bins="auto", color="#4c78a8")
    axes[0, 1].set(xlabel="Mean site energy (eV)", ylabel="Frames")
    axes[0, 2].hist(gap, bins="auto", color="#f58518")
    axes[0, 2].set(xlabel=r"$E_A-E_B$ (eV)", ylabel="Frames")
    bins_j = np.histogram_bin_edges(jscr, bins="auto")
    axes[1, 0].hist(jscr, bins=bins_j, color="#54a24b")
    axes[1, 0].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 0].set(xlabel=r"Signed screened $J$ (cm$^{-1}$)", ylabel="Frames")
    axes[1, 1].hist(abs_j, bins="auto", color="#e45756")
    axes[1, 1].set(xlabel=r"$|J|$ (cm$^{-1}$)", ylabel="Frames")
    axes[1, 2].hist(splitting, bins="auto", color="#b279a2")
    axes[1, 2].set(xlabel=r"Davydov splitting $2|J|$ (cm$^{-1}$)", ylabel="Frames")
    figure.suptitle("Rapid TDDFT/TDC ensemble: 50 frames across 100 ns")
    figure.tight_layout()
    figure.savefig(OUTDIR / "energy_coupling_histograms.png", dpi=220)
    figure.savefig(OUTDIR / "energy_coupling_histograms.pdf")
    plt.close(figure)

    wavelength_figure, wavelength_axes = plt.subplots(1, 2, figsize=(10, 4))
    raw_nm = np.concatenate([HC_EV_NM / ea, HC_EV_NM / eb])
    corrected_nm = np.concatenate([HC_EV_NM / eac, HC_EV_NM / ebc])
    wavelength_axes[0].hist(raw_nm, bins="auto", color="#4c78a8")
    wavelength_axes[0].set(xlabel="Raw photon wavelength (nm)", ylabel="Site-frame samples")
    wavelength_axes[1].hist(corrected_nm, bins="auto", color="#f58518")
    wavelength_axes[1].set(
        xlabel="Same-geometry STEOM-shifted wavelength (nm)",
        ylabel="Site-frame samples",
    )
    wavelength_figure.tight_layout()
    wavelength_figure.savefig(OUTDIR / "wavelength_histograms.png", dpi=220)
    wavelength_figure.savefig(OUTDIR / "wavelength_histograms.pdf")
    plt.close(wavelength_figure)

    failure_audit = {
        "requested_frames": 50,
        "complete_frames": len(complete),
        "failed_frames": len(failed),
        "no_interpolation": True,
        "failures": [
            {
                "frame_index": int(row["frame_index"]),
                "time_ps": float(row["time_ps"]),
                "failure_class": row["failure_class"],
                "retry_count": int(row["retry_count"]),
                "notes": row["notes"],
            }
            for row in failed
        ],
    }
    (OUTDIR / "failure_audit.json").write_text(
        json.dumps(failure_audit, indent=2) + "\n"
    )
    report = [
        "# Corrected 50-frame TDDFT sampling over 100 ns",
        "",
        f"- Complete frames: {len(complete)}/50",
        f"- Classified failures: {len(failed)}/50",
        "- Physical times: 2, 4, ..., 100 ns",
        "- Link caps: instantaneous severed-bond directions; links excluded from TDC",
        f"- Rapid TDDFT method: {calibration['method']}/{calibration['basis']}",
        f"- Same-geometry STEOM reference correction: {delta:+.8f} eV applied equally to both sites",
        (
            f"- Raw site wavelengths: A "
            f"{result['distributions']['site_A_raw_nm']['mean']:.3f} ± "
            f"{result['distributions']['site_A_raw_nm']['sample_sd']:.3f} nm; B "
            f"{result['distributions']['site_B_raw_nm']['mean']:.3f} ± "
            f"{result['distributions']['site_B_raw_nm']['sample_sd']:.3f} nm"
        ),
        (
            f"- STEOM-shifted site wavelengths: A "
            f"{result['distributions']['site_A_STEOM_shifted_nm']['mean']:.3f} ± "
            f"{result['distributions']['site_A_STEOM_shifted_nm']['sample_sd']:.3f} nm; B "
            f"{result['distributions']['site_B_STEOM_shifted_nm']['mean']:.3f} ± "
            f"{result['distributions']['site_B_STEOM_shifted_nm']['sample_sd']:.3f} nm"
        ),
        (
            f"- Screened signed TDC: "
            f"{result['distributions']['J_screened_cm-1']['mean']:.3f} ± "
            f"{result['distributions']['J_screened_cm-1']['sample_sd']:.3f} cm^-1"
        ),
        (
            "- Autocorrelation-aware uncertainty: circular block-bootstrap "
            "95% confidence intervals with per-observable positive-sequence "
            "autocorrelation block lengths (stored in summary.json)"
        ),
        f"- Constant-shift A-B gap invariance: {gap_invariance:.3e} eV",
        "",
        "The 2 ns sampling interval cannot resolve femtosecond decoherence; these",
        "frames quantify slow conformational and inhomogeneous disorder only.",
        "",
    ]
    (OUTDIR / "NUMERICAL_REPORT.md").write_text("\n".join(report))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

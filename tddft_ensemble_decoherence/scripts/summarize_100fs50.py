#!/usr/bin/env python3
"""Analyze a contiguous 100 fs paired-site TDDFT/TDC decoherence window."""

from __future__ import annotations

import argparse
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
MANIFEST = ROOT / "manifests/short100fs50_frames.csv"
OUTDIR = ROOT / "rapid_100fs50"
HBAR_EV_FS = 0.6582119569
EV_TO_CM = 8065.544005
HC_EV_NM = 1239.8419843320026


def distribution(values: np.ndarray) -> dict:
    return {
        "n": int(len(values)),
        "mean": float(np.mean(values)),
        "sample_sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "median": float(np.median(values)),
        "iqr": float(np.quantile(values, 0.75) - np.quantile(values, 0.25)),
        "q025": float(np.quantile(values, 0.025)),
        "q05": float(np.quantile(values, 0.05)),
        "q95": float(np.quantile(values, 0.95)),
        "q975": float(np.quantile(values, 0.975)),
    }


def acf_unbiased(values: np.ndarray) -> np.ndarray:
    centered = values - np.mean(values)
    n = len(centered)
    return np.asarray(
        [
            np.dot(centered[: n - lag], centered[lag:]) / (n - lag)
            for lag in range(n)
        ],
        dtype=float,
    )


def cumulant(values: np.ndarray, dt_fs: float) -> dict:
    correlation = acf_unbiased(values)
    time_fs = np.arange(len(values), dtype=float) * dt_fs
    g = np.zeros(len(values), dtype=float)
    for index in range(1, len(values)):
        tau = time_fs[: index + 1]
        integrand = (time_fs[index] - tau) * correlation[: index + 1]
        g[index] = float(np.trapezoid(integrand, tau) / HBAR_EV_FS**2)
    coherence = np.exp(-g)
    return {
        "time_fs": time_fs,
        "correlation_eV2": correlation,
        "g": g,
        "coherence": coherence,
        "T2star_fs": first_crossing(time_fs, coherence),
    }


def first_crossing(time_fs: np.ndarray, coherence: np.ndarray) -> float | None:
    target = np.exp(-1.0)
    indices = np.flatnonzero(coherence <= target)
    if not len(indices):
        return None
    index = int(indices[0])
    if index == 0:
        return float(time_fs[0])
    x0, x1 = time_fs[index - 1], time_fs[index]
    y0, y1 = coherence[index - 1], coherence[index]
    if y1 == y0:
        return float(x1)
    return float(x0 + (target - y0) * (x1 - x0) / (y1 - y0))


def positive_correlation_time(correlation: np.ndarray, dt_fs: float) -> float:
    if correlation[0] <= 0.0:
        return 0.0
    stop = len(correlation)
    nonpositive = np.flatnonzero(correlation[1:] <= 0.0)
    if len(nonpositive):
        stop = int(nonpositive[0]) + 1
    normalized = correlation[: stop + 1] / correlation[0]
    return float(np.trapezoid(normalized, dx=dt_fs))


def moving_block_bootstrap_t2(
    values: np.ndarray,
    dt_fs: float,
    block_samples: int,
    replicates: int = 1000,
) -> dict:
    """Autocorrelation-aware circular moving-block uncertainty for T2*."""
    rng = np.random.default_rng(206)
    n = len(values)
    estimates = []
    censored = 0
    blocks_needed = int(np.ceil(n / block_samples))
    offsets = np.arange(block_samples)
    for _ in range(replicates):
        starts = rng.integers(0, n, size=blocks_needed)
        indices = np.concatenate([(start + offsets) % n for start in starts])[:n]
        estimate = cumulant(values[indices], dt_fs)["T2star_fs"]
        if estimate is None:
            censored += 1
        else:
            estimates.append(float(estimate))
    finite = np.asarray(estimates, dtype=float)
    return {
        "method": "circular moving-block bootstrap",
        "random_seed": 206,
        "replicates": replicates,
        "block_samples": block_samples,
        "block_duration_fs": block_samples * dt_fs,
        "finite_estimate_count": int(len(finite)),
        "right_censored_count": censored,
        "finite_percentile_95_CI_fs": (
            [float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975))]
            if len(finite)
            else [None, None]
        ),
        "finite_median_fs": float(np.median(finite)) if len(finite) else None,
    }


def direct_coherence(
    ea_eV: np.ndarray,
    eb_eV: np.ndarray,
    j_cm: np.ndarray,
    dt_fs: float,
) -> np.ndarray:
    n = len(ea_eV)
    accumulated = np.zeros(n, dtype=np.complex128)
    counts = np.zeros(n, dtype=int)
    initial = np.asarray([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
    for origin in range(n):
        psi = initial.copy()
        accumulated[0] += psi[0] * np.conj(psi[1])
        counts[0] += 1
        for step in range(origin, n - 1):
            h = np.asarray(
                [
                    [ea_eV[step], j_cm[step] / EV_TO_CM],
                    [j_cm[step] / EV_TO_CM, eb_eV[step]],
                ],
                dtype=float,
            )
            h -= np.eye(2) * np.trace(h) / 2.0
            eigenvalues, eigenvectors = np.linalg.eigh(h)
            phase = np.exp(-1j * eigenvalues * dt_fs / HBAR_EV_FS)
            psi = eigenvectors @ (phase * (eigenvectors.T @ psi))
            lag = step - origin + 1
            accumulated[lag] += psi[0] * np.conj(psi[1])
            counts[lag] += 1
    return np.abs(accumulated / counts) / 0.5


def normalized_cross_correlation(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    a = left - np.mean(left)
    b = right - np.mean(right)
    denominator = np.std(a) * np.std(b)
    if denominator == 0.0:
        return np.zeros(len(a))
    n = len(a)
    return np.asarray(
        [
            np.dot(a[: n - lag], b[lag:]) / (n - lag) / denominator
            for lag in range(n)
        ]
    )


def main() -> None:
    global MANIFEST, OUTDIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--window",
        choices=(
            "short100fs50",
            "short100fs50_fixedmm",
            "short100fs100_wb97xd3",
            # 2026-08-31 T2 campaign: 1 fs cadence on the corrected 511-aa
            # construct with full-system embedding. Everything above predates
            # both fixes and runs at 2-4 fs.
            "v6t2a",
            "v6t2b",
            "v6t2c",
        ),
        default="short100fs50",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    MANIFEST = ROOT / "manifests" / f"{args.window}_frames.csv"
    OUTDIR = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else ROOT
        / (
            "rapid_100fs50_fixedmm"
            if args.window == "short100fs50_fixedmm"
            else "wb97xd3_100fs100"
            if args.window == "short100fs100_wb97xd3"
            else f"t2_{args.window}"
            if args.window.startswith("v6t2")
            else "rapid_100fs50"
        )
    )
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected_frames = (
        100
        if args.window == "short100fs100_wb97xd3" or args.window.startswith("v6t2")
        else 50
    )
    if len(rows) != expected_frames:
        raise RuntimeError(
            f"Expected {expected_frames} manifest rows, found {len(rows)}"
        )
    complete = [row for row in rows if row["convergence_status"] == "complete"]
    failed = [row for row in rows if row["convergence_status"] == "failed"]
    unresolved = [
        row for row in rows if row["convergence_status"] not in {"complete", "failed"}
    ]
    if unresolved:
        raise RuntimeError(
            "Campaign still has unclassified frames: "
            + ",".join(row["frame_index"] for row in unresolved)
        )
    if failed:
        raise RuntimeError(
            "Decoherence analysis requires the contiguous measured series; "
            "failed frames are not interpolated: "
            + ",".join(row["frame_index"] for row in failed)
        )

    time_fs = np.asarray([float(row["time_fs"]) for row in complete])
    dt_fs = float(np.median(np.diff(time_fs)))
    if not np.allclose(np.diff(time_fs), dt_fs) or dt_fs <= 0.0:
        raise RuntimeError("The complete series does not have a contiguous cadence")
    ea = np.asarray([float(row["site_A_energy_eV"]) for row in complete])
    eb = np.asarray([float(row["site_B_energy_eV"]) for row in complete])
    eac = np.asarray([float(row["site_A_corrected_eV"]) for row in complete])
    ebc = np.asarray([float(row["site_B_corrected_eV"]) for row in complete])
    jvac = np.asarray([float(row["J_vacuum_cm-1"]) for row in complete])
    jscr = np.asarray([float(row["J_screened_cm-1"]) for row in complete])
    jpda = np.asarray([float(row["J_pda_screened_cm-1"]) for row in complete])
    separation = np.asarray([float(row["separation_A"]) for row in complete])
    gap = ea - eb
    corrected_gap = eac - ebc
    gap_invariance = float(np.max(np.abs(gap - corrected_gap)))
    if gap_invariance > 2.0e-9:
        raise RuntimeError(f"Constant-shift gap invariance failed: {gap_invariance}")

    raw = cumulant(gap, dt_fs)
    corrected = cumulant(corrected_gap, dt_fs)
    if not np.allclose(raw["coherence"], corrected["coherence"], atol=1.0e-12):
        raise RuntimeError("Constant STEOM shift changed cumulant coherence")
    stride_results = {}
    for stride in (1, 2, 4):
        result = cumulant(gap[::stride], dt_fs * stride)
        stride_results[str(stride)] = {
            "sampling_interval_fs": dt_fs * stride,
            "samples": int(len(gap[::stride])),
            "T2star_fs": result["T2star_fs"],
        }
    segment_results = {}
    midpoint = len(gap) // 2
    for name, indices in {
        "first_half": slice(0, midpoint),
        "last_half": slice(midpoint, len(gap)),
    }.items():
        segment_samples = len(gap[indices])
        result = cumulant(gap[indices], dt_fs)
        segment_results[name] = {
            "samples": int(segment_samples),
            "T2star_fs": result["T2star_fs"],
            "censored_at_fs": (
                (segment_samples - 1) * dt_fs
                if result["T2star_fs"] is None
                else None
            ),
        }

    direct_full = direct_coherence(ea, eb, jscr, dt_fs)
    direct_gap_only = direct_coherence(ea, eb, np.zeros_like(jscr), dt_fs)
    lag_fs = np.arange(len(gap), dtype=float) * dt_fs
    direct_full_t2 = first_crossing(lag_fs, direct_full)
    direct_gap_t2 = first_crossing(lag_fs, direct_gap_only)
    coupling_effect = float(np.max(np.abs(direct_full - direct_gap_only)))

    correlation_time = positive_correlation_time(raw["correlation_eV2"], dt_fs)
    bootstrap_block_samples = max(
        3, int(np.ceil(2.0 * correlation_time / dt_fs))
    )
    t2_bootstrap = moving_block_bootstrap_t2(
        gap, dt_fs, bootstrap_block_samples
    )
    correlations = {
        "pearson_zero_lag": {
            "EA_EB": float(np.corrcoef(ea, eb)[0, 1]),
            "EA_J": float(np.corrcoef(ea, jscr)[0, 1]),
            "EB_J": float(np.corrcoef(eb, jscr)[0, 1]),
            "gap_J": float(np.corrcoef(gap, jscr)[0, 1]),
        },
        "lag_fs": lag_fs.tolist(),
        "EA_autocorrelation_normalized": (
            acf_unbiased(ea) / acf_unbiased(ea)[0]
        ).tolist(),
        "EB_autocorrelation_normalized": (
            acf_unbiased(eb) / acf_unbiased(eb)[0]
        ).tolist(),
        "gap_autocorrelation_normalized": (
            raw["correlation_eV2"] / raw["correlation_eV2"][0]
        ).tolist(),
        "J_autocorrelation_normalized": (
            acf_unbiased(jscr) / acf_unbiased(jscr)[0]
        ).tolist(),
        "gap_to_J_cross_correlation_normalized": normalized_cross_correlation(
            gap, jscr
        ).tolist(),
    }
    calibration = json.loads((OUTDIR / "calibration.json").read_text())
    result = {
        "status": "complete",
        "requested_frames": expected_frames,
        "complete_frames": expected_frames,
        "failed_frames": 0,
        "sampling_interval_fs": dt_fs,
        "coverage_fs": [float(time_fs[0]), float(time_fs[-1])],
        "method": calibration["method"],
        "basis": calibration["basis"],
        "production_roots": calibration.get("production_roots", calibration["roots"]),
        "steom_constant_shift_eV": calibration["steom_constant_shift_eV"],
        "constant_shift_gap_invariance_max_abs_eV": gap_invariance,
        "cumulant": {
            "T2star_raw_TDDFT_fs": raw["T2star_fs"],
            "T2star_constant_shifted_STEOM_fs": corrected["T2star_fs"],
            "bath_positive_correlation_time_fs": correlation_time,
            "time_fs": raw["time_fs"].tolist(),
            "C_gap_eV2": raw["correlation_eV2"].tolist(),
            "g": raw["g"].tolist(),
            "coherence": raw["coherence"].tolist(),
            "stride_sensitivity": stride_results,
            "segment_sensitivity": segment_results,
            "moving_block_bootstrap": t2_bootstrap,
        },
        "direct_two_site_propagation": {
            "T2_e_minus_1_full_signed_J_fs": direct_full_t2,
            "T2_e_minus_1_gap_only_fs": direct_gap_t2,
            "maximum_absolute_coherence_difference_due_to_J": coupling_effect,
            "time_fs": lag_fs.tolist(),
            "coherence_full_signed_J": direct_full.tolist(),
            "coherence_gap_only": direct_gap_only.tolist(),
        },
        "distributions": {
            "site_A_raw_eV": distribution(ea),
            "site_B_raw_eV": distribution(eb),
            "site_A_corrected_eV": distribution(eac),
            "site_B_corrected_eV": distribution(ebc),
            "site_A_corrected_nm": distribution(HC_EV_NM / eac),
            "site_B_corrected_nm": distribution(HC_EV_NM / ebc),
            "mean_site_raw_eV": distribution((ea + eb) / 2.0),
            "differential_gap_eV": distribution(gap),
            "J_vacuum_cm-1": distribution(jvac),
            "J_screened_cm-1": distribution(jscr),
            "absolute_J_screened_cm-1": distribution(np.abs(jscr)),
            "Davydov_splitting_screened_cm-1": distribution(2.0 * np.abs(jscr)),
            "J_PDA_screened_cm-1": distribution(jpda),
            "separation_A": distribution(separation),
        },
        "selected_root_counts": {
            "site_A": dict(Counter(row["site_A_selected_root"] for row in complete)),
            "site_B": dict(Counter(row["site_B_selected_root"] for row in complete)),
        },
        "correlations": correlations,
        "fixed_charge_comparison_fs": {
            "direct_minimum_image": 21.8,
            "independent_PME": 24.0,
            "one_ps_block_median": 24.5,
            "one_ps_block_SD": 3.7,
            "manuscript_phenomenological_sensitivity": 60.0,
        },
        "interpretation": {
            "bath_correlation_time": (
                "Positive-sequence integral of the normalized classical gap ACF."
            ),
            "pure_dephasing": (
                "Classical Gaussian second-cumulant estimate from geometry-dependent "
                "TDDFT differential-gap fluctuations."
            ),
            "inhomogeneous_broadening": (
                "Not estimable from this 100 fs window; requires a long ensemble."
            ),
            "population_relaxation": "Not determined by this calculation.",
        },
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MANIFEST, OUTDIR / "frame_manifest.csv")
    (OUTDIR / "summary.json").write_text(json.dumps(result, indent=2) + "\n")

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(lag_fs, raw["correlation_eV2"], color="#4c78a8")
    axes[0].axhline(0.0, color="black", linewidth=0.7)
    axes[0].set(xlabel="Lag (fs)", ylabel=r"$C_X(t)$ (eV$^2$)")
    axes[1].plot(lag_fs, raw["coherence"], label="Gaussian cumulant")
    axes[1].plot(lag_fs, direct_gap_only, label="direct, gap only")
    axes[1].plot(lag_fs, direct_full, label="direct, signed $J(t)$")
    axes[1].axhline(np.exp(-1.0), color="black", linestyle="--", linewidth=0.8)
    axes[1].set(xlabel="Lag (fs)", ylabel="Normalized coherence", ylim=(0.0, 1.08))
    axes[1].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(OUTDIR / "correlation_and_coherence.png", dpi=220)
    figure.savefig(OUTDIR / "correlation_and_coherence.pdf")
    plt.close(figure)

    figure, axes = plt.subplots(2, 3, figsize=(13.5, 7.5))
    shared_energy_bins = np.histogram_bin_edges(np.concatenate([ea, eb]), bins="auto")
    axes[0, 0].hist(ea, bins=shared_energy_bins, alpha=0.62, label="A")
    axes[0, 0].hist(eb, bins=shared_energy_bins, alpha=0.62, label="B")
    axes[0, 0].legend(frameon=False)
    axes[0, 0].set(xlabel="Raw excitation energy (eV)", ylabel="Frames")
    axes[0, 1].hist((ea + eb) / 2.0, bins="auto")
    axes[0, 1].set(xlabel="Mean site energy (eV)", ylabel="Frames")
    axes[0, 2].hist(gap, bins="auto")
    axes[0, 2].set(xlabel=r"$E_A-E_B$ (eV)", ylabel="Frames")
    axes[1, 0].hist(jscr, bins="auto")
    axes[1, 0].axvline(0.0, color="black", linewidth=0.7)
    axes[1, 0].set(xlabel=r"Signed screened $J$ (cm$^{-1}$)", ylabel="Frames")
    axes[1, 1].hist(np.abs(jscr), bins="auto")
    axes[1, 1].set(xlabel=r"$|J|$ (cm$^{-1}$)", ylabel="Frames")
    axes[1, 2].hist(2.0 * np.abs(jscr), bins="auto")
    axes[1, 2].set(xlabel=r"$2|J|$ (cm$^{-1}$)", ylabel="Frames")
    figure.tight_layout()
    figure.savefig(OUTDIR / "energy_coupling_histograms.png", dpi=220)
    figure.savefig(OUTDIR / "energy_coupling_histograms.pdf")
    plt.close(figure)

    t2 = raw["T2star_fs"]
    report = [
        "# Rapid 100 fs TDDFT/TDC decoherence report",
        "",
        (
            f"- Frames: {expected_frames}/{expected_frames} at "
            f"{dt_fs:g} fs cadence; no interpolation"
        ),
        f"- Method: {calibration['method']}/{calibration['basis']}",
        (
            f"- Same-geometry STEOM shift: "
            f"{float(calibration['steom_constant_shift_eV']):+.8f} eV"
        ),
        (
            f"- Cumulant T2*: {t2:.3f} fs"
            if t2 is not None
            else f"- Cumulant T2*: >{lag_fs[-1]:g} fs (no crossing)"
        ),
        (
            f"- Positive-sequence bath correlation time: "
            f"{correlation_time:.3f} fs"
        ),
        (
            f"- Direct propagation e^-1 crossing: full signed J "
            f"{direct_full_t2 if direct_full_t2 is not None else f'>{lag_fs[-1]:g}'} fs; "
            f"gap-only {direct_gap_t2 if direct_gap_t2 is not None else f'>{lag_fs[-1]:g}'} fs"
        ),
        f"- Constant-shift gap invariance: {gap_invariance:.3e} eV",
        (
            "- Previous classical estimates: 21.8 fs minimum-image, 24.0 fs PME, "
            "24.5 ± 3.7 fs across 1 ps blocks; manuscript sensitivity value 60 fs"
        ),
        "",
        "This 100 fs window probes ultrafast correlations and dephasing but cannot",
        "characterize long-time inhomogeneous broadening. Population relaxation is",
        "not determined by this calculation.",
        "",
    ]
    (OUTDIR / "NUMERICAL_REPORT.md").write_text("\n".join(report))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

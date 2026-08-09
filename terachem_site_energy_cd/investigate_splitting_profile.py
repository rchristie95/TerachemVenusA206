#!/usr/bin/env python3
"""Numerical identifiability study for the Venus CD splitting profile.

This program launches no electronic-structure calculations.  It combines the
audited corrected-coupling/transition-dipole geometry archive with a general
nondegenerate two-site Hamiltonian and the published constrained dVenus
Lorentzian parameters.  Its purpose is to determine which combinations of
coupling and site-energy detuning can reproduce the *profile*, and which
quantities remain unidentifiable without frame-resolved site energies or an
absolute CD calibration.

The archived J/geometry ensemble is used only as an ensemble sensitivity
input.  No TDDFT trajectory is joined to it here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = ROOT / "coupling_nvt_production_cr2_1000_20260721" / "coupling_geometry.npz"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "results" / "splitting_profile_investigation"

EV_TO_CM = 8065.544005
C_CM_PER_S = 2.99792458e10

# Nguyen et al., Biophysical Journal 124, 4293--4309 (2025), Table S3.
EXP_CENTER_CM = 19322.0
EXP_HALF_SPLITTING_CM = 130.79
EXP_FULL_SPLITTING_CM = 2.0 * EXP_HALF_SPLITTING_CM
EXP_LOW_AMPLITUDE = 3.59
EXP_LOW_HWHM_CM = 293.51
EXP_HIGH_AMPLITUDE = -2.27
EXP_HIGH_HWHM_CM = 439.25
EXP_THIRD_CENTER_CM = 20774.0
EXP_THIRD_AMPLITUDE = -0.34
EXP_THIRD_HWHM_CM = 565.79

STEOM_SITE_CM = 19088.2
STEOM_SITE_EV = 2.366633
TDDFT_REFERENCE_EV = 3.43441587
COMMON_STEOM_SHIFT_EV = STEOM_SITE_EV - TDDFT_REFERENCE_EV


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def lorentzian_area(offset_cm: np.ndarray, hwhm_cm: float) -> np.ndarray:
    return (hwhm_cm / np.pi) / (offset_cm * offset_cm + hwhm_cm * hwhm_cm)


def lorentzian_peak(
    grid_cm: np.ndarray, amplitude: float, center_cm: float, hwhm_cm: float
) -> np.ndarray:
    return amplitude / (1.0 + ((grid_cm - center_cm) / hwhm_cm) ** 2)


def load_ensemble(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as archive:
        required = {"frame", "J_cm", "mu_A", "mu_B", "r_A", "r_B"}
        missing = required.difference(archive.files)
        if missing:
            raise ValueError(f"{path} lacks arrays: {', '.join(sorted(missing))}")
        data = {name: np.asarray(archive[name], dtype=float) for name in required}
    n = len(data["J_cm"])
    expected = {
        "frame": (n,),
        "J_cm": (n,),
        "mu_A": (n, 3),
        "mu_B": (n, 3),
        "r_A": (n, 3),
        "r_B": (n, 3),
    }
    for name, shape in expected.items():
        if data[name].shape != shape:
            raise ValueError(f"{name} has shape {data[name].shape}; expected {shape}")
        if not np.all(np.isfinite(data[name])):
            raise ValueError(f"Non-finite values in {path}:{name}")
    data["triple_product"] = np.einsum(
        "ij,ij->i",
        data["r_B"] - data["r_A"],
        np.cross(data["mu_A"], data["mu_B"]),
    )
    return data


def experimental_profile(grid_cm: np.ndarray, include_third: bool = True) -> np.ndarray:
    profile = lorentzian_peak(
        grid_cm,
        EXP_LOW_AMPLITUDE,
        EXP_CENTER_CM - EXP_HALF_SPLITTING_CM,
        EXP_LOW_HWHM_CM,
    )
    profile += lorentzian_peak(
        grid_cm,
        EXP_HIGH_AMPLITUDE,
        EXP_CENTER_CM + EXP_HALF_SPLITTING_CM,
        EXP_HIGH_HWHM_CM,
    )
    if include_third:
        profile += lorentzian_peak(
            grid_cm,
            EXP_THIRD_AMPLITUDE,
            EXP_THIRD_CENTER_CM,
            EXP_THIRD_HWHM_CM,
        )
    return profile


def model_components(
    grid_cm: np.ndarray,
    ensemble: dict[str, np.ndarray],
    detuning_cm: float | np.ndarray,
    *,
    j_scale: float = 1.0,
    center_cm: float = EXP_CENTER_CM,
    low_hwhm_cm: float = EXP_LOW_HWHM_CM,
    high_hwhm_cm: float = EXP_HIGH_HWHM_CM,
    chunk_size: int = 128,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Return low/high interaction-CD components and per-frame diagnostics."""
    j_cm = j_scale * ensemble["J_cm"]
    detuning = np.broadcast_to(np.asarray(detuning_cm, dtype=float), j_cm.shape)
    omega = np.hypot(detuning, 2.0 * j_cm)
    low_energy = center_cm - 0.5 * omega
    high_energy = center_cm + 0.5 * omega

    # For the two-site exciton-chirality term, c_A*c_B = +/- J/Omega.
    # The common absolute-unit prefactor is intentionally omitted.
    geometric = ensemble["triple_product"]
    low_rotation = +np.pi * low_energy * (j_cm / omega) * geometric
    high_rotation = -np.pi * high_energy * (j_cm / omega) * geometric

    low = np.zeros_like(grid_cm)
    high = np.zeros_like(grid_cm)
    for start in range(0, len(j_cm), chunk_size):
        stop = min(start + chunk_size, len(j_cm))
        sl = slice(start, stop)
        low += np.sum(
            low_rotation[sl, None]
            * lorentzian_area(grid_cm[None, :] - low_energy[sl, None], low_hwhm_cm),
            axis=0,
        )
        high += np.sum(
            high_rotation[sl, None]
            * lorentzian_area(grid_cm[None, :] - high_energy[sl, None], high_hwhm_cm),
            axis=0,
        )
    n = len(j_cm)
    return low / n, high / n, {
        "omega_cm": omega,
        "mixing": 2.0 * np.abs(j_cm) / omega,
        "low_energy_cm": low_energy,
        "high_energy_cm": high_energy,
        "low_rotation_relative": low_rotation,
        "high_rotation_relative": high_rotation,
    }


def normalized_extrema(grid_cm: np.ndarray, profile: np.ndarray) -> dict[str, float]:
    mask = (grid_cm >= 18000.0) & (grid_cm <= 20500.0)
    indices = np.flatnonzero(mask)
    i_max = int(indices[np.argmax(profile[mask])])
    i_min = int(indices[np.argmin(profile[mask])])
    return {
        "maximum_cm-1": float(grid_cm[i_max]),
        "minimum_cm-1": float(grid_cm[i_min]),
        "extrema_separation_cm-1": float(abs(grid_cm[i_max] - grid_cm[i_min])),
        "maximum_nm": float(1.0e7 / grid_cm[i_max]),
        "minimum_nm": float(1.0e7 / grid_cm[i_min]),
    }


def fit_reconstructed_profile(
    grid_cm: np.ndarray,
    target: np.ndarray,
    ensemble: dict[str, np.ndarray],
    j_scale: float,
) -> dict[str, float]:
    """Fit detuning while profiling linear scale/background nuisance terms.

    This fits the Table-S3 reconstruction, not raw experimental observations;
    the residual is therefore a shape diagnostic and has no statistical-error
    interpretation.
    """
    third = lorentzian_peak(grid_cm, 1.0, EXP_THIRD_CENTER_CM, EXP_THIRD_HWHM_CM)
    constant = np.ones_like(grid_cm)
    slope = (grid_cm - EXP_CENTER_CM) / 1000.0
    target_range = float(np.ptp(target))

    def evaluate(detuning_cm: float, details: bool = False):
        low, high, diagnostics = model_components(
            grid_cm, ensemble, detuning_cm, j_scale=j_scale
        )
        model = low + high
        design = np.column_stack([model, third, constant, slope])
        coefficients, _, _, _ = np.linalg.lstsq(design, target, rcond=None)
        residual = target - design @ coefficients
        normalized_rms = float(np.sqrt(np.mean(residual * residual)) / target_range)
        if not details:
            return normalized_rms
        omega = diagnostics["omega_cm"]
        return normalized_rms, coefficients, model, residual, diagnostics, omega

    result = minimize_scalar(
        evaluate,
        method="bounded",
        bounds=(0.0, 700.0),
        options={"xatol": 0.01},
    )
    rms, coefficients, model, residual, diagnostics, omega = evaluate(result.x, True)
    extrema = normalized_extrema(grid_cm, coefficients[0] * model)
    return {
        "J_scale": float(j_scale),
        "J_mean_cm-1": float(np.mean(j_scale * ensemble["J_cm"])),
        "best_detuning_cm-1": float(result.x),
        "best_gap_mean_cm-1": float(np.mean(omega)),
        "best_gap_sample_sd_cm-1": float(np.std(omega, ddof=1)),
        "best_mean_mixing": float(np.mean(diagnostics["mixing"])),
        "best_mean_participation_ratio": float(
            np.mean(1.0 / (1.0 - 0.5 * diagnostics["mixing"] ** 2))
        ),
        "normalized_rms_profile_residual": float(rms),
        "profile_scale": float(coefficients[0]),
        "third_band_amplitude": float(coefficients[1]),
        "baseline_constant": float(coefficients[2]),
        "baseline_slope_per_1000cm-1": float(coefficients[3]),
        **extrema,
    }


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def hwhm_from_t2(t2_fs: float) -> float:
    return 1.0 / (2.0 * np.pi * C_CM_PER_S * t2_fs * 1e-15)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ensemble = load_ensemble(args.archive)
    j_cm = ensemble["J_cm"]

    required_delta = np.sqrt(np.maximum(EXP_FULL_SPLITTING_CM**2 - 4.0 * j_cm**2, 0.0))
    unconstrained_gap = 1.0e7 / 511.8 - 1.0e7 / 521.7
    required_delta_unconstrained = np.sqrt(
        np.maximum(unconstrained_gap**2 - 4.0 * j_cm**2, 0.0)
    )

    profile_grid = np.linspace(18000.0, 21500.0, 1401)
    target = experimental_profile(profile_grid, include_third=True)
    target_couplet = experimental_profile(profile_grid, include_third=False)

    scenario_definitions = [
        ("resonant_corrected_J", 0.0),
        ("100_cm-1_detuning", 100.0),
        ("constrained_fit_gap", float(np.mean(required_delta))),
        ("unconstrained_fit_gap", float(np.mean(required_delta_unconstrained))),
        ("corrected_charge_smoke_diagnostic", 1915.98),
    ]
    scenario_rows: list[dict[str, float | str]] = []
    profiles: dict[str, np.ndarray] = {}
    for name, detuning in scenario_definitions:
        low, high, diagnostics = model_components(profile_grid, ensemble, detuning)
        # Nguyen's experimental difference is TDX - TD.  The uncoupled
        # interaction-induced exciton-chirality contribution is zero in this
        # reduced model, so its corresponding difference is -(low + high).
        profile = -(low + high)
        profiles[name] = profile
        mixing = diagnostics["mixing"]
        omega = diagnostics["omega_cm"]
        extrema = normalized_extrema(profile_grid, profile)
        scenario_rows.append(
            {
                "scenario": name,
                "detuning_cm-1": detuning,
                "detuning_eV": detuning / EV_TO_CM,
                "gap_mean_cm-1": float(np.mean(omega)),
                "gap_sample_sd_cm-1": float(np.std(omega, ddof=1)),
                "mixing_mean": float(np.mean(mixing)),
                "minor_site_weight_mean": float(
                    np.mean(0.5 * (1.0 - np.sqrt(1.0 - mixing * mixing)))
                ),
                "participation_ratio_mean": float(
                    np.mean(1.0 / (1.0 - 0.5 * mixing * mixing))
                ),
                "relative_CD_area_vs_resonance": float(
                    np.mean(np.abs(diagnostics["low_rotation_relative"]))
                    / np.mean(
                        np.abs(
                            model_components(profile_grid[:2], ensemble, 0.0)[2][
                                "low_rotation_relative"
                            ]
                        )
                    )
                ),
                **extrema,
            }
        )

    j_scales = [0.25, 0.5, 1.0, 2.0, 3.0, 3.9, 4.1]
    ridge_rows = [
        fit_reconstructed_profile(profile_grid, target, ensemble, scale)
        for scale in j_scales
    ]

    # Exact energy-gap ridge implied by component centers alone.
    center_ridge_rows: list[dict[str, float]] = []
    for assumed_j in np.linspace(0.0, EXP_HALF_SPLITTING_CM, 132):
        detuning = np.sqrt(max(EXP_FULL_SPLITTING_CM**2 - 4.0 * assumed_j**2, 0.0))
        mixing = 0.0 if EXP_FULL_SPLITTING_CM == 0 else 2.0 * assumed_j / EXP_FULL_SPLITTING_CM
        center_ridge_rows.append(
            {
                "assumed_J_cm-1": float(assumed_j),
                "required_detuning_cm-1": float(detuning),
                "mixing": float(mixing),
                "minor_site_weight": float(0.5 * (1.0 - np.sqrt(1.0 - mixing * mixing))),
            }
        )

    # Differential-disorder sensitivity.  Gauss-Hermite nodes give a
    # deterministic normal average without Monte Carlo noise.
    nodes, weights = np.polynomial.hermite.hermgauss(15)
    weights = weights / np.sqrt(np.pi)
    disorder_rows: list[dict[str, float | str]] = []
    base_delta = float(np.mean(required_delta))
    for sigma_delta in [0.0, 25.0, 50.0, 100.0, 200.0, 400.0]:
        all_omega = []
        all_mixing = []
        branch_means = []
        for node, weight in zip(nodes, weights):
            delta = base_delta + np.sqrt(2.0) * sigma_delta * node
            omega = np.hypot(delta, 2.0 * j_cm)
            all_omega.append((weight, omega))
            all_mixing.append((weight, 2.0 * np.abs(j_cm) / omega))
            branch_means.append((weight, 0.5 * omega))
        mean_omega = sum(weight * np.mean(value) for weight, value in all_omega)
        second_omega = sum(weight * np.mean(value * value) for weight, value in all_omega)
        mean_mixing = sum(weight * np.mean(value) for weight, value in all_mixing)
        mean_branch = sum(weight * np.mean(value) for weight, value in branch_means)
        second_branch = sum(weight * np.mean(value * value) for weight, value in branch_means)
        disorder_rows.append(
            {
                "distribution": "normal_about_required_mean",
                "detuning_mean_cm-1": base_delta,
                "detuning_normal_sd_cm-1": sigma_delta,
                "gap_mean_cm-1": float(mean_omega),
                "gap_rms_sd_cm-1": float(np.sqrt(max(second_omega - mean_omega**2, 0.0))),
                "mixing_mean": float(mean_mixing),
                "single_branch_rms_sd_cm-1": float(
                    np.sqrt(max(second_branch - mean_branch**2, 0.0))
                ),
            }
        )

    # A nominal homodimer need not have a persistent signed bias.  Solve the
    # alternative in which Delta is normally distributed about zero and its
    # variance alone produces the observed mean gap.
    zero_nodes, zero_weights = np.polynomial.hermite.hermgauss(80)
    zero_weights = zero_weights / np.sqrt(np.pi)

    def zero_mean_disorder_stats(sigma_delta: float) -> tuple[float, float, float]:
        delta = np.sqrt(2.0) * sigma_delta * zero_nodes[:, None]
        omega = np.hypot(delta, 2.0 * j_cm[None, :])
        mean_omega = float(np.sum(zero_weights[:, None] * omega) / len(j_cm))
        second_omega = float(
            np.sum(zero_weights[:, None] * omega * omega) / len(j_cm)
        )
        mean_mixing = float(
            np.sum(
                zero_weights[:, None] * (2.0 * np.abs(j_cm)[None, :] / omega)
            )
            / len(j_cm)
        )
        return mean_omega, np.sqrt(max(second_omega - mean_omega**2, 0.0)), mean_mixing

    zero_mean_sigma = brentq(
        lambda sigma: zero_mean_disorder_stats(sigma)[0] - EXP_FULL_SPLITTING_CM,
        1.0,
        1000.0,
    )
    zero_mean_gap, zero_mean_gap_sd, zero_mean_mixing = zero_mean_disorder_stats(
        zero_mean_sigma
    )
    disorder_rows.append(
        {
            "distribution": "zero_mean_normal_target_gap",
            "detuning_mean_cm-1": 0.0,
            "detuning_normal_sd_cm-1": float(zero_mean_sigma),
            "gap_mean_cm-1": zero_mean_gap,
            "gap_rms_sd_cm-1": zero_mean_gap_sd,
            "mixing_mean": zero_mean_mixing,
            "single_branch_rms_sd_cm-1": 0.5 * zero_mean_gap_sd,
        }
    )

    exp_extrema = normalized_extrema(profile_grid, target)
    exp_couplet_extrema = normalized_extrema(profile_grid, target_couplet)
    fit_corrected = next(row for row in ridge_rows if row["J_scale"] == 1.0)

    summary = {
        "status": "numerical_identifiability_study_complete_no_electronic_structure_launched",
        "scope": {
            "manuscripts_modified": False,
            "terachem_jobs_launched": 0,
            "TDDFT_J_frame_join_performed": False,
            "absolute_molar_CD_claimed": False,
        },
        "corrected_coupling_ensemble": {
            "frames": int(len(j_cm)),
            "J_mean_cm-1": float(np.mean(j_cm)),
            "J_sample_sd_cm-1": float(np.std(j_cm, ddof=1)),
            "triple_product_sign_positive_frames": int(
                np.sum(ensemble["triple_product"] > 0.0)
            ),
            "triple_product_sign_negative_frames": int(
                np.sum(ensemble["triple_product"] < 0.0)
            ),
            "archive": str(args.archive.relative_to(ROOT)),
            "archive_sha256": sha256(args.archive),
        },
        "published_constrained_fit_reconstruction": {
            "center_cm-1": EXP_CENTER_CM,
            "latent_component_full_separation_cm-1": EXP_FULL_SPLITTING_CM,
            "low_center_cm-1": EXP_CENTER_CM - EXP_HALF_SPLITTING_CM,
            "high_center_cm-1": EXP_CENTER_CM + EXP_HALF_SPLITTING_CM,
            "low_hwhm_cm-1": EXP_LOW_HWHM_CM,
            "high_hwhm_cm-1": EXP_HIGH_HWHM_CM,
            "low_to_high_absolute_peak_amplitude_ratio": abs(
                EXP_LOW_AMPLITUDE / EXP_HIGH_AMPLITUDE
            ),
            "low_to_high_absolute_integrated_area_ratio": abs(
                EXP_LOW_AMPLITUDE
                * EXP_LOW_HWHM_CM
                / (EXP_HIGH_AMPLITUDE * EXP_HIGH_HWHM_CM)
            ),
            "full_profile_extrema": exp_extrema,
            "couplet_only_extrema": exp_couplet_extrema,
            "note": "Table-S3 reconstruction; raw observations and fitted baseline/slope unavailable.",
        },
        "detuning_required_with_corrected_J": {
            "constrained_gap_cm-1": EXP_FULL_SPLITTING_CM,
            "detuning_mean_cm-1": float(np.mean(required_delta)),
            "detuning_sample_sd_cm-1": float(np.std(required_delta, ddof=1)),
            "detuning_mean_eV": float(np.mean(required_delta) / EV_TO_CM),
            "mixing_mean": float(np.mean(2.0 * np.abs(j_cm) / EXP_FULL_SPLITTING_CM)),
            "minor_site_weight_mean": float(
                np.mean(
                    0.5
                    * (
                        1.0
                        - np.sqrt(
                            1.0 - (2.0 * np.abs(j_cm) / EXP_FULL_SPLITTING_CM) ** 2
                        )
                    )
                )
            ),
            "unconstrained_gap_cm-1": float(unconstrained_gap),
            "unconstrained_required_detuning_mean_cm-1": float(
                np.mean(required_delta_unconstrained)
            ),
            "unconstrained_required_detuning_mean_eV": float(
                np.mean(required_delta_unconstrained) / EV_TO_CM
            ),
            "zero_mean_normal_disorder_alternative": {
                "detuning_signed_mean_cm-1": 0.0,
                "detuning_sd_needed_for_target_mean_gap_cm-1": float(zero_mean_sigma),
                "gap_mean_cm-1": zero_mean_gap,
                "gap_sd_cm-1": zero_mean_gap_sd,
                "single_branch_sd_cm-1": 0.5 * zero_mean_gap_sd,
                "mixing_mean": zero_mean_mixing,
            },
        },
        "profile_fit_with_corrected_J": fit_corrected,
        "STEOM_shift": {
            "TDDFT_reference_eV": TDDFT_REFERENCE_EV,
            "STEOM_reference_eV": STEOM_SITE_EV,
            "common_additive_shift_eV": COMMON_STEOM_SHIFT_EV,
            "common_additive_shift_cm-1": COMMON_STEOM_SHIFT_EV * EV_TO_CM,
            "gap_effect_cm-1": 0.0,
            "common_shift_needed_to_align_STEOM_origin_to_fit_center_cm-1": EXP_CENTER_CM
            - STEOM_SITE_CM,
            "interpretation": "A one-point common additive correction changes only spectral origin. Site-dependent or affine calibration requires additional matched STEOM points.",
        },
        "linewidth_context": {
            "T2_60fs_hwhm_cm-1": hwhm_from_t2(60.0),
            "T2_128fs_hwhm_cm-1": hwhm_from_t2(128.0),
            "effective_T2_from_low_fitted_hwhm_fs": 1.0
            / (2.0 * np.pi * C_CM_PER_S * EXP_LOW_HWHM_CM)
            * 1e15,
            "effective_T2_from_high_fitted_hwhm_fs": 1.0
            / (2.0 * np.pi * C_CM_PER_S * EXP_HIGH_HWHM_CM)
            * 1e15,
            "caution": "Fitted widths include inhomogeneous/multiband effects and must not be interpreted directly as pure-dephasing times.",
        },
        "existing_100fs_TDDFT_diagnostic_not_joined": {
            "mean_gap_eV": 0.1830506208,
            "sample_sd_gap_eV": 0.1109196636,
            "mean_gap_cm-1": 0.1830506208 * EV_TO_CM,
            "sample_sd_gap_cm-1": 0.1109196636 * EV_TO_CM,
            "common_shift_gap_invariance_max_abs_eV": 4.440892098500626e-16,
            "limitations": [
                "100-fs nonequilibrium window, not the 1-ns NVT ensemble",
                "opposite CR2 RESP charge not represented by the corrected production convention",
                "different basis and non-joinable frames",
            ],
        },
        "conclusions": [
            "The constrained component centers determine sqrt(Delta^2+4J^2), not J alone.",
            "Corrected J can reproduce the center separation if a roughly 0.0314-eV detuning is persistent.",
            "That solution is weakly mixed and approximately 98.4% localized per exciton eigenstate.",
            "A common TDDFT-to-STEOM additive shift cannot create or remove site detuning.",
            "Normalized profile fitting has a J/Delta/amplitude ridge; absolute CD and matched framewise site energies are required to break it.",
            "The third experimental band is outside an S1-only two-state model.",
        ],
        "outputs": {
            "scenario_csv": "scenario_sensitivity.csv",
            "profile_ridge_csv": "profile_fit_ridge.csv",
            "center_ridge_csv": "component_center_identifiability_ridge.csv",
            "disorder_csv": "detuning_disorder_sensitivity.csv",
            "profile_csv": "profile_comparison.csv",
            "figure": "splitting_profile_investigation.pdf",
        },
    }

    write_csv(args.output_dir / "scenario_sensitivity.csv", scenario_rows)
    write_csv(args.output_dir / "profile_fit_ridge.csv", ridge_rows)
    write_csv(args.output_dir / "component_center_identifiability_ridge.csv", center_ridge_rows)
    write_csv(args.output_dir / "detuning_disorder_sensitivity.csv", disorder_rows)

    profile_rows = []
    normalized_target = target / np.max(np.abs(target))
    for index, wavenumber in enumerate(profile_grid):
        row: dict[str, float] = {
            "wavenumber_cm-1": float(wavenumber),
            "wavelength_nm": float(1.0e7 / wavenumber),
            "published_fit_reconstruction_normalized": float(normalized_target[index]),
        }
        for name, profile in profiles.items():
            row[name] = float(profile[index] / np.max(np.abs(profile)))
        profile_rows.append(row)
    write_csv(args.output_dir / "profile_comparison.csv", profile_rows)

    with (args.output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(1.0e7 / profile_grid, normalized_target, color="black", lw=2.0, label="Published fit reconstruction")
    for name, color in [
        ("resonant_corrected_J", "#377eb8"),
        ("constrained_fit_gap", "#e41a1c"),
        ("unconstrained_fit_gap", "#4daf4a"),
    ]:
        profile = profiles[name]
        ax.plot(
            1.0e7 / profile_grid,
            profile / np.max(np.abs(profile)),
            lw=1.5,
            color=color,
            label=name.replace("_", " "),
        )
    ax.set_xlim(550, 475)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized relative CD")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Profile is not the same as latent center splitting")

    ax = axes[0, 1]
    ridge_j = np.asarray([row["assumed_J_cm-1"] for row in center_ridge_rows])
    ridge_delta = np.asarray([row["required_detuning_cm-1"] for row in center_ridge_rows])
    ax.plot(ridge_j, ridge_delta, color="#984ea3", lw=2.0)
    ax.scatter([np.mean(j_cm)], [np.mean(required_delta)], color="black", zorder=3)
    ax.annotate("corrected J", (np.mean(j_cm), np.mean(required_delta)), xytext=(45, 225), arrowprops={"arrowstyle": "->"})
    ax.set_xlabel(r"Assumed $J$ (cm$^{-1}$)")
    ax.set_ylabel(r"Required $|\Delta|$ (cm$^{-1}$)")
    ax.set_title(r"Exact ridge: $\Omega^2=\Delta^2+4J^2$")

    ax = axes[1, 0]
    ax.plot(
        [row["J_mean_cm-1"] for row in ridge_rows],
        [row["best_detuning_cm-1"] for row in ridge_rows],
        "o-",
        color="#ff7f00",
        label="best profile detuning",
    )
    ax.axvline(np.mean(j_cm), color="black", ls="--", lw=1.0)
    ax.set_xlabel(r"Mean assumed $J$ (cm$^{-1}$)")
    ax.set_ylabel(r"Best $|\Delta|$ (cm$^{-1}$)")
    ax.set_title("Normalized profile retains a coupling-detuning ridge")

    ax = axes[1, 1]
    persistent_disorder_rows = [
        row for row in disorder_rows if row["distribution"] == "normal_about_required_mean"
    ]
    ax.plot(
        [row["detuning_normal_sd_cm-1"] for row in persistent_disorder_rows],
        [row["gap_rms_sd_cm-1"] for row in persistent_disorder_rows],
        "o-",
        color="#a65628",
        label=r"SD($\Omega$)",
    )
    ax.plot(
        [row["detuning_normal_sd_cm-1"] for row in persistent_disorder_rows],
        [row["single_branch_rms_sd_cm-1"] for row in persistent_disorder_rows],
        "s--",
        color="#f781bf",
        label="single-branch SD",
    )
    ax.set_xlabel(r"Detuning disorder $\sigma_\Delta$ (cm$^{-1}$)")
    ax.set_ylabel(r"Induced spectral SD (cm$^{-1}$)")
    ax.set_title("Differential disorder broadens the inferred branches")
    ax.legend(frameon=False)

    fig.savefig(args.output_dir / "splitting_profile_investigation.pdf")
    fig.savefig(args.output_dir / "splitting_profile_investigation.png", dpi=180)
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

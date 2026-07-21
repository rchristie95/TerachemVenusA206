#!/usr/bin/env python3
r"""Generate Nguyen-style CD panels from the numerical Venus ensemble.

Nguyen et al. plot (i) the two parent CD traces and (ii) their difference,
``TDX - TD``, together with a three-Lorentzian fit.  The numerical data in this
repository contain the interaction-induced tandem CD but not the intrinsic
monomer/TDX magnetic-dipole contribution.  The closest supported analogue is
therefore

* an uncoupled exciton-chirality reference (the same two sites with ``J=0``),
* the coupled tandem exciton contribution evaluated over all 1000 frames, and
* ``uncoupled - coupled``, which has the same subtraction convention as Nguyen.

The Nguyen three-Lorentzian functional form is fitted to that calculated
difference.  The third basis function is retained at the published 481.37-nm
position and width, but only its amplitude is fitted.  A negligible amplitude
therefore records the fact that the present S1-only two-state model does not
predict Nguyen's third experimental band.

No wavelength-density Jacobian is applied: the calculated response ordinate is
re-expressed against photon wavelength in the same way as an experimental CD
trace.  All ordinates are normalized because the available transition-dipole
geometry does not define absolute molar ellipticity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_DATA_DIR = ROOT / "coupling_nvt_production_cr2_1000_20260721"
DEFAULT_VALIDATION = ROOT / "reference" / "orca_validation.json"
C_CM_PER_S = 2.99792458e10

# Nguyen et al. constrained Table S3 third component.  Its amplitude is fitted
# to our numerical spectrum; its centre and width are retained only to test
# whether our S1-only model contains a corresponding feature.
NGUYEN_THIRD_CENTER_CM = 20774.0
NGUYEN_THIRD_HWHM_CM = 565.79


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--validation-json", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    parser.add_argument(
        "--t2-fs",
        type=float,
        default=60.0,
        help="homogeneous pure-dephasing time (default: 60 fs)",
    )
    parser.add_argument(
        "--npts",
        type=int,
        default=8000,
        help="uniform wavenumber-grid points over 400--600 nm",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def portable_path(path: Path) -> str:
    """Record a repository-relative path without leaking workstation layout."""
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return f"external/{resolved.name}"


def load_site_energy(path: Path) -> tuple[float, float, float]:
    with path.open(encoding="utf-8") as handle:
        record = json.load(handle)["orca_steom_bright_state"]
    if not record.get("spectrum_converged", False):
        raise ValueError(f"STEOM spectrum is not marked converged in {path}")
    e0_cm = float(record["wavenumber_cm-1"])
    wavelength_nm = 1.0e7 / e0_cm
    oscillator_strength = float(record["oscillator_strength"])
    if e0_cm <= 0.0 or oscillator_strength <= 0.0:
        raise ValueError(f"Invalid STEOM bright-state record in {path}")
    return e0_cm, wavelength_nm, oscillator_strength


def read_samples(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No coupling samples in {path}")
    required = {"frame", "J_cm"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} lacks columns: {', '.join(sorted(missing))}")
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def load_ensemble(data_dir: Path) -> tuple[dict[str, np.ndarray], Path, Path]:
    csv_path = data_dir / "coupling_samples.csv"
    npz_path = data_dir / "coupling_geometry.npz"
    if not csv_path.is_file() or not npz_path.is_file():
        raise FileNotFoundError(
            f"Expected {csv_path} and {npz_path}; the production ensemble is missing"
        )
    samples = read_samples(csv_path)
    with np.load(npz_path) as archive:
        required = {"frame", "J_cm", "mu_A", "mu_B", "r_A", "r_B"}
        missing = required.difference(archive.files)
        if missing:
            raise ValueError(f"{npz_path} lacks arrays: {', '.join(sorted(missing))}")
        ensemble = {key: np.asarray(archive[key], dtype=float) for key in required}
    n = len(samples["frame"])
    expected = {
        "frame": (n,),
        "J_cm": (n,),
        "mu_A": (n, 3),
        "mu_B": (n, 3),
        "r_A": (n, 3),
        "r_B": (n, 3),
    }
    for key, shape in expected.items():
        if ensemble[key].shape != shape:
            raise ValueError(f"{key} has shape {ensemble[key].shape}; expected {shape}")
        if not np.all(np.isfinite(ensemble[key])):
            raise ValueError(f"Non-finite entries in {npz_path}:{key}")
    if not np.array_equal(ensemble["frame"].astype(int), samples["frame"].astype(int)):
        raise ValueError("Frame numbers differ between coupling CSV and geometry NPZ")
    if not np.allclose(ensemble["J_cm"], samples["J_cm"], atol=1e-8, rtol=0.0):
        raise ValueError("J values differ between coupling CSV and geometry NPZ")
    return ensemble, csv_path, npz_path


def homogeneous_hwhm_cm(t2_fs: float) -> float:
    return 1.0 / (2.0 * np.pi * C_CM_PER_S * t2_fs * 1e-15)


def lorentzian_area(offset_cm: np.ndarray, gamma_cm: float) -> np.ndarray:
    return (gamma_cm / np.pi) / (offset_cm * offset_cm + gamma_cm * gamma_cm)


def numerical_cd_components(
    grid_cm: np.ndarray,
    ensemble: dict[str, np.ndarray],
    e0_cm: float,
    gamma_cm: float,
    chunk_size: int = 128,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return TD high/low-energy components and the J=0 reference CD."""
    mu_a = ensemble["mu_A"]
    mu_b = ensemble["mu_B"]
    displacement = ensemble["r_B"] - ensemble["r_A"]
    base = (np.pi * e0_cm / 2.0) * np.einsum(
        "ij,ij->i", displacement, np.cross(mu_a, mu_b)
    )
    coupling = ensemble["J_cm"]
    td_high = np.zeros_like(grid_cm)
    td_low = np.zeros_like(grid_cm)
    uncoupled = np.zeros_like(grid_cm)
    n = len(coupling)
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        sl = slice(start, stop)
        high_shape = lorentzian_area(
            grid_cm[None, :] - (e0_cm + coupling[sl, None]), gamma_cm
        )
        low_shape = lorentzian_area(
            grid_cm[None, :] - (e0_cm - coupling[sl, None]), gamma_cm
        )
        zero_shape = lorentzian_area(grid_cm[None, :] - e0_cm, gamma_cm)
        td_high += np.sum((-base[sl, None]) * high_shape, axis=0)
        td_low += np.sum((+base[sl, None]) * low_shape, axis=0)
        # At J=0 the equal and opposite exciton-chirality contributions cancel.
        uncoupled += np.sum(
            (-base[sl, None]) * zero_shape + (+base[sl, None]) * zero_shape,
            axis=0,
        )
    return td_high / n, td_low / n, uncoupled / n


def peak_lorentzian(
    grid_cm: np.ndarray, amplitude: float, center_cm: float, hwhm_cm: float
) -> np.ndarray:
    return amplitude / (1.0 + ((grid_cm - center_cm) / hwhm_cm) ** 2)


def fit_nguyen_form(
    grid_cm: np.ndarray,
    difference_norm: np.ndarray,
    e0_cm: float,
    mean_j: float,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """Fit the Table-S3 functional form with a fixed third-band basis."""
    # Parameter order: A_long, A_short, nu_long, nu_short, gamma_long,
    # gamma_short, A3, C, m_scaled.  The slope uses (nu-E0)/1000 so that its
    # fitted coefficient is well scaled numerically.
    p0 = np.asarray(
        [
            +1.0,
            -1.0,
            e0_cm - mean_j,
            e0_cm + mean_j,
            100.0,
            100.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )
    lower = np.asarray(
        [
            0.05,
            -2.0,
            e0_cm - mean_j - 70.0,
            e0_cm + mean_j - 70.0,
            20.0,
            20.0,
            -0.25,
            -0.10,
            -0.10,
        ]
    )
    upper = np.asarray(
        [
            2.0,
            -0.05,
            e0_cm - mean_j + 70.0,
            e0_cm + mean_j + 70.0,
            400.0,
            400.0,
            +0.25,
            +0.10,
            +0.10,
        ]
    )

    wavelength_nm = 1.0e7 / grid_cm
    fit_mask = (wavelength_nm >= 400.0) & (wavelength_nm <= 600.0)
    weights = 1.0 + 8.0 * np.abs(difference_norm)

    def components(parameters: np.ndarray) -> dict[str, np.ndarray]:
        (
            a_long,
            a_short,
            center_long,
            center_short,
            gamma_long,
            gamma_short,
            a3,
            constant,
            slope_scaled,
        ) = parameters
        long_band = peak_lorentzian(grid_cm, a_long, center_long, gamma_long)
        short_band = peak_lorentzian(grid_cm, a_short, center_short, gamma_short)
        third_band = peak_lorentzian(
            grid_cm, a3, NGUYEN_THIRD_CENTER_CM, NGUYEN_THIRD_HWHM_CM
        )
        baseline = constant + slope_scaled * (grid_cm - e0_cm) / 1000.0
        total = long_band + short_band + third_band + baseline
        return {
            "long": long_band,
            "short": short_band,
            "third": third_band,
            "baseline": baseline,
            "total": total,
        }

    def residual(parameters: np.ndarray) -> np.ndarray:
        model = components(parameters)["total"]
        return (model[fit_mask] - difference_norm[fit_mask]) * weights[fit_mask]

    result = least_squares(
        residual,
        p0,
        bounds=(lower, upper),
        x_scale="jac",
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
        max_nfev=10000,
    )
    if not result.success:
        raise RuntimeError(f"Nguyen-form fit failed: {result.message}")
    model_components = components(result.x)
    unweighted = model_components["total"][fit_mask] - difference_norm[fit_mask]
    rms = float(np.sqrt(np.mean(unweighted * unweighted)))
    max_abs = float(np.max(np.abs(unweighted)))
    names = (
        "amplitude_long",
        "amplitude_short",
        "center_long_cm-1",
        "center_short_cm-1",
        "hwhm_long_cm-1",
        "hwhm_short_cm-1",
        "amplitude_third",
        "baseline_constant",
        "baseline_slope_per_1000cm-1",
    )
    fit = {name: float(value) for name, value in zip(names, result.x)}
    fit.update(
        {
            "center_long_nm": 1.0e7 / fit["center_long_cm-1"],
            "center_short_nm": 1.0e7 / fit["center_short_cm-1"],
            "third_center_cm-1_fixed": NGUYEN_THIRD_CENTER_CM,
            "third_center_nm_fixed": 1.0e7 / NGUYEN_THIRD_CENTER_CM,
            "third_hwhm_cm-1_fixed": NGUYEN_THIRD_HWHM_CM,
            "rms_residual_normalized": rms,
            "max_abs_residual_normalized": max_abs,
            "nfev": int(result.nfev),
        }
    )
    return fit, model_components


def wavelength_order(grid_cm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    wavelength_nm = 1.0e7 / grid_cm
    order = np.argsort(wavelength_nm)
    return wavelength_nm[order], order


def mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.9,
            "lines.linewidth": 1.6,
        }
    )
    return plt


def style_nguyen_axis(ax, ylabel: str) -> None:
    ax.axhline(0.0, color="0.35", lw=0.8, ls=":", zorder=0)
    ax.set_xlim(400.0, 600.0)
    ax.set_xticks([400, 440, 480, 520, 560, 600])
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(ylabel)
    ax.tick_params(direction="out", length=4, width=0.8)


def save_figure(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight")


def write_csv(
    path: Path,
    grid_cm: np.ndarray,
    wavelength_nm: np.ndarray,
    order: np.ndarray,
    arrays: dict[str, np.ndarray],
) -> None:
    columns = [
        "wavenumber_cm-1",
        "wavelength_nm",
        "uncoupled_exciton_cd",
        "coupled_td_cd",
        "difference_tdx_minus_td",
        "difference_normalized",
        "fit_total",
        "fit_long_component",
        "fit_short_component",
        "fit_third_component",
        "fit_baseline",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for idx_out, idx_grid in enumerate(order):
            writer.writerow(
                [
                    f"{grid_cm[idx_grid]:.8f}",
                    f"{wavelength_nm[idx_out]:.8f}",
                    *[
                        f"{arrays[key][idx_grid]:.12e}"
                        for key in (
                            "uncoupled",
                            "td",
                            "difference",
                            "difference_norm",
                            "fit_total",
                            "fit_long",
                            "fit_short",
                            "fit_third",
                            "fit_baseline",
                        )
                    ],
                ]
            )


def main() -> None:
    args = parse_args()
    if args.t2_fs <= 0.0 or args.npts < 2000:
        raise ValueError("Require --t2-fs > 0 and --npts >= 2000")
    data_dir = args.data_dir.resolve()
    validation_path = args.validation_json.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ensemble, csv_source, npz_source = load_ensemble(data_dir)
    e0_cm, e0_nm, oscillator_strength = load_site_energy(validation_path)
    gamma_cm = homogeneous_hwhm_cm(args.t2_fs)
    grid_cm = np.linspace(1.0e7 / 600.0, 1.0e7 / 400.0, args.npts)
    td_high, td_low, uncoupled = numerical_cd_components(
        grid_cm, ensemble, e0_cm, gamma_cm
    )
    td = td_high + td_low
    # Nguyen's observable is TDX - TD.  In the supported reduced model the
    # uncoupled TDX exciton-chirality term is zero, so the difference is -TD.
    difference = uncoupled - td
    scale = float(np.max(np.abs(difference)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("Calculated difference CD has zero or invalid amplitude")
    difference_norm = difference / scale
    td_norm = td / scale
    uncoupled_norm = uncoupled / scale

    coupling = ensemble["J_cm"]
    mean_j = float(np.mean(coupling))
    std_j = float(np.std(coupling, ddof=1))
    mean_short_nm = 1.0e7 / (e0_cm + mean_j)
    mean_long_nm = 1.0e7 / (e0_cm - mean_j)
    fit, fit_components = fit_nguyen_form(grid_cm, difference_norm, e0_cm, mean_j)

    wavelength_nm, order = wavelength_order(grid_cm)
    plt = mpl()

    # Nguyen Fig. 4(c) analogue.  These are interaction-induced contributions,
    # not complete absolute TD/TDX molar-ellipticity spectra.
    fig, ax = plt.subplots(figsize=(4.35, 3.55))
    ax.plot(
        wavelength_nm,
        uncoupled_norm[order],
        color="#3659d9",
        label=r"uncoupled reference ($J=0$)",
    )
    ax.plot(
        wavelength_nm,
        td_norm[order],
        color="#e32620",
        label="coupled tandem (1000 frames)",
    )
    style_nguyen_axis(ax, "Exciton-induced CD (normalized)")
    ax.set_ylim(-1.12, 1.12)
    ax.legend(loc="lower left", frameon=False)
    fig.tight_layout()
    traces_stem = out_dir / "Fig_Spectra_NguyenStyle_Traces"
    save_figure(fig, traces_stem)
    plt.close(fig)

    # Nguyen Fig. 4(d) analogue: the numerical difference and the same
    # three-Lorentzian functional form (with a data-fitted third amplitude).
    fig, ax = plt.subplots(figsize=(4.35, 3.55))
    ax.plot(
        wavelength_nm,
        difference_norm[order],
        color="black",
        label=r"calculated TDX$-$TD analogue",
    )
    ax.plot(
        wavelength_nm,
        fit_components["total"][order],
        color="#e32620",
        ls="--",
        label="Nguyen-form 3-Lorentzian fit",
    )
    style_nguyen_axis(ax, r"Calculated $\Delta$CD (normalized)")
    ax.set_ylim(-1.12, 1.12)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    difference_stem = out_dir / "Fig_Spectra_NguyenStyle_Difference"
    save_figure(fig, difference_stem)
    plt.close(fig)

    # Supplementary diagnostic in the style of Nguyen Fig. S3: components and
    # residual.  It makes the absent third S1 band explicit rather than hiding it.
    fig, (ax, residual_ax) = plt.subplots(
        2,
        1,
        figsize=(5.2, 4.4),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.08},
    )
    ax.plot(wavelength_nm, difference_norm[order], color="black", lw=1.2, label="numerical")
    ax.plot(
        wavelength_nm,
        fit_components["long"][order],
        color="#d73027",
        ls=":",
        label="long-wavelength component",
    )
    ax.plot(
        wavelength_nm,
        fit_components["short"][order],
        color="#2855b6",
        ls=":",
        label="short-wavelength component",
    )
    ax.plot(
        wavelength_nm,
        fit_components["third"][order],
        color="#7b3294",
        ls=":",
        label="481-nm test component",
    )
    ax.plot(
        wavelength_nm,
        fit_components["total"][order],
        color="#8b1a1a",
        lw=1.8,
        label="sum",
    )
    ax.axhline(0.0, color="0.35", lw=0.8, ls=":")
    ax.set_ylabel(r"Calculated $\Delta$CD (normalized)")
    ax.set_ylim(-1.12, 1.12)
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    residual = difference_norm - fit_components["total"]
    residual_ax.plot(wavelength_nm, residual[order], color="black", lw=1.0)
    residual_ax.axhline(0.0, color="0.35", lw=0.8, ls=":")
    residual_ax.set_xlim(400.0, 600.0)
    residual_ax.set_xticks([400, 440, 480, 520, 560, 600])
    residual_ax.set_xlabel("Wavelength (nm)")
    residual_ax.set_ylabel("residual")
    max_residual = max(0.01, 1.15 * float(np.max(np.abs(residual))))
    residual_ax.set_ylim(-max_residual, max_residual)
    for axis in (ax, residual_ax):
        axis.tick_params(direction="out", length=4, width=0.8)
    components_stem = out_dir / "Fig_Spectra_NguyenStyle_Components"
    save_figure(fig, components_stem)
    plt.close(fig)

    output_arrays = {
        "uncoupled": uncoupled,
        "td": td,
        "difference": difference,
        "difference_norm": difference_norm,
        "fit_total": fit_components["total"],
        "fit_long": fit_components["long"],
        "fit_short": fit_components["short"],
        "fit_third": fit_components["third"],
        "fit_baseline": fit_components["baseline"],
    }
    # Keep ad-hoc sensitivity runs from overwriting the canonical manuscript
    # audit package.  Only the default manuscript output writes into data_dir;
    # custom --out-dir runs keep their CSV and audit beside their figures.
    canonical_run = out_dir == HERE.resolve()
    record_dir = data_dir if canonical_run else out_dir
    data_csv = record_dir / "spectra_ensemble_nguyen_convention.csv"
    write_csv(data_csv, grid_cm, wavelength_nm, order, output_arrays)

    audit = {
        "description": "Nguyen-style numerical CD reconstruction from the 1000-frame ensemble",
        "observable": "uncoupled exciton-chirality reference minus coupled tandem CD",
        "subtraction_convention": "TDX - TD",
        "limitations": [
            "interaction-induced CD only; intrinsic monomer/TDX CD is unavailable",
            "normalized ordinate; not absolute molar ellipticity",
            "S1-only two-state model; no independently predicted third band",
            "481.37-nm third-component centre and width retained from Nguyen Table S3 only as a zero-signal test basis",
        ],
        "frames": int(len(coupling)),
        "site_energy_cm-1": e0_cm,
        "site_wavelength_nm_from_wavenumber": e0_nm,
        "oscillator_strength": oscillator_strength,
        "T2_star_fs": args.t2_fs,
        "homogeneous_hwhm_cm-1": gamma_cm,
        "J_mean_cm-1": mean_j,
        "J_sample_sd_cm-1": std_j,
        "full_splitting_mean_cm-1": 2.0 * mean_j,
        "full_splitting_sample_sd_cm-1": 2.0 * std_j,
        "mean_short_wavelength_nm": mean_short_nm,
        "mean_long_wavelength_nm": mean_long_nm,
        "mean_wavelength_separation_nm": mean_long_nm - mean_short_nm,
        "fit": fit,
        "sources": {
            "coupling_csv": portable_path(csv_source),
            "coupling_csv_sha256": sha256(csv_source),
            "geometry_npz": portable_path(npz_source),
            "geometry_npz_sha256": sha256(npz_source),
            "orca_validation": portable_path(validation_path),
            "orca_validation_sha256": sha256(validation_path),
        },
        "outputs": {
            "traces_pdf": portable_path(traces_stem.with_suffix(".pdf")),
            "difference_pdf": portable_path(difference_stem.with_suffix(".pdf")),
            "components_pdf": portable_path(components_stem.with_suffix(".pdf")),
            "grid_csv": portable_path(data_csv),
        },
    }
    audit_path = record_dir / "spectra_nguyen_style_audit.json"
    with audit_path.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2)
        handle.write("\n")

    print(f"frames: {len(coupling)}")
    print(f"STEOM site: {e0_cm:.4f} cm^-1 = {e0_nm:.4f} nm; f={oscillator_strength:.9f}")
    print(f"J: {mean_j:.6f} +/- {std_j:.6f} cm^-1")
    print(f"homogeneous HWHM: {gamma_cm:.6f} cm^-1 at T2*={args.t2_fs:.1f} fs")
    print(
        f"mean exciton centers: {mean_short_nm:.4f}, {mean_long_nm:.4f} nm; "
        f"delta-lambda={mean_long_nm - mean_short_nm:.4f} nm"
    )
    print(
        "fit centers: "
        f"{fit['center_short_nm']:.4f}, {fit['center_long_nm']:.4f} nm; "
        f"HWHM={fit['hwhm_short_cm-1']:.3f}, {fit['hwhm_long_cm-1']:.3f} cm^-1"
    )
    print(
        f"third amplitude at {fit['third_center_nm_fixed']:.4f} nm: "
        f"{fit['amplitude_third']:.6e}; RMS residual={fit['rms_residual_normalized']:.6e}"
    )
    print(f"wrote: {traces_stem.with_suffix('.pdf').name}")
    print(f"wrote: {difference_stem.with_suffix('.pdf').name}")
    print(f"wrote: {components_stem.with_suffix('.pdf').name}")
    print(f"wrote: {data_csv}")
    print(f"wrote: {audit_path}")


if __name__ == "__main__":
    main()

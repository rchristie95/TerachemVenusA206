#!/usr/bin/env python3
r"""Generate the three ensemble-spectroscopy panels used in Figure 5.

The production mutual transition-density couplings and placed transition-dipole
geometries are read from the 1000-frame NVT output.  The script therefore uses
the complete sampled ensemble rather than hard-coded coupling or geometry
values.  Excitation axes are plotted as photon wavelength.

The CD panel also evaluates the constrained three-Lorentzian dVenus model from
Table S3 of Nguyen et al., Biophysical Journal (2025),
doi:10.1016/j.bpj.2025.10.022.  That published trace is the difference molar
ellipticity ``dVenus-TDX minus dVenus-TD``.  Its independently normalized overlay
is consequently a comparison of band positions and splitting only, not of
absolute amplitude or sign with the calculated tandem CD spectrum.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = HERE.parent / "coupling_nvt_production_cr2_1000_20260721"
C_CM_PER_S = 2.99792458e10

# Nguyen et al. Table S3 constrained dVenus fit (wavenumbers in cm^-1).
NGUYEN_G1 = 19322.0
NGUYEN_DELTA = 130.79
NGUYEN_A1 = +3.59
NGUYEN_GAMMA1 = 293.51
NGUYEN_A2 = -2.27
NGUYEN_GAMMA2 = 439.25
NGUYEN_PEAK3 = 20774.0
NGUYEN_A3 = -0.34
NGUYEN_GAMMA3 = 565.79


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="production coupling directory (default: 1000-frame CR2 run)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=HERE,
        help="output directory for panel PDFs/PNGs (default: manuscript directory)",
    )
    parser.add_argument(
        "--e0-nm",
        type=float,
        default=515.0,
        help="monomer absorption center in nm (default: 515)",
    )
    parser.add_argument(
        "--t2-fs",
        type=float,
        default=60.0,
        help="homogeneous pure-dephasing time in fs (default: 60)",
    )
    parser.add_argument(
        "--npts",
        type=int,
        default=6000,
        help="number of wavenumber-grid points (default: 6000)",
    )
    return parser.parse_args()


def read_sample_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No samples found in {path}")
    required = {"frame", "J_cm", "separation_A", "angle_deg"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def load_ensemble(data_dir: Path) -> dict[str, np.ndarray]:
    csv_path = data_dir / "coupling_samples.csv"
    geometry_path = data_dir / "coupling_geometry.npz"
    missing = [str(path) for path in (csv_path, geometry_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Production coupling data are not available yet:\n  "
            + "\n  ".join(missing)
            + "\nRun again after the 1000-frame GPU coupling calculation has completed, "
            "or supply --data-dir."
        )

    samples = read_sample_csv(csv_path)
    with np.load(geometry_path) as archive:
        required = {"frame", "J_cm", "mu_A", "mu_B", "r_A", "r_B"}
        missing_keys = required.difference(archive.files)
        if missing_keys:
            raise ValueError(
                f"{geometry_path} is missing arrays: {', '.join(sorted(missing_keys))}"
            )
        ensemble = {key: np.asarray(archive[key]) for key in required}

    n = len(samples["frame"])
    expected_shapes = {
        "frame": (n,),
        "J_cm": (n,),
        "mu_A": (n, 3),
        "mu_B": (n, 3),
        "r_A": (n, 3),
        "r_B": (n, 3),
    }
    for key, shape in expected_shapes.items():
        if ensemble[key].shape != shape:
            raise ValueError(f"{key} has shape {ensemble[key].shape}; expected {shape}")
        if not np.all(np.isfinite(ensemble[key])):
            raise ValueError(f"Non-finite values in {geometry_path}:{key}")

    if not np.array_equal(ensemble["frame"].astype(int), samples["frame"].astype(int)):
        raise ValueError("Frame numbers disagree between the coupling CSV and geometry NPZ")
    if not np.allclose(ensemble["J_cm"], samples["J_cm"], atol=1e-8, rtol=0.0):
        raise ValueError("Couplings disagree between the coupling CSV and geometry NPZ")

    ensemble["separation_A"] = samples["separation_A"]
    ensemble["angle_deg"] = samples["angle_deg"]
    ensemble["source_csv"] = np.asarray(str(csv_path))
    ensemble["source_geometry"] = np.asarray(str(geometry_path))
    return ensemble


def homogeneous_hwhm_cm(t2_fs: float) -> float:
    return 1.0 / (2.0 * np.pi * C_CM_PER_S * t2_fs * 1e-15)


def lorentzian(offset_cm: np.ndarray, gamma_cm: float) -> np.ndarray:
    return (gamma_cm / np.pi) / (offset_cm * offset_cm + gamma_cm * gamma_cm)


def band_properties(ensemble: dict[str, np.ndarray], e0_cm: float):
    mu_a = ensemble["mu_A"]
    mu_b = ensemble["mu_B"]
    displacement = ensemble["r_B"] - ensemble["r_A"]
    mu_plus = (mu_a + mu_b) / np.sqrt(2.0)
    mu_minus = (mu_a - mu_b) / np.sqrt(2.0)
    intensity_plus = np.einsum("ij,ij->i", mu_plus, mu_plus)
    intensity_minus = np.einsum("ij,ij->i", mu_minus, mu_minus)
    scalar_triple = np.einsum("ij,ij->i", displacement, np.cross(mu_a, mu_b))
    base = (np.pi * e0_cm / 2.0) * scalar_triple
    coupling = ensemble["J_cm"]
    return {
        "nu_plus": e0_cm + coupling,
        "nu_minus": e0_cm - coupling,
        "intensity_plus": intensity_plus,
        "intensity_minus": intensity_minus,
        "rotation_plus": -base,
        "rotation_minus": +base,
    }


def ensemble_spectra(
    grid_cm: np.ndarray,
    properties: dict[str, np.ndarray],
    gamma_cm: float,
    chunk_size: int = 128,
) -> tuple[np.ndarray, np.ndarray]:
    absorption = np.zeros_like(grid_cm)
    cd = np.zeros_like(grid_cm)
    n = len(properties["nu_plus"])
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        sl = slice(start, stop)
        plus_shape = lorentzian(
            grid_cm[None, :] - properties["nu_plus"][sl, None], gamma_cm
        )
        minus_shape = lorentzian(
            grid_cm[None, :] - properties["nu_minus"][sl, None], gamma_cm
        )
        absorption += np.sum(
            properties["intensity_plus"][sl, None] * plus_shape
            + properties["intensity_minus"][sl, None] * minus_shape,
            axis=0,
        )
        cd += np.sum(
            properties["rotation_plus"][sl, None] * plus_shape
            + properties["rotation_minus"][sl, None] * minus_shape,
            axis=0,
        )
    return absorption / n, cd / n


def wavenumber_density_to_wavelength(
    grid_cm: np.ndarray, signal_per_cm: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Transform a spectral density from reciprocal centimeters to per nm."""
    wavelength_nm = 1.0e7 / grid_cm
    order = np.argsort(wavelength_nm)
    wavelength_nm = wavelength_nm[order]
    jacobian = 1.0e7 / wavelength_nm**2
    return wavelength_nm, signal_per_cm[order] * jacobian


def nguyen_triple_lorentzian(grid_cm: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return Table S3 components and zero-baseline total for constrained dVenus."""
    center1 = NGUYEN_G1 - NGUYEN_DELTA
    center2 = NGUYEN_G1 + NGUYEN_DELTA
    first = NGUYEN_A1 / (1.0 + ((grid_cm - center1) / NGUYEN_GAMMA1) ** 2)
    second = NGUYEN_A2 / (1.0 + ((grid_cm - center2) / NGUYEN_GAMMA2) ** 2)
    third = NGUYEN_A3 / (1.0 + ((grid_cm - NGUYEN_PEAK3) / NGUYEN_GAMMA3) ** 2)
    return first, second, third, first + second + third


def map_signal_to_wavelength(
    grid_cm: np.ndarray, signal: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Re-express an observed ordinate versus wavelength (no density Jacobian)."""
    wavelength_nm = 1.0e7 / grid_cm
    order = np.argsort(wavelength_nm)
    return wavelength_nm[order], signal[order]


def normalize(signal: np.ndarray) -> np.ndarray:
    scale = float(np.max(np.abs(signal)))
    return signal / scale if scale > 0.0 else signal.copy()


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
            "legend.fontsize": 8.3,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.formatter.use_mathtext": True,
        }
    )
    return plt


C_ABS = "#1f6aa5"
C_POS = "#c0392b"
C_NEG = "#1f6aa5"
C_EXP = "#3a7d44"
C_EXP_LIGHT = "#72a879"
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


def save_panel(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")


def main() -> None:
    args = parse_args()
    if args.e0_nm <= 0.0 or args.t2_fs <= 0.0 or args.npts < 1000:
        raise ValueError("Require positive --e0-nm/--t2-fs and --npts >= 1000")

    data_dir = args.data_dir.resolve()
    ensemble = load_ensemble(data_dir)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    coupling = ensemble["J_cm"]
    mean_j = float(np.mean(coupling))
    nframes = len(coupling)
    std_j = float(np.std(coupling, ddof=1)) if nframes > 1 else 0.0
    mean_sep = float(np.mean(ensemble["separation_A"]))
    std_sep = (
        float(np.std(ensemble["separation_A"], ddof=1)) if nframes > 1 else 0.0
    )
    mean_angle = float(np.mean(ensemble["angle_deg"]))
    std_angle = float(np.std(ensemble["angle_deg"], ddof=1)) if nframes > 1 else 0.0
    e0_cm = 1.0e7 / args.e0_nm
    gamma_cm = homogeneous_hwhm_cm(args.t2_fs)

    # Wide enough to retain Nguyen's third band near 481 nm.
    grid_cm = np.linspace(1.0e7 / 540.0, 1.0e7 / 460.0, args.npts)
    properties = band_properties(ensemble, e0_cm)
    absorption_cm, cd_cm = ensemble_spectra(grid_cm, properties, gamma_cm)
    wavelength_nm, absorption_nm = wavenumber_density_to_wavelength(
        grid_cm, absorption_cm
    )
    wavelength_cd_nm, cd_nm = wavenumber_density_to_wavelength(grid_cm, cd_cm)

    nguyen_components = nguyen_triple_lorentzian(grid_cm)
    nguyen_wavelength_nm, nguyen_total = map_signal_to_wavelength(
        grid_cm, nguyen_components[-1]
    )

    plt = mpl()

    # Panel (a): homogeneous linewidth versus the sampled Davydov splitting.
    splitting = 2.0 * mean_j
    splitting_sd = 2.0 * std_j
    experimental_splitting = 2.0 * NGUYEN_DELTA
    t2 = np.linspace(20.0, 250.0, 500)
    fwhm = 1.0 / (np.pi * C_CM_PER_S * t2 * 1e-15)
    t_cross = 1.0 / (np.pi * C_CM_PER_S * splitting) * 1e15
    f_at_t2 = 1.0 / (np.pi * C_CM_PER_S * args.t2_fs * 1e-15)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axvspan(t_cross, t2.max(), color=C_EXP, alpha=0.10)
    ax.plot(t2, fwhm, color=C_ABS, lw=2.2, label=r"homog. FWHM $1/(\pi cT_2^{*})$")
    ax.axhspan(
        splitting - splitting_sd,
        splitting + splitting_sd,
        color=C_POS,
        alpha=0.12,
        lw=0,
    )
    ax.axhline(
        splitting,
        color=C_POS,
        lw=1.8,
        ls="--",
        label=fr"ensemble $2|J|={splitting:.0f}\pm{splitting_sd:.0f}$ cm$^{{-1}}$",
    )
    ax.axhline(
        experimental_splitting,
        color=C_EXP,
        lw=1.6,
        ls="-.",
        label=fr"Nguyen constrained $2\delta={experimental_splitting:.0f}$ cm$^{{-1}}$",
    )
    ax.axvline(args.t2_fs, color=C_MEAN, lw=1.0, ls=":")
    ax.plot([args.t2_fs], [f_at_t2], marker="o", ms=5.5, color=C_MEAN, zorder=5)
    ax.annotate(
        fr"$T_2^*={args.t2_fs:.0f}$ fs",
        xy=(args.t2_fs, f_at_t2),
        xytext=(args.t2_fs + 18.0, f_at_t2 - 40.0),
        fontsize=9,
        color=C_MEAN,
        arrowprops=dict(arrowstyle="-", color=C_MEAN, lw=0.8),
    )
    ax.text(
        0.97,
        0.28,
        "resolved region",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=8.5,
        color=C_EXP,
    )
    ax.set_xlim(20.0, 250.0)
    ax.set_ylim(0.0, max(400.0, experimental_splitting + 50.0))
    ax.set_xlabel(r"$T_2^{*}$ (fs)")
    ax.set_ylabel(r"linewidth or splitting (cm$^{-1}$)")
    style(ax)
    leg(ax, loc="upper right")
    fig.tight_layout()
    resolvability_stem = out_dir / "Fig_Spectra_Resolvability"
    save_panel(fig, resolvability_stem)
    plt.close(fig)

    # Panel (b): ensemble absorption in photon wavelength.
    fig, ax = plt.subplots(figsize=FIGSIZE)
    absorption_norm = absorption_nm / (float(np.max(absorption_nm)) + 1e-30)
    ax.plot(wavelength_nm, absorption_norm, color=C_ABS, lw=2.0)
    theory_short_nm = 1.0e7 / (e0_cm + mean_j)
    theory_long_nm = 1.0e7 / (e0_cm - mean_j)
    for center in (theory_short_nm, theory_long_nm):
        ax.axvline(center, color="0.50", ls=":", lw=1.1)
    ax.annotate(
        "",
        xy=(theory_short_nm, 1.04),
        xytext=(theory_long_nm, 1.04),
        arrowprops=dict(arrowstyle="<->", color=C_MEAN, lw=1.2),
    )
    ax.text(
        (theory_short_nm + theory_long_nm) / 2.0,
        1.09,
        fr"$\Delta\lambda={theory_long_nm - theory_short_nm:.2f}$ nm",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    ax.set_xlim(500.0, 530.0)
    ax.set_ylim(0.0, 1.22)
    ax.set_xlabel("photon wavelength (nm)")
    ax.set_ylabel("absorption (normalized)")
    style(ax)
    fig.tight_layout()
    absorption_stem = out_dir / "Fig_Spectra_Absorption"
    save_panel(fig, absorption_stem)
    plt.close(fig)

    # Panel (c): calculated tandem CD and the published constrained dVenus
    # difference-molar-ellipticity model.  Separate normalization is deliberate:
    # the overlay compares only band positions and splitting.
    fig, ax = plt.subplots(figsize=FIGSIZE)
    cd_norm = normalize(cd_nm)
    nguyen_norm = normalize(nguyen_total)
    ax.fill_between(
        wavelength_cd_nm,
        cd_norm,
        0.0,
        where=(cd_norm >= 0.0),
        color=C_POS,
        alpha=0.24,
        lw=0,
    )
    ax.fill_between(
        wavelength_cd_nm,
        cd_norm,
        0.0,
        where=(cd_norm < 0.0),
        color=C_NEG,
        alpha=0.24,
        lw=0,
    )
    ax.plot(
        wavelength_cd_nm,
        cd_norm,
        color=C_MEAN,
        lw=1.5,
        label="calculated tandem CD",
    )
    ax.plot(
        nguyen_wavelength_nm,
        nguyen_norm,
        color=C_EXP,
        lw=1.9,
        ls="--",
        label="Nguyen constrained 3-Lorentzian fit",
    )
    experimental_long_nm = 1.0e7 / (NGUYEN_G1 - NGUYEN_DELTA)
    experimental_short_nm = 1.0e7 / (NGUYEN_G1 + NGUYEN_DELTA)
    third_nm = 1.0e7 / NGUYEN_PEAK3
    for center in (experimental_short_nm, experimental_long_nm):
        ax.axvline(center, color=C_EXP_LIGHT, ls=":", lw=1.0)
    ax.axvline(third_nm, color=C_EXP_LIGHT, ls=":", lw=0.8, alpha=0.65)
    ax.axhline(0.0, color="k", lw=0.7)
    ax.text(
        0.03,
        0.05,
        "positions/splitting only;\nindependent normalization",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.0,
        color="0.25",
    )
    ax.set_xlim(470.0, 535.0)
    ax.set_xlabel("photon wavelength (nm)")
    ax.set_ylabel("normalized signal (separate scales)")
    style(ax)
    # Keep the long comparison labels outside the data region so neither peak is hidden.
    leg(ax, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    fig.tight_layout()
    cd_stem = out_dir / "Fig_Spectra_CD"
    save_panel(fig, cd_stem)
    plt.close(fig)

    print(f"source CSV: {ensemble['source_csv'].item()}")
    print(f"source geometry: {ensemble['source_geometry'].item()}")
    print(f"frames: {len(coupling)}")
    print(f"J: {mean_j:.6f} +/- {std_j:.6f} cm^-1 (sample SD)")
    print(f"separation: {mean_sep:.6f} +/- {std_sep:.6f} A")
    print(f"inter-dipole angle: {mean_angle:.6f} +/- {std_angle:.6f} degrees")
    print(
        f"calculated exciton centers: {theory_short_nm:.4f} and "
        f"{theory_long_nm:.4f} nm; delta-lambda={theory_long_nm-theory_short_nm:.4f} nm"
    )
    print(
        f"Nguyen constrained centers: {experimental_short_nm:.4f} and "
        f"{experimental_long_nm:.4f} nm; delta-lambda="
        f"{experimental_long_nm-experimental_short_nm:.4f} nm; "
        f"third band={third_nm:.4f} nm"
    )
    print(
        "IMPORTANT: Nguyen's published signal is dVenus-TDX minus dVenus-TD "
        "difference molar ellipticity. The normalized overlay compares positions and "
        "splitting only, not absolute amplitude or sign."
    )
    for stem in (resolvability_stem, absorption_stem, cd_stem):
        print(f"wrote {stem.with_suffix('.pdf').name} and {stem.with_suffix('.png').name}")


if __name__ == "__main__":
    main()

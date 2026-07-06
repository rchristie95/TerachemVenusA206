#!/usr/bin/env python3
r"""
absorption_cd_spectra.py  --  Excitonic absorption + circular-dichroism (CD) lineshapes.

Reviewer item 2 (raised by three referees): the open-quantum-systems section is
"illustrative" and disconnected from the quantum chemistry, and no spectroscopic
signature is ever computed. This script closes that gap by turning the computed
Davydov coupling J (and its fluctuations) plus the dephasing into an absorption
and a CD lineshape, and overlaying the experimental Davydov splitting window.

Experimental window: Nguyen et al. infer an apparent excitonic COUPLING
U = 131-186 cm^-1 from the dVenus tandem-dimer CD spectrum. Under their
Delta E = 2U convention this is a Davydov SPLITTING of 262-372 cm^-1, which is
the like-for-like comparison for our splitting 2|J|. (The earlier default that
shaded 131-186 as a splitting compared U against 2|J| and is fixed here.)

Physics (degenerate excitonic dimer, |1> = |e1 g2>, |2> = |g1 e2>):
  - Eigenstates  |+-> = (|1> +- |2>)/sqrt(2)  at energies  nu_+- = E0 +- J.
  - Transition dipoles  mu_+- = (mu_A +- mu_B)/sqrt(2);  dipole strength D = |mu|^2.
  - Rotational strength of the exciton couplet (Rosenfeld / DeVoe exciton theory):
        R_+- = -+ (pi * nu0 / 2) * R_AB . (mu_A x mu_B)
    giving a conservative bisignate CD couplet whose sign flips with the sign of J.
  - Broadening:
        * homogeneous  (Lorentzian HWHM)  gamma_cm = 1 / (2 pi c T2*)   from dephasing
        * inhomogeneous (Gaussian sigma)   sigma_cm = std(J)             from the MD
          coupling distribution (coupling_ensemble.py -> coupling_distribution.json)
    combined as a Voigt profile per band. The two reviewer additions reinforce
    each other: the spread in J from item 1 IS the inhomogeneous broadening here.

Geometry of the two transition dipoles:
  - Best (used for the manuscript figure): pass --geometry-json with explicit
    mu_A, mu_B (3-vectors) and r_A, r_B (positions, Angstrom), exported from the
    STEOM transition density placed at both chromophore sites by
    export_dipole_geometry.py (same Kabsch/`super` placement as the coupling
    pipeline; gives |r_A - r_B| ~ 25.5 A, inter-dipole angle ~ 104 deg).
  - Default (schematic, clearly labelled): built from --separation, --angle,
    --skew (out-of-plane chirality angle) and --dipole-debye, seeded with the
    STEOM values (25.5 A, 104 deg). Use this only for a quick look.

Outputs (in --out, default `lineshape_out/`), styled to match the paper figures:
  Fig_Spectra_Coupling.pdf    panel (a): J distribution over the NVT ensemble
  Fig_Spectra_Absorption.pdf  panel (b): Davydov doublet, bands at +/- J, 2|J|
  Fig_Spectra_CD.pdf          panel (c): bisignate CD couplet + exp. window
  Fig_Spectra.pdf             the three panels composed (preview)
  lineshape_data.csv          the raw absorption/CD grid
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import voigt_profile

# Speed of light in cm/s (for the dephasing -> linewidth conversion).
C_CM_PER_S = 2.99792458e10


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def build_default_geometry(separation_a, angle_deg, skew_deg, dipole_debye):
    """
    Schematic but adjustable dimer geometry that yields a non-zero CD couplet.

    Inter-monomer axis is along z; mu_A lies along x; mu_B is rotated by the
    inter-dipole angle and tilted out of plane by `skew` so that
    R_AB . (mu_A x mu_B) != 0 (a coplanar arrangement would give zero CD).
    """
    theta = np.radians(angle_deg)
    chi = np.radians(skew_deg)
    m = dipole_debye
    mu_A = m * np.array([1.0, 0.0, 0.0])
    mu_B = m * np.array([np.cos(theta),
                         np.sin(theta) * np.cos(chi),
                         np.sin(theta) * np.sin(chi)])
    r_A = np.array([0.0, 0.0, -separation_a / 2.0])
    r_B = np.array([0.0, 0.0, +separation_a / 2.0])
    return mu_A, mu_B, r_A, r_B


def load_geometry_json(path):
    with open(path) as f:
        g = json.load(f)
    mu_A = np.asarray(g["mu_A"], float)
    mu_B = np.asarray(g["mu_B"], float)
    r_A = np.asarray(g["r_A"], float)
    r_B = np.asarray(g["r_B"], float)
    return mu_A, mu_B, r_A, r_B


# --------------------------------------------------------------------------- #
# Exciton band parameters
# --------------------------------------------------------------------------- #
def exciton_bands(E0, J, mu_A, mu_B, r_A, r_B):
    """
    Return per-band (energy, dipole strength, rotational strength) for the
    symmetric (+) and antisymmetric (-) excitons.
    """
    mu_plus = (mu_A + mu_B) / np.sqrt(2.0)
    mu_minus = (mu_A - mu_B) / np.sqrt(2.0)
    D_plus = float(np.dot(mu_plus, mu_plus))
    D_minus = float(np.dot(mu_minus, mu_minus))

    R_AB = r_B - r_A
    cross = np.cross(mu_A, mu_B)
    # Conservative couplet: R_+ = -(pi nu0 / 2) R_AB.(muA x muB); R_- = +that.
    base = (np.pi * E0 / 2.0) * float(np.dot(R_AB, cross))
    R_plus = -base
    R_minus = +base

    return {
        "plus": {"nu": E0 + J, "D": D_plus, "R": R_plus},
        "minus": {"nu": E0 - J, "D": D_minus, "R": R_minus},
    }


def homogeneous_hwhm_cm(t2_star_fs):
    """Lorentzian HWHM (cm^-1) from a pure-dephasing time T2* (fs)."""
    t2_s = t2_star_fs * 1e-15
    return 1.0 / (2.0 * np.pi * C_CM_PER_S * t2_s)


# --------------------------------------------------------------------------- #
# Spectra
# --------------------------------------------------------------------------- #
def build_spectra(grid_cm, bands, sigma_cm, gamma_cm):
    """Sum Voigt-broadened bands into absorption and CD spectra (arb. units)."""
    absorption = np.zeros_like(grid_cm)
    cd = np.zeros_like(grid_cm)
    # voigt_profile needs sigma>0; clamp a tiny floor for the pure-Lorentzian case.
    sig = max(sigma_cm, 1e-6)
    for band in bands.values():
        shape = voigt_profile(grid_cm - band["nu"], sig, gamma_cm)
        absorption += band["D"] * shape
        cd += band["R"] * shape
    return absorption, cd


def load_ensemble_geometry(path):
    """Per-frame geometry npz from coupling_ensemble.py: mu_A, mu_B, r_A, r_B, J_cm."""
    d = np.load(path)
    return {k: np.asarray(d[k], float) for k in ("mu_A", "mu_B", "r_A", "r_B", "J_cm")}


def build_spectra_ensemble(grid_cm, geo, E0, gamma_cm, sigma_h_cm=0.0):
    """
    Absorption + CD as an explicit sum over the MD ensemble (hard data): for every
    frame, place the two Davydov bands at E0 +/- J(frame) with dipole strengths and
    rotational strengths from THAT frame's actual STEOM dipoles/centroids, then
    broaden each band by the HOMOGENEOUS width only (Lorentzian gamma from T2*,
    plus an optional small Gaussian sigma_h). The inhomogeneous envelope and the
    CD couplet amplitude/sign then emerge from the real sampled disorder rather
    than an assumed Gaussian applied to a single geometry.

    Returns (absorption, cd, diag) with per-ensemble mean band diagnostics.
    """
    mu_A, mu_B, r_A, r_B, J = geo["mu_A"], geo["mu_B"], geo["r_A"], geo["r_B"], geo["J_cm"]
    n = len(J)
    absorption = np.zeros_like(grid_cm)
    cd = np.zeros_like(grid_cm)
    sig = max(sigma_h_cm, 1e-6)
    D_plus = D_minus = R_couplet = 0.0
    for i in range(n):
        mp = (mu_A[i] + mu_B[i]) / np.sqrt(2.0)
        mm = (mu_A[i] - mu_B[i]) / np.sqrt(2.0)
        dP = float(np.dot(mp, mp))
        dM = float(np.dot(mm, mm))
        base = (np.pi * E0 / 2.0) * float(np.dot(r_B[i] - r_A[i], np.cross(mu_A[i], mu_B[i])))
        nu_p, nu_m = E0 + J[i], E0 - J[i]
        sh_p = voigt_profile(grid_cm - nu_p, sig, gamma_cm)
        sh_m = voigt_profile(grid_cm - nu_m, sig, gamma_cm)
        absorption += dP * sh_p + dM * sh_m
        cd += (-base) * sh_p + (+base) * sh_m
        D_plus += dP; D_minus += dM; R_couplet += base
    absorption /= n
    cd /= n
    diag = {"n": n, "D_plus": D_plus / n, "D_minus": D_minus / n,
            "R_couplet": R_couplet / n, "J_mean": float(np.mean(J)),
            "J_std": float(np.std(J, ddof=1)) if n > 1 else 0.0}
    return absorption, cd, diag


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--J", type=float, default=96.38,
                   help="Davydov coupling J (cm^-1); STEOM thermal mean. "
                        "Overridden by --distribution mean when given.")
    p.add_argument("--E0", type=float, default=18437.0, help="Monomer site energy (cm^-1).")
    p.add_argument("--t2-star-fs", type=float, default=60.0, help="Pure-dephasing time T2* (fs).")
    p.add_argument("--distribution", type=Path, default=None,
                   help="coupling_distribution.json from coupling_ensemble.py "
                        "(sets J = mean and the inhomogeneous Gaussian width = std).")
    p.add_argument("--sigma-cm", type=float, default=None,
                   help="Override inhomogeneous Gaussian std (cm^-1).")
    p.add_argument("--geometry-json", type=Path, default=None,
                   help="JSON with mu_A, mu_B, r_A, r_B for a SINGLE geometry (from the QM density).")
    p.add_argument("--ensemble-geometry", type=Path, default=None,
                   help="coupling_geometry.npz (per-frame mu/r/J from coupling_ensemble.py): "
                        "build absorption+CD as an explicit sum over the MD ensemble, so the "
                        "broadening and CD couplet are grounded in the real sampled disorder.")
    p.add_argument("--hist-binwidth", type=float, default=0.4,
                   help="Fixed histogram bin width for panel (a) (cm^-1); keeps the bin width "
                        "constant as the sample count grows so the histogram just converges.")
    p.add_argument("--separation", type=float, default=25.5,
                   help="Centroid separation (Angstrom); STEOM value (schematic geometry only).")
    p.add_argument("--angle", type=float, default=104.0,
                   help="Inter-dipole angle (deg); STEOM value (schematic geometry only).")
    p.add_argument("--skew", type=float, default=45.0, help="Out-of-plane chirality angle (deg).")
    p.add_argument("--dipole-debye", type=float, default=9.6,
                   help="Monomer transition-dipole magnitude (D); STEOM |mu| (schematic geometry only).")
    p.add_argument("--exp-splitting", type=float, nargs=2, default=[262.0, 372.0],
                   metavar=("LO", "HI"),
                   help="Experimental Davydov splitting window (cm^-1); Nguyen U=131-186 -> 2U=262-372.")
    p.add_argument("--window", type=float, default=900.0, help="Half-width of the energy axis (cm^-1).")
    p.add_argument("--npts", type=int, default=4000, help="Energy-grid resolution.")
    p.add_argument("--out", type=Path, default=Path("lineshape_out"))
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    J = args.J
    sigma_cm = args.sigma_cm if args.sigma_cm is not None else 0.0
    samples = None
    if args.distribution is not None and args.distribution.exists():
        with open(args.distribution) as f:
            dist = json.load(f)
        J = float(dist.get("mean", J))
        if args.sigma_cm is None:
            sigma_cm = float(dist.get("std", 0.0))
        samples = np.asarray(dist.get("samples", []), float)
        print(f"[*] From {args.distribution}: J(mean)={J:.2f} cm^-1, std={sigma_cm:.2f} cm^-1, "
              f"n={samples.size} samples")

    gamma_cm = homogeneous_hwhm_cm(args.t2_star_fs)
    print(f"[*] Homogeneous HWHM from T2*={args.t2_star_fs:.1f} fs: {gamma_cm:.2f} cm^-1")
    print(f"[*] Inhomogeneous Gaussian sigma: {sigma_cm:.2f} cm^-1")

    grid = np.linspace(args.E0 - args.window, args.E0 + args.window, args.npts)

    if args.ensemble_geometry is not None and args.ensemble_geometry.exists():
        # HARD-DATA path: absorption/CD summed over the real MD ensemble.
        geo = load_ensemble_geometry(args.ensemble_geometry)
        absorption, cd, diag = build_spectra_ensemble(grid, geo, args.E0, gamma_cm)
        J = diag["J_mean"]
        if samples is None:
            samples = geo["J_cm"]
        print(f"[*] Ensemble lineshape over {diag['n']} MD frames "
              f"({args.ensemble_geometry.name}):")
        print(f"    - <D(+)>={diag['D_plus']:.3f}  <D(-)>={diag['D_minus']:.3f}  "
              f"<R couplet>={diag['R_couplet']:.3e}")
        print(f"    - J = {diag['J_mean']:.2f} +/- {diag['J_std']:.2f} cm^-1, "
              f"2|J| = {2*abs(J):.2f} cm^-1  (inhomogeneous width from real disorder)")
    else:
        if args.geometry_json is not None:
            mu_A, mu_B, r_A, r_B = load_geometry_json(args.geometry_json)
            geom_note = f"geometry-json ({args.geometry_json.name})"
        else:
            mu_A, mu_B, r_A, r_B = build_default_geometry(
                args.separation, args.angle, args.skew, args.dipole_debye)
            geom_note = "schematic default geometry"
        print(f"[*] Using {geom_note} (single geometry + Gaussian inhomogeneous width)")

        bands = exciton_bands(args.E0, J, mu_A, mu_B, r_A, r_B)
        print(f"    - band(+): nu={bands['plus']['nu']:.1f}  D={bands['plus']['D']:.3f}  R={bands['plus']['R']:.3e}")
        print(f"    - band(-): nu={bands['minus']['nu']:.1f}  D={bands['minus']['D']:.3f}  R={bands['minus']['R']:.3e}")
        print(f"    - computed Davydov splitting 2|J| = {2*abs(J):.2f} cm^-1")
        absorption, cd = build_spectra(grid, bands, sigma_cm, gamma_cm)

    # ----- write data -----
    csv_path = args.out / "lineshape_data.csv"
    with open(csv_path, "w") as f:
        f.write("wavenumber_cm,absorption,cd\n")
        for x, a, c in zip(grid, absorption, cd):
            f.write(f"{x:.4f},{a:.8e},{c:.8e}\n")

    # ----- plots -----
    paths = _plot(args, grid, absorption, cd, J, sigma_cm, samples)

    two_j = 2 * abs(J)
    lo, hi = args.exp_splitting
    within = lo <= two_j <= hi
    rel = "within" if within else ("below" if two_j < lo else "above")
    print("\n" + "=" * 60)
    print(f"  computed 2|J| = {two_j:.1f} cm^-1   "
          f"experiment = {lo:.0f}-{hi:.0f} cm^-1  ({rel} the window)")
    print("=" * 60)
    print(f"  data        : {csv_path}")
    for k, v in paths.items():
        print(f"  {k:11s} : {v}")


# --------------------------------------------------------------------------- #
# Plotting (styled to match the paper's other matplotlib figures)
# --------------------------------------------------------------------------- #
# Palette shared with open_quantum_dynamics.py.
_C_ABS = "#1f6aa5"   # absorption / distribution (blue)
_C_POS = "#c0392b"   # CD positive lobe (red)
_C_NEG = "#1f6aa5"   # CD negative lobe (blue)
_C_EXP = "#3a7d44"   # experimental window (green)
_C_MEAN = "#222222"


def _mpl():
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


def _style(ax):
    ax.grid(alpha=0.25, lw=0.6)
    for s in ax.spines.values():
        s.set_linewidth(0.8)


def _leg(ax, **kw):
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False, **kw)


FIGSIZE = (4.2, 3.4)


def _hist_bins(samples, binwidth):
    """Fixed-width bin edges spanning the samples (constant bin width as n grows)."""
    lo = np.floor(samples.min() / binwidth) * binwidth
    hi = np.ceil(samples.max() / binwidth) * binwidth + binwidth
    return np.arange(lo, hi, binwidth)


def _panel_coupling(plt, out, samples, J, sigma_cm, binwidth=0.4):
    """Panel (a): J distribution over the NVT ensemble."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    if samples is not None and samples.size:
        ax.hist(samples, bins=_hist_bins(samples, binwidth), color=_C_ABS,
                alpha=0.80, edgecolor="white", lw=0.4)
        ax.set_ylabel("MD snapshots")
    else:
        # Fall back to the implied Gaussian if no per-frame samples are present.
        xs = np.linspace(J - 4 * sigma_cm, J + 4 * sigma_cm, 400)
        g = np.exp(-0.5 * ((xs - J) / max(sigma_cm, 1e-6)) ** 2)
        ax.plot(xs, g, color=_C_ABS, lw=2)
        ax.set_ylabel("density (arb.)")
    ax.axvline(J, color=_C_MEAN, lw=1.8, label=fr"$\bar J = {J:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.axvspan(J - sigma_cm, J + sigma_cm, color=_C_MEAN, alpha=0.10,
               label=fr"$\pm\sigma = {sigma_cm:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.set_xlabel(r"Davydov coupling $J$ (cm$^{-1}$)")
    _style(ax)
    _leg(ax, loc="upper right")
    fig.tight_layout()
    path = out / "Fig_Spectra_Coupling.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def _panel_absorption(plt, out, rel, absorption, J):
    """Panel (b): excitonic absorption Davydov doublet."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    y = absorption / (absorption.max() + 1e-30)
    ax.plot(rel, y, color=_C_ABS, lw=2)
    for x in (-J, +J):
        ax.axvline(x, color="0.55", ls=":", lw=1.1)
    # 2|J| splitting annotation between the two bands.
    ax.annotate("", xy=(-J, 1.04), xytext=(+J, 1.04),
                arrowprops=dict(arrowstyle="<->", color=_C_MEAN, lw=1.2))
    ax.text(0.0, 1.09, fr"$2|J| = {2*abs(J):.0f}\,\mathrm{{cm^{{-1}}}}$",
            ha="center", va="bottom", fontsize=10)
    ax.set_ylim(top=1.20)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel("absorption (norm.)")
    _style(ax)
    fig.tight_layout()
    path = out / "Fig_Spectra_Absorption.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def _panel_cd(plt, out, rel, cd, exp_splitting):
    """Panel (c): bisignate CD couplet with the experimental splitting window."""
    lo, hi = exp_splitting
    fig, ax = plt.subplots(figsize=FIGSIZE)
    y = cd / (np.max(np.abs(cd)) + 1e-30)
    # Colour the two lobes of the conservative couplet.
    ax.fill_between(rel, y, 0, where=(y >= 0), color=_C_POS, alpha=0.30, lw=0)
    ax.fill_between(rel, y, 0, where=(y < 0), color=_C_NEG, alpha=0.30, lw=0)
    ax.plot(rel, y, color=_C_MEAN, lw=1.6, label="computed CD")
    ax.axhline(0, color="k", lw=0.8)
    # Experimental Davydov splitting: bands sit at +/- (splitting/2); shade that window.
    for sgn in (-1, +1):
        ax.axvspan(sgn * lo / 2, sgn * hi / 2, color=_C_EXP, alpha=0.20,
                   label=fr"exp. $\Delta E/2$" if sgn == 1 else None)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel(r"$\Delta\varepsilon$ (norm.)")
    _style(ax)
    _leg(ax, loc="upper right")
    fig.tight_layout()
    path = out / "Fig_Spectra_CD.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def _plot(args, grid, absorption, cd, J, sigma_cm, samples):
    plt = _mpl()
    rel = grid - args.E0
    bw = args.hist_binwidth
    p_a = _panel_coupling(plt, args.out, samples, J, sigma_cm, binwidth=bw)
    p_b = _panel_absorption(plt, args.out, rel, absorption, J)
    p_c = _panel_cd(plt, args.out, rel, cd, args.exp_splitting)

    # Composed 3-panel preview (the manuscript uses the three PDFs as subfigures).
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.4))
    for ax, sub in zip(axes, ("a", "b", "c")):
        ax.set_title(f"({sub})", loc="left", fontsize=12)
    # (a)
    if samples is not None and samples.size:
        axes[0].hist(samples, bins=_hist_bins(samples, bw), color=_C_ABS,
                     alpha=0.80, edgecolor="white", lw=0.4)
        axes[0].set_ylabel("MD snapshots")
    axes[0].axvline(J, color=_C_MEAN, lw=1.8, label=fr"$\bar J = {J:.0f}$")
    axes[0].axvspan(J - sigma_cm, J + sigma_cm, color=_C_MEAN, alpha=0.10,
                    label=fr"$\pm\sigma = {sigma_cm:.0f}$")
    axes[0].set_xlabel(r"$J$ (cm$^{-1}$)")
    _style(axes[0]); _leg(axes[0], loc="upper right")
    # (b)
    yb = absorption / (absorption.max() + 1e-30)
    axes[1].plot(rel, yb, color=_C_ABS, lw=2)
    for x in (-J, +J):
        axes[1].axvline(x, color="0.55", ls=":", lw=1.1)
    axes[1].annotate("", xy=(-J, 1.04), xytext=(+J, 1.04),
                     arrowprops=dict(arrowstyle="<->", color=_C_MEAN, lw=1.2))
    axes[1].text(0.0, 1.09, fr"$2|J|={2*abs(J):.0f}$", ha="center", va="bottom", fontsize=10)
    axes[1].set_ylim(top=1.20)
    axes[1].set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)"); axes[1].set_ylabel("abs. (norm.)")
    _style(axes[1])
    # (c)
    lo, hi = args.exp_splitting
    yc = cd / (np.max(np.abs(cd)) + 1e-30)
    axes[2].fill_between(rel, yc, 0, where=(yc >= 0), color=_C_POS, alpha=0.30, lw=0)
    axes[2].fill_between(rel, yc, 0, where=(yc < 0), color=_C_NEG, alpha=0.30, lw=0)
    axes[2].plot(rel, yc, color=_C_MEAN, lw=1.6)
    axes[2].axhline(0, color="k", lw=0.8)
    for sgn in (-1, +1):
        axes[2].axvspan(sgn * lo / 2, sgn * hi / 2, color=_C_EXP, alpha=0.20,
                        label=fr"exp. $\Delta E/2$" if sgn == 1 else None)
    axes[2].set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)"); axes[2].set_ylabel(r"$\Delta\varepsilon$ (norm.)")
    _style(axes[2]); _leg(axes[2], loc="upper right")
    fig.tight_layout()
    composed = args.out / "Fig_Spectra.pdf"
    fig.savefig(composed)
    plt.close(fig)

    return {"coupling(a)": p_a, "absorption(b)": p_b, "CD(c)": p_c, "composed": composed}


if __name__ == "__main__":
    main()

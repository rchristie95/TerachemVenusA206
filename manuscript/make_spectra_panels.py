#!/usr/bin/env python3
r"""
Generate the three excitonic-spectroscopy panel PDFs (Figure 5):
  (a) Fig_Spectra_Resolvability.pdf -- homogeneous linewidth vs Davydov splitting vs T2*
  (b) Fig_Spectra_Absorption.pdf    -- excitonic absorption Davydov doublet at +/- J
  (c) Fig_Spectra_CD.pdf            -- bisignate CD couplet with the experimental splitting window

Panel (a) is NOT the coupling distribution (that is Figure 4c, Fig_Tandem_Histogram.pdf, and
re-plotting it here would duplicate it). Instead it shows spectroscopic resolvability: the
homogeneous Lorentzian FWHM 1/(pi c T2*) crosses below the fixed Davydov splitting 2|J| once
T2* exceeds ~48 fs, so the exciton doublet/couplet is spectrally resolvable across the
physically relevant dephasing range -- and the adopted 60 fs sits just inside that regime.

This is a lightweight, self-contained companion to the released `absorption_cd_spectra.py`
that renders the manuscript panels directly from the reported ensemble statistics. If the
real ensemble outputs are present (coupling_distribution.json / coupling_geometry.npz), pass
them to `absorption_cd_spectra.py` instead for the production figure.

Numbers are the manuscript values:
  J = 111 +/- 3 cm^-1 (tandem NVT ensemble)   E0 = 515 nm (experimental absorption max)
  T2* = 60 fs (homogeneous pure-dephasing)     geometry: 25.5 A separation, 104 deg inter-dipole
  experimental Davydov splitting 2U = 262-372 cm^-1 (Nguyen et al.)
"""

from pathlib import Path
import numpy as np
from scipy.special import voigt_profile

OUT_DIR = Path(__file__).parent
C_CM_PER_S = 2.99792458e10

# ---- manuscript parameters ----
J_MEAN = 111.0          # cm^-1, tandem ensemble mean
J_STD = 3.0             # cm^-1, tandem ensemble std (inhomogeneous width)
T2_STAR_FS = 60.0       # fs, homogeneous pure-dephasing time
E0_CM = 1.0e7 / 515.0   # cm^-1, experimental absorption max (~19417 cm^-1)
SEP_A = 25.5            # Angstrom, chromophore centroid separation
ANGLE_DEG = 104.0       # inter-dipole angle
SKEW_DEG = 35.0         # out-of-plane chirality angle (representative; sign/|R| need placed geometry)
DIPOLE_D = 9.6          # Debye, monomer transition-dipole magnitude
EXP_SPLIT = (262.0, 372.0)   # cm^-1, experimental Davydov splitting window (Nguyen 2U)
WINDOW = 500.0          # cm^-1, half-width of energy axis
NPTS = 4000


def homogeneous_hwhm_cm(t2_fs):
    return 1.0 / (2.0 * np.pi * C_CM_PER_S * t2_fs * 1e-15)


def geometry(sep_a, angle_deg, skew_deg, dip_d):
    th, ch = np.radians(angle_deg), np.radians(skew_deg)
    mu_A = dip_d * np.array([1.0, 0.0, 0.0])
    mu_B = dip_d * np.array([np.cos(th), np.sin(th) * np.cos(ch), np.sin(th) * np.sin(ch)])
    r_A = np.array([0.0, 0.0, -sep_a / 2.0])
    r_B = np.array([0.0, 0.0, +sep_a / 2.0])
    return mu_A, mu_B, r_A, r_B


def bands(E0, J, mu_A, mu_B, r_A, r_B):
    mu_p, mu_m = (mu_A + mu_B) / np.sqrt(2), (mu_A - mu_B) / np.sqrt(2)
    base = (np.pi * E0 / 2.0) * float(np.dot(r_B - r_A, np.cross(mu_A, mu_B)))
    return {"plus": {"nu": E0 + J, "D": float(mu_p @ mu_p), "R": -base},
            "minus": {"nu": E0 - J, "D": float(mu_m @ mu_m), "R": +base}}


def spectra(grid, bnds, sigma_cm, gamma_cm):
    absn = np.zeros_like(grid)
    cd = np.zeros_like(grid)
    sig = max(sigma_cm, 1e-6)
    for b in bnds.values():
        shape = voigt_profile(grid - b["nu"], sig, gamma_cm)
        absn += b["D"] * shape
        cd += b["R"] * shape
    return absn, cd


def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.labelsize": 12, "axes.titlesize": 12, "legend.fontsize": 9,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.formatter.use_mathtext": True,
    })
    return plt


C_ABS = "#1f6aa5"
C_POS = "#c0392b"
C_NEG = "#1f6aa5"
C_EXP = "#3a7d44"
C_EXP_DK = "#2f6b39"
C_MEAN = "#222222"
FIGSIZE = (4.2, 3.4)


def style(ax):
    ax.grid(alpha=0.25, lw=0.6)
    for s in ax.spines.values():
        s.set_linewidth(0.8)


def leg(ax, **kw):
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False, **kw)


def main():
    plt = mpl()
    gamma = homogeneous_hwhm_cm(T2_STAR_FS)
    print(f"[*] homogeneous HWHM (T2*={T2_STAR_FS:.0f} fs) = {gamma:.1f} cm^-1 "
          f"(FWHM {2*gamma:.0f} cm^-1); inhomogeneous sigma = {J_STD:.0f} cm^-1")
    grid = np.linspace(E0_CM - WINDOW, E0_CM + WINDOW, NPTS)
    rel = grid - E0_CM

    mu_A, mu_B, r_A, r_B = geometry(SEP_A, ANGLE_DEG, SKEW_DEG, DIPOLE_D)
    bnds = bands(E0_CM, J_MEAN, mu_A, mu_B, r_A, r_B)
    absn, cd = spectra(grid, bnds, J_STD, gamma)

    # --- Panel (a): spectroscopic resolvability of the Davydov splitting ---
    splitting = 2.0 * J_MEAN                                        # cm^-1
    t2 = np.linspace(20.0, 250.0, 500)                             # fs
    fwhm = 1.0 / (np.pi * C_CM_PER_S * t2 * 1e-15)                 # cm^-1, homogeneous FWHM
    t_cross = 1.0 / (np.pi * C_CM_PER_S * splitting) * 1e15        # fs where FWHM = 2|J|
    f60 = 1.0 / (np.pi * C_CM_PER_S * T2_STAR_FS * 1e-15)          # cm^-1 at T2* = 60 fs
    print(f"[*] resolvability: FWHM = 2|J| = {splitting:.0f} cm^-1 at T2* = {t_cross:.0f} fs; "
          f"FWHM(60 fs) = {f60:.0f} cm^-1")
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axvspan(t_cross, 250.0, color=C_EXP, alpha=0.12)
    ax.text(0.97, 0.06, "splitting\nresolved", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=9, color=C_EXP_DK)
    ax.plot(t2, fwhm, color=C_ABS, lw=2.2, label=r"homog. FWHM $1/(\pi c T_2^{*})$")
    ax.axhline(splitting, color=C_POS, lw=1.8, ls="--",
               label=fr"$2|J| = {splitting:.0f}\,\mathrm{{cm^{{-1}}}}$")
    ax.axvline(T2_STAR_FS, color=C_MEAN, lw=1.0, ls=":")
    ax.plot([T2_STAR_FS], [f60], marker="o", ms=6, color=C_MEAN, zorder=5)
    ax.annotate(r"$T_2^{*}=60$ fs", xy=(T2_STAR_FS, f60), xytext=(T2_STAR_FS + 10, f60 + 55),
                fontsize=9, color=C_MEAN)
    ax.set_xlim(20, 250)
    ax.set_ylim(0, 400)
    ax.set_xlabel(r"$T_2^{*}$ (fs)")
    ax.set_ylabel(r"linewidth, splitting (cm$^{-1}$)")
    style(ax); leg(ax, loc="upper right")
    fig.tight_layout()
    p_r = OUT_DIR / "Fig_Spectra_Resolvability.pdf"
    fig.savefig(p_r); plt.close(fig); print(f"  wrote {p_r.name}")

    # --- Panel (b): absorption Davydov doublet ---
    fig, ax = plt.subplots(figsize=FIGSIZE)
    y = absn / (absn.max() + 1e-30)
    ax.plot(rel, y, color=C_ABS, lw=2)
    for x in (-J_MEAN, +J_MEAN):
        ax.axvline(x, color="0.55", ls=":", lw=1.1)
    ax.annotate("", xy=(-J_MEAN, 1.04), xytext=(+J_MEAN, 1.04),
                arrowprops=dict(arrowstyle="<->", color=C_MEAN, lw=1.2))
    ax.text(0.0, 1.09, fr"$2|J| = {2*J_MEAN:.0f}\,\mathrm{{cm^{{-1}}}}$",
            ha="center", va="bottom", fontsize=10)
    ax.set_ylim(top=1.22)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel("absorption (norm.)")
    style(ax)
    fig.tight_layout()
    p_a = OUT_DIR / "Fig_Spectra_Absorption.pdf"
    fig.savefig(p_a); plt.close(fig); print(f"  wrote {p_a.name}")

    # --- Panel (c): bisignate CD couplet ---
    fig, ax = plt.subplots(figsize=FIGSIZE)
    yc = cd / (np.max(np.abs(cd)) + 1e-30)
    ax.fill_between(rel, yc, 0, where=(yc >= 0), color=C_POS, alpha=0.30, lw=0)
    ax.fill_between(rel, yc, 0, where=(yc < 0), color=C_NEG, alpha=0.30, lw=0)
    ax.plot(rel, yc, color=C_MEAN, lw=1.6, label="computed CD")
    ax.axhline(0, color="k", lw=0.8)
    lo, hi = EXP_SPLIT
    for sgn in (-1, +1):
        ax.axvspan(sgn * lo / 2, sgn * hi / 2, color=C_EXP, alpha=0.20,
                   label=r"exp. $\Delta E/2$" if sgn == 1 else None)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel(r"$\Delta\varepsilon$ (norm.)")
    style(ax); leg(ax, loc="upper right")
    fig.tight_layout()
    p_b = OUT_DIR / "Fig_Spectra_CD.pdf"
    fig.savefig(p_b); plt.close(fig); print(f"  wrote {p_b.name}")

    two_j = 2 * J_MEAN
    within = lo <= two_j <= hi
    rel_txt = "within" if within else ("below" if two_j < lo else "above")
    print(f"\n  computed 2|J| = {two_j:.0f} cm^-1 ; experiment {lo:.0f}-{hi:.0f} cm^-1 "
          f"({rel_txt} the window)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""
absorption_cd_spectra.py  --  Excitonic absorption + circular-dichroism (CD) lineshapes.

Reviewer item 2 (raised by three referees): the open-quantum-systems section is
"illustrative" and disconnected from the quantum chemistry, and no spectroscopic
signature is ever computed. This script closes that gap by turning the computed
Davydov coupling J (and its fluctuations) plus the dephasing into an absorption
and a CD lineshape, and overlaying the experimental Davydov splitting window
(Nguyen et al.: 131-186 cm^-1 from the dVenus tandem-dimer CD spectrum).

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
  - Best: pass --geometry-json with explicit mu_A, mu_B (3-vectors) and r_A, r_B
    (positions, Angstrom), e.g. exported from the QM transition density.
  - Default (schematic, clearly labelled): built from --separation, --angle
    (inter-dipole angle, default 92.85 deg from the paper), --skew (out-of-plane
    chirality angle) and --dipole-debye. Use this for a quick look; quote the
    geometry-json result in the manuscript.

Outputs (in --out, default `lineshape_out/`):
  Fig_Absorption_Spectrum.pdf, Fig_CD_Spectrum.pdf, lineshape_data.csv
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--J", type=float, default=74.38, help="Davydov coupling J (cm^-1).")
    p.add_argument("--E0", type=float, default=18437.0, help="Monomer site energy (cm^-1).")
    p.add_argument("--t2-star-fs", type=float, default=60.0, help="Pure-dephasing time T2* (fs).")
    p.add_argument("--distribution", type=Path, default=None,
                   help="coupling_distribution.json from coupling_ensemble.py "
                        "(sets J = mean and the inhomogeneous Gaussian width = std).")
    p.add_argument("--sigma-cm", type=float, default=None,
                   help="Override inhomogeneous Gaussian std (cm^-1).")
    p.add_argument("--geometry-json", type=Path, default=None,
                   help="JSON with mu_A, mu_B, r_A, r_B (preferred; from the QM density).")
    p.add_argument("--separation", type=float, default=27.6, help="Centroid separation (Angstrom).")
    p.add_argument("--angle", type=float, default=92.85, help="Inter-dipole angle (deg).")
    p.add_argument("--skew", type=float, default=45.0, help="Out-of-plane chirality angle (deg).")
    p.add_argument("--dipole-debye", type=float, default=10.0, help="Monomer transition-dipole magnitude.")
    p.add_argument("--exp-splitting", type=float, nargs=2, default=[131.0, 186.0],
                   metavar=("LO", "HI"), help="Experimental Davydov splitting window (cm^-1).")
    p.add_argument("--window", type=float, default=900.0, help="Half-width of the energy axis (cm^-1).")
    p.add_argument("--npts", type=int, default=4000, help="Energy-grid resolution.")
    p.add_argument("--out", type=Path, default=Path("lineshape_out"))
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    J = args.J
    sigma_cm = args.sigma_cm if args.sigma_cm is not None else 0.0
    if args.distribution is not None and args.distribution.exists():
        with open(args.distribution) as f:
            dist = json.load(f)
        J = float(dist.get("mean", J))
        if args.sigma_cm is None:
            sigma_cm = float(dist.get("std", 0.0))
        print(f"[*] From {args.distribution}: J(mean)={J:.2f} cm^-1, std={sigma_cm:.2f} cm^-1")

    gamma_cm = homogeneous_hwhm_cm(args.t2_star_fs)
    print(f"[*] Homogeneous HWHM from T2*={args.t2_star_fs:.1f} fs: {gamma_cm:.2f} cm^-1")
    print(f"[*] Inhomogeneous Gaussian sigma: {sigma_cm:.2f} cm^-1")

    if args.geometry_json is not None:
        mu_A, mu_B, r_A, r_B = load_geometry_json(args.geometry_json)
        geom_note = f"geometry-json ({args.geometry_json.name})"
    else:
        mu_A, mu_B, r_A, r_B = build_default_geometry(
            args.separation, args.angle, args.skew, args.dipole_debye)
        geom_note = "schematic default geometry"
    print(f"[*] Using {geom_note}")

    bands = exciton_bands(args.E0, J, mu_A, mu_B, r_A, r_B)
    print(f"    - band(+): nu={bands['plus']['nu']:.1f}  D={bands['plus']['D']:.3f}  R={bands['plus']['R']:.3e}")
    print(f"    - band(-): nu={bands['minus']['nu']:.1f}  D={bands['minus']['D']:.3f}  R={bands['minus']['R']:.3e}")
    print(f"    - computed Davydov splitting 2|J| = {2*abs(J):.2f} cm^-1")

    grid = np.linspace(args.E0 - args.window, args.E0 + args.window, args.npts)
    absorption, cd = build_spectra(grid, bands, sigma_cm, gamma_cm)

    # ----- write data -----
    csv_path = args.out / "lineshape_data.csv"
    with open(csv_path, "w") as f:
        f.write("wavenumber_cm,absorption,cd\n")
        for x, a, c in zip(grid, absorption, cd):
            f.write(f"{x:.4f},{a:.8e},{c:.8e}\n")

    # ----- plots -----
    abs_pdf, cd_pdf = _plot(args, grid, absorption, cd, J, gamma_cm, sigma_cm)

    print("\n" + "=" * 52)
    print(f"  computed 2|J| = {2*abs(J):.1f} cm^-1   "
          f"experiment = {args.exp_splitting[0]:.0f}-{args.exp_splitting[1]:.0f} cm^-1")
    print("=" * 52)
    print(f"  data        : {csv_path}")
    print(f"  absorption  : {abs_pdf}")
    print(f"  CD          : {cd_pdf}")


def _plot(args, grid, absorption, cd, J, gamma_cm, sigma_cm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rel = grid - args.E0
    lo, hi = args.exp_splitting

    # Absorption
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(rel, absorption / (absorption.max() + 1e-30), color="#4C72B0", lw=2)
    ax.axvline(-J, color="0.5", ls=":", lw=1)
    ax.axvline(+J, color="0.5", ls=":", lw=1)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel("Absorption (norm.)")
    ax.set_title(fr"Excitonic absorption ($2|J|={2*abs(J):.0f}$ cm$^{{-1}}$)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    abs_pdf = args.out / "Fig_Absorption_Spectrum.pdf"
    fig.savefig(abs_pdf)
    plt.close(fig)

    # CD with experimental splitting overlay
    fig, ax = plt.subplots(figsize=(6, 4.2))
    norm = np.max(np.abs(cd)) + 1e-30
    ax.plot(rel, cd / norm, color="#C44E52", lw=2, label="computed CD")
    ax.axhline(0, color="k", lw=0.8)
    # Experimental Davydov splitting window: shade |nu-E0| in [lo/2, hi/2] on both sides.
    for sgn in (-1, +1):
        ax.axvspan(sgn * lo / 2, sgn * hi / 2, color="#55A868", alpha=0.18,
                   label="exp. splitting/2" if sgn == 1 else None)
    ax.set_xlabel(r"$\nu - E_0$ (cm$^{-1}$)")
    ax.set_ylabel(r"$\Delta\varepsilon$ (norm.)")
    ax.set_title(r"Excitonic CD couplet vs experimental splitting")
    ax.legend(fontsize=9, frameon=True, loc="best")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    cd_pdf = args.out / "Fig_CD_Spectrum.pdf"
    fig.savefig(cd_pdf)
    plt.close(fig)
    return abs_pdf, cd_pdf


if __name__ == "__main__":
    main()

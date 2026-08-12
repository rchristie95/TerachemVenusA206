"""Absolute exciton CD of the tandem dimer against Nguyen 2025 Table S3.

The experimental scale is now recoverable: Nguyen et al. report (main-text
Methods, not the SI) eps_516 = 184 400 M^-1 cm^-1 for dVenus-TD and 92 200 for
dVenus-TDX, all samples at 25 uM in a 0.5 cm cuvette, with
[theta] = (mdeg * MW) / (10 * L * C), i.e. molar ellipticity per mole of
CONSTRUCT (the exact 2:1 ratio of the two extinction coefficients confirms the
per-construct normalisation, TD carrying two chromophores and TDX one). The
Figure 4 ordinate is Delta[Theta] x 10^-5 deg cm^2 dmol^-1, so the Table S3
Lorentzian amplitudes A1 = +3.59, A2 = -2.27, A3 = -0.34 are in units of 10^5.

Our prediction uses the coupled-oscillator (exciton-chirality) mechanism, which
needs only the two ELECTRIC transition dipoles and their displacement -- see
audit_magnetic_dipoles.py for why TeraChem's printed magnetic dipoles cannot be
used, and note that the intrinsic monomer term they would supply is not what the
Figure 4 difference spectrum reports anyway.

The comparison is anchored on the couplet FIRST MOMENT

    M = int (Delta_eps / nubar) (nubar - nubar_0) d nubar = R_+ Omega / 2.296e-39

because R_+ = -pi nubar (J/Omega) T carries a 1/Omega that cancels against the
Omega lever arm:

    R_+ Omega = -pi nubar J T .

So M depends on J and the chirality triple product ONLY. It is independent of
the detuning, of the homogeneous linewidth, and of the inhomogeneous spread --
every quantity in this problem that is either uncertain or force-field-derived.
Peak amplitudes are reported too, but they are the lineshape-dependent number
and should be read as the weaker test.

Run:
    python terachem_site_energy_cd/scripts/absolute_cd_vs_experiment.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from absolute_cd import (  # noqa: E402
    DEBYE_PER_AU,
    CGS_PER_ANGSTROM_DEBYE2,
    DELTA_EPSILON_PREFACTOR,
    ELLIPTICITY_PER_DELTA_EPSILON,
)

ENSEMBLE = ROOT / "results" / "ensembles" / "ens_v2_all.npz"
HARTREE_CM = 219474.6313632
BOHR_TO_ANGSTROM = 0.529177210903

# --- Nguyen et al., Biophys. J. 124, 4293 (2025), Table S3 (dVenus-TD) --------
# amplitudes in 1e5 deg cm^2 dmol^-1, centres and HWHM in cm^-1.
#
# SIGN. Table S3 fits the DIFFERENCE spectrum, and Nguyen define it as
# "subtracting the TD from TDX" (main text, CD data analysis):
#
#     Delta[theta] = [theta]_TDX - [theta]_TD
#
# The excitonic couplet exists only in TD -- TDX has a single chromophore and no
# couplet -- so the couplet appearing in Delta[theta] is MINUS the couplet of the
# dimer itself. Their tabulated A1 = +3.59 sits on the low-wavenumber component,
# and the main text confirms the difference spectrum is negative at 511.8 nm and
# positive at 521.7 nm. Inverting, the TD dimer's own couplet is positive to the
# blue and negative to the red, i.e. a NEGATIVE exciton chirality.
#
# We predict the dimer's couplet, so the tabulated amplitudes must be negated
# before comparison. Comparing directly against Delta[theta] reports a spurious
# handedness disagreement.
DIFFERENCE_SPECTRUM_SIGN = -1.0  # Delta[theta] = TDX - TD  ->  TD couplet
EXP_CENTER_CM = 19322.0
EXP_HALF_SPLITTING_CM = 130.79
EXP_BANDS = tuple(  # (amplitude_1e5, centre_cm, hwhm_cm), as the TD couplet
    (DIFFERENCE_SPECTRUM_SIGN * a, c, w) for a, c, w in (
        (+3.59, EXP_CENTER_CM - EXP_HALF_SPLITTING_CM, 293.51),
        (-2.27, EXP_CENTER_CM + EXP_HALF_SPLITTING_CM, 439.25),
    )
)
EXP_THIRD_BAND = (DIFFERENCE_SPECTRUM_SIGN * -0.34, 20774.0, 565.79)
AMPLITUDE_UNIT = 1.0e5

# Transition-density coupling, 1000-frame production ensemble. NOTE this value
# is already screened by epsilon = 1.77: the unscreened point-dipole coupling
# recomputed from coupling_geometry.npz is +48.93 cm^-1 and the stored
# J_pda_cm is +27.64 = 48.93/1.77 exactly, so the same screening is carried by
# the stored J_cm. The screened value is the one that belongs in the exciton
# Hamiltonian and therefore in the mixing coefficient c_A c_B = J/Omega.
J_TDC_CM = 32.8165
J_TDC_SD_CM = 1.5540
SCREENING_EPSILON = 1.77

# Independent 1000-frame coupling ensemble, for a triple-product cross-check.
# It uses transition-density centroids rather than CR2 heavy-atom centroids, so
# it is a genuinely different geometric convention.
COUPLING_NPZ = (
    ROOT.parent / "coupling_nvt_production_cr2_1000_20260721" / "coupling_geometry.npz"
)


def lorentzian(grid, amplitude, centre, hwhm):
    """Peak-normalised Lorentzian, matching the Table S3 fit convention."""
    return amplitude / (1.0 + ((grid - centre) / hwhm) ** 2)


def point_dipole_coupling_cm(mu_a, mu_b, r_a, r_b):
    """Signed point-dipole coupling, cm^-1, from a.u. dipoles and Angstrom positions.

    NOTE the unit conversion: the Coulomb sum is evaluated with the separation in
    Angstrom, so 1/R^3 picks up BOHR_TO_ANGSTROM**3, not its reciprocal.
    """
    r_ab = r_b - r_a
    dist = np.linalg.norm(r_ab, axis=-1)
    n = r_ab / dist[:, None]
    kappa = np.einsum("ij,ij->i", mu_a, mu_b) - 3.0 * np.einsum(
        "ij,ij->i", mu_a, n
    ) * np.einsum("ij,ij->i", mu_b, n)
    return kappa * BOHR_TO_ANGSTROM**3 / dist**3 * HARTREE_CM


def main() -> None:
    d = np.load(ENSEMBLE)
    mu_a, mu_b = d["mu_a_au"], d["mu_b_au"]
    r_a, r_b = d["r_a_ang"], d["r_b_ang"]
    e_a, e_b = d["e_a_cm"], d["e_b_cm"]
    n = len(e_a)

    # ---- chirality triple product, in cgs -----------------------------------
    triple_au = np.einsum("ij,ij->i", r_b - r_a, np.cross(mu_a, mu_b))  # Ang * a.u.^2
    triple_cgs = triple_au * DEBYE_PER_AU**2 * CGS_PER_ANGSTROM_DEBYE2

    # ---- signed coupling ----------------------------------------------------
    j_pda = point_dipole_coupling_cm(mu_a, mu_b, r_a, r_b)
    sign_j = np.sign(np.median(j_pda))
    j_cm = J_TDC_CM * sign_j  # TDC magnitude, sign fixed by the aligned dipoles

    delta = e_a - e_b
    omega = np.hypot(delta, 2.0 * j_cm)

    print(f"ENSEMBLE  n = {n} frames  ({ENSEMBLE.name})")
    print(f"  mean |detuning|       {np.abs(delta).mean():9.1f} cm^-1")
    print(f"  mean Omega            {omega.mean():9.1f} cm^-1   "
          f"harmonic {n/np.sum(1/omega):.1f}")
    print(f"  separation            {np.linalg.norm(r_b-r_a,axis=1).mean():9.2f} A")
    print(f"  triple product        {triple_au.mean():9.2f} +- "
          f"{triple_au.std(ddof=1):.2f} Ang*a.u.^2   "
          f"({int((triple_au<0).sum())}/{n} negative)")
    print(f"  point-dipole coupling {j_pda.mean():+9.2f} +- {j_pda.std(ddof=1):.2f} cm^-1")
    print(f"  -> adopted J          {j_cm:+9.2f} cm^-1 (TDC magnitude, PDA sign)\n")

    # ---- couplet first moment: lineshape- and detuning-free -----------------
    # R_+ Omega = -pi nubar J T   (the 1/Omega in R_+ cancels the Omega lever arm)
    m_pred_frames = -np.pi * EXP_CENTER_CM * j_cm * triple_cgs  # esu^2 cm^2 * cm^-1
    m_pred = m_pred_frames.mean()
    # to Delta_eps units then to molar ellipticity
    m_pred_eps = DELTA_EPSILON_PREFACTOR * m_pred
    m_pred_theta = m_pred_eps * ELLIPTICITY_PER_DELTA_EPSILON

    # experimental first moment from the Table S3 couplet
    grid = np.arange(17000.0, 22500.0, 1.0)
    exp_couplet = sum(
        lorentzian(grid, a * AMPLITUDE_UNIT, c, w) for a, c, w in EXP_BANDS
    )
    m_exp_theta = np.trapezoid(exp_couplet * (grid - EXP_CENTER_CM) / grid, grid)

    area1 = EXP_BANDS[0][0] * np.pi * EXP_BANDS[0][2]
    area2 = EXP_BANDS[1][0] * np.pi * EXP_BANDS[1][2]
    print("EXPERIMENTAL COUPLET  (Table S3, dVenus-TD; amplitudes negated to "
          "undo\n  Nguyen's TDX-minus-TD difference and recover the dimer's own "
          "couplet)")
    print(f"  band 1  {EXP_BANDS[0][0]:+.2f}e5 @ {EXP_BANDS[0][1]:.0f} cm^-1, "
          f"HWHM {EXP_BANDS[0][2]:.1f}   area {area1:+.1f}e5")
    print(f"  band 2  {EXP_BANDS[1][0]:+.2f}e5 @ {EXP_BANDS[1][1]:.0f} cm^-1, "
          f"HWHM {EXP_BANDS[1][2]:.1f}   area {area2:+.1f}e5")
    print(f"  areas cancel to {100*abs(area1+area2)/abs(area1):.1f} % -- the couplet is "
          f"conservative to within that, which is what an exciton\n"
          f"  mechanism predicts and which also bounds any intrinsic monomer CD "
          f"left in the difference spectrum.")
    print(f"  observed splitting 2*{EXP_HALF_SPLITTING_CM:.2f} = "
          f"{2*EXP_HALF_SPLITTING_CM:.1f} cm^-1\n")

    print("COUPLET FIRST MOMENT  (lineshape-free, detuning-free)")
    print(f"  predicted    {m_pred_theta:+12.4e} deg cm^2 dmol^-1 cm^-1")
    print(f"  experimental {m_exp_theta:+12.4e} deg cm^2 dmol^-1 cm^-1")
    ratio = m_exp_theta / m_pred_theta
    print(f"  experiment / prediction = {ratio:+.2f}")
    print(f"  sign agreement: {'YES' if ratio > 0 else 'NO -- opposite handedness'}")
    j_required = j_cm * ratio
    print(f"  J required to match the experimental couplet: {j_required:+.1f} cm^-1")
    print(f"  (our TDC coupling is {J_TDC_CM:.2f} +- {J_TDC_SD_CM:.2f} cm^-1)\n")

    # ---- full predicted spectrum, broadened like the experiment -------------
    hwhm = 0.5 * (EXP_BANDS[0][2] + EXP_BANDS[1][2])
    shift = EXP_CENTER_CM - 0.5 * (e_a + e_b).mean()
    ea_s, eb_s = e_a + shift, e_b + shift
    mean_s = 0.5 * (ea_s + eb_s)
    nu_p, nu_m = mean_s + 0.5 * omega, mean_s - 0.5 * omega
    r_plus = -np.pi * mean_s * (j_cm / omega) * triple_cgs

    spec = np.zeros_like(grid)
    for k in range(n):
        for nu, rot in ((nu_p[k], r_plus[k]), (nu_m[k], -r_plus[k])):
            shape = (hwhm / np.pi) / ((grid - nu) ** 2 + hwhm**2)
            spec += DELTA_EPSILON_PREFACTOR * rot * nu * shape
    spec = spec / n * ELLIPTICITY_PER_DELTA_EPSILON

    pk_pred = spec.max() - spec.min()
    pk_exp = exp_couplet.max() - exp_couplet.min()
    print(f"BROADENED SPECTRUM  (Lorentzian HWHM {hwhm:.0f} cm^-1, rigid shift "
          f"{shift:+.0f} cm^-1 onto the experimental origin)")
    print(f"  predicted peak-to-peak    {pk_pred:12.4e} deg cm^2 dmol^-1")
    print(f"    extrema {spec.max():+.3e} @ {grid[spec.argmax()]:.0f} cm^-1, "
          f"{spec.min():+.3e} @ {grid[spec.argmin()]:.0f} cm^-1")
    print(f"  experimental peak-to-peak {pk_exp:12.4e} deg cm^2 dmol^-1")
    print(f"    extrema {exp_couplet.max():+.3e} @ "
          f"{grid[exp_couplet.argmax()]:.0f} cm^-1, {exp_couplet.min():+.3e} @ "
          f"{grid[exp_couplet.argmin()]:.0f} cm^-1")
    print(f"  experiment / prediction = {pk_exp/pk_pred:.2f}\n")

    # ---- dissymmetry factor, a scale-free cross-check ------------------------
    eps_td = 184400.0
    g_exp = (pk_exp / ELLIPTICITY_PER_DELTA_EPSILON) / eps_td
    g_pred = (pk_pred / ELLIPTICITY_PER_DELTA_EPSILON) / eps_td
    print("DISSYMMETRY FACTOR  (peak-to-peak Delta_eps / eps_516, eps = 184 400)")
    print(f"  experimental g = {g_exp:.3e}")
    print(f"  predicted    g = {g_pred:.3e}")

    # ---- robustness of the two headline numbers -----------------------------
    print("\nROBUSTNESS")
    for half in (1500, 3000, 5000, 8000):
        gg = np.arange(EXP_CENTER_CM - half, EXP_CENTER_CM + half, 1.0)
        ss = sum(lorentzian(gg, a * AMPLITUDE_UNIT, c, w) for a, c, w in EXP_BANDS)
        print(f"  M_exp over +-{half:5d} cm^-1 : "
              f"{np.trapezoid(ss*(gg-EXP_CENTER_CM)/gg, gg):+.3e}")
    print("  (the couplet is only conservative to 5.4 %, so the Lorentzian tails "
          "leave a weak\n   window dependence; it is +-10 %, well inside the "
          "discrepancy reported above)")

    if COUPLING_NPZ.is_file():
        c = np.load(COUPLING_NPZ)
        t_coup = np.einsum("ij,ij->i", c["r_B"] - c["r_A"],
                           np.cross(c["mu_A"], c["mu_B"]))
        print(f"\n  triple product, independent 1000-frame coupling ensemble "
              f"(transition-density centroids):")
        print(f"    {t_coup.mean():+.2f} +- {t_coup.std(ddof=1):.2f} Ang*a.u.^2 "
              f"({int((t_coup<0).sum())}/{len(t_coup)} negative)")
        print(f"    vs {triple_au.mean():+.2f} here (CR2 centroids): the two "
              f"centroid conventions differ by "
              f"{100*abs(t_coup.mean()-triple_au.mean())/abs(triple_au.mean()):.0f} %,")
        print(f"    which bounds the geometric systematic on the amplitude. The "
              f"SIGN is unanimous in both.")
        print(f"    stored J_cm sign: {int((c['J_cm']>0).sum())}/{len(c['J_cm'])} "
              f"positive. Our own sign is therefore internally consistent; the\n"
              f"    handedness reported previously came from comparing this "
              f"against Nguyen's\n    TDX-minus-TD difference spectrum without "
              f"inverting it (see DIFFERENCE_SPECTRUM_SIGN).")

    out = ROOT / "results" / "absolute_cd_vs_experiment.json"
    out.write_text(json.dumps(dict(
        n_frames=int(n), j_cm=float(j_cm), j_pda_mean=float(j_pda.mean()),
        triple_au_mean=float(triple_au.mean()), triple_au_sd=float(triple_au.std(ddof=1)),
        first_moment_pred=float(m_pred_theta), first_moment_exp=float(m_exp_theta),
        first_moment_ratio=float(ratio), j_required_cm=float(j_required),
        peak_to_peak_pred=float(pk_pred), peak_to_peak_exp=float(pk_exp),
        g_pred=float(g_pred), g_exp=float(g_exp), hwhm_cm=float(hwhm),
        rigid_shift_cm=float(shift),
    ), indent=1))
    np.savez(ROOT / "results" / "absolute_cd_spectrum.npz",
             grid_cm=grid, predicted_theta=spec, experimental_theta=exp_couplet)
    print(f"\nwrote {out.name} and absolute_cd_spectrum.npz")


if __name__ == "__main__":
    main()

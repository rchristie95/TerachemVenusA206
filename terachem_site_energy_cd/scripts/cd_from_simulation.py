#!/usr/bin/env python3
"""Forward-predict the Venus tandem CD couplet from the simulation alone.

WHY. The manuscript's CD argument is entirely negative: it shows that the
couplet SEPARATION cannot be read as a Hamiltonian gap, and asserts that the
amplitude and first moment are the observables that carry J. It never actually
draws the predicted couplet. Every existing CD script in this directory either
runs OUR spectrum through THEIR fit (nguyen_fit_on_simulation.py), or fits a
Lorentzian triple to a shape and compares widths (cd_width_with_computed_
disorder.py, vibronic_cd_vs_experiment.py), or compares a single scalar moment
(absolute_cd_vs_experiment.py). None of them produces Delta_epsilon(lambda) on
an absolute molar axis next to the measured spectrum.

They also all run on ens_v2_all.npz -- 95 frames of the superseded TeraChem
CAM-B3LYP ensemble, with a single fixed J = 32.82 cm^-1 and rigidly-fitted
transition densities. This runs on the 500-frame v5 campaign, where J, both
transition dipoles and both site energies are recomputed for every frame.

WHAT IS AND IS NOT TAKEN FROM EXPERIMENT
  computed here          per-frame E_A, E_B, J, mu_A, mu_B, R_AB (v5camp)
  computed elsewhere     HWHM 340 cm^-1 (cumulant on the measured gap ACF),
                         sigma_site,slow 159 cm^-1 (1 ns minus the 592 fs
                         window), S = 0.253 and omega = 1369.7 cm^-1
                         (CAM-B3LYP frequencies), shift -0.9745 eV (8 matched
                         STEOM/TDDFT site pairs)
  taken from experiment  NOTHING in the `abinitio` origin case; the band centre
                         only, in the `anchored` case.

THE DOUBLE-COUNTING TRAP. The 340 cm^-1 cumulant width already contains every
site-energy fluctuation faster than 592 fs. The raw v5 frame-to-frame detuning
spread (sd 762 cm^-1) contains those same fluctuations, because each frame is an
instantaneous snapshot. Using both charges the same nuclear motion twice. The
primary case therefore rescales the ensemble detuning down to the measured SLOW
component (sigma_gap = sqrt(2) * 159 = 225 cm^-1), preserving each frame's
correlation between detuning and geometry, and leaves the fast part in the
lineshape. The raw case is reported alongside as the (over-broad) upper bound.

ABSOLUTE SCALE. Two standard cgs relations fix it with no free normalisation:
    R = 2.296e-39 int (Delta_eps / nubar) d nubar
    D = 9.184e-39 int (eps      / nubar) d nubar
so both Delta_eps and eps come out in M^-1 cm^-1 and their ratio is the
dissymmetry g_abs, which is independent of |mu| to leading order and is
therefore the comparison that survives the known transition-dipole systematic.
Note that vibronic_exciton.vibronic_spectra omits the nubar_k weight and so is
NOT on this absolute axis; the per-state assembly below applies it.

Read-only apart from results/cd_from_simulation.{json,npz,pdf}.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from absolute_cd import ELLIPTICITY_PER_DELTA_EPSILON            # noqa: E402
from vibronic_exciton import frame_vibronic_observables          # noqa: E402

ENSEMBLE = ROOT / "results" / "ensembles" / "ens_v5camp.npz"
OUT = ROOT / "results"

# --- computed inputs, all from elsewhere in the repo -------------------------
HWHM_CM = 340.0            # cumulant on site A's 592 fs QM/MM gap ACF
SIGMA_SITE_SLOW = 159.0    # what 1 ns adds beyond the 592 fs window
HUANG_RHYS = 0.253         # CAM-B3LYP/6-31G* excited-state frequencies
OMEGA_VIB = 1369.7
STEOM_SHIFT_EV = -0.9745050587499999      # 8 matched STEOM/TDDFT site pairs
CM_PER_EV = 8065.543937

# --- experiment, used only for the overlay and the anchored origin ----------
EXP_CENTER_CM = 19322.0
EXP_HALF_SPLIT = 130.79
DIFF_SIGN = -1.0           # Delta[theta] = TDX - TD  ->  TD couplet
EXP_BANDS = tuple((DIFF_SIGN * a * 1e5, c, w) for a, c, w in (
    (+3.59, EXP_CENTER_CM - EXP_HALF_SPLIT, 293.51),
    (-2.27, EXP_CENTER_CM + EXP_HALF_SPLIT, 439.25)))
EXP_THIRD = (DIFF_SIGN * -0.34e5, 20774.0, 565.79)
EPS_TD_EXP = 184400.0      # tandem-dimer molar extinction, 2 x Venus monomer

# --- cgs conventions --------------------------------------------------------
R_TO_DELTA_EPS = 1.0 / 2.296e-39
D_TO_EPS = 1.0 / 9.184e-39
DEBYE_ESU_CM = 1.0e-18

GRID = np.arange(16500.0, 23500.0, 2.0)
N_MAX = 5


def lorentzian(grid, amplitude, centre, hwhm):
    return amplitude / (1.0 + ((grid - centre) / hwhm) ** 2)


def unit_lorentzian(grid, centre, hwhm):
    """Normalised to unit area over wavenumber, so the cgs prefactors hold."""
    return hwhm / np.pi / ((grid - centre) ** 2 + hwhm ** 2)


def build(e_a, e_b, j, triple, cos_alpha, mu2_cgs, hwhm):
    """Absolute Delta_eps and eps (M^-1 cm^-1) for one ensemble."""
    d_eps = np.zeros_like(GRID)
    eps = np.zeros_like(GRID)
    for ea, eb, jj, tt, ca, m2 in zip(e_a, e_b, j, triple, cos_alpha, mu2_cgs):
        obs = frame_vibronic_observables(ea, eb, jj, tt, OMEGA_VIB,
                                         HUANG_RHYS, N_MAX, cos_alpha=ca)
        for nu, dstr, rot in zip(obs["energies_cm"], obs["dipole_strengths"],
                                 obs["rotational_strengths_relative"]):
            shape = unit_lorentzian(GRID, nu, hwhm)
            d_eps += R_TO_DELTA_EPS * rot * nu * shape
            eps += D_TO_EPS * dstr * m2 * nu * shape
    n = len(e_a)
    return d_eps / n, eps / n


def fwhm(grid, spectrum):
    half = spectrum.max() / 2.0
    idx = np.where(spectrum >= half)[0]
    return grid[idx[-1]] - grid[idx[0]] if idx.size else float("nan")


def couplet_metrics(grid, spectrum):
    i_hi, i_lo = int(np.argmax(spectrum)), int(np.argmin(spectrum))
    return dict(
        peak_to_peak=float(spectrum.max() - spectrum.min()),
        nu_max=float(grid[i_hi]), nu_min=float(grid[i_lo]),
        separation_cm=abs(float(grid[i_hi] - grid[i_lo])),
        lobe_ratio=float(abs(spectrum.max() / spectrum.min()))
        if spectrum.min() else float("nan"),
        first_moment=float(np.trapezoid(
            spectrum * (grid - EXP_CENTER_CM) / grid, grid)),
        area=float(np.trapezoid(spectrum, grid)),
    )


def main() -> None:
    d = np.load(ENSEMBLE)
    e_a_raw, e_b_raw = d["e_a"], d["e_b"]
    j, triple, cos_alpha = d["j"], d["triple"], d["cos_alpha"]
    mu2_cgs = (d["mu_a"] * d["mu_b"]) * DEBYE_ESU_CM ** 2     # |mu|^2, esu^2cm^2

    delta_raw = e_a_raw - e_b_raw
    mean_raw = 0.5 * (e_a_raw + e_b_raw)
    sd_delta_raw = float(delta_raw.std(ddof=1))
    sd_mean_raw = float(mean_raw.std(ddof=1))

    sigma_gap_slow = np.sqrt(2.0) * SIGMA_SITE_SLOW
    sigma_mean_slow = SIGMA_SITE_SLOW / np.sqrt(2.0)
    k_delta = sigma_gap_slow / sd_delta_raw
    k_mean = sigma_mean_slow / sd_mean_raw

    print("=== v5camp CD forward model ===")
    print(f"  frames                {len(j)}")
    print(f"  <J>                   {j.mean():.2f} +/- {j.std(ddof=1):.2f} cm^-1")
    print(f"  <|Delta|> raw         {np.abs(delta_raw).mean():.1f} cm^-1 "
          f"(sd {sd_delta_raw:.1f})")
    print(f"  <T>                   {triple.mean():+.3e} cgs, "
          f"one sign in {100*max(np.mean(triple>0), np.mean(triple<0)):.0f}% of frames")
    print(f"  |mu| A/B              {d['mu_a'].mean():.2f} / {d['mu_b'].mean():.2f} D")
    print(f"  cos(alpha)            {cos_alpha.mean():+.4f}  "
          f"(alpha = {np.degrees(np.arccos(cos_alpha.mean())):.1f} deg)")
    print(f"  fast/slow rescale     detuning x {k_delta:.4f}, mean x {k_mean:.4f}")

    # ---- energy origin -----------------------------------------------------
    shift_abinitio = STEOM_SHIFT_EV * CM_PER_EV
    centre_abinitio = mean_raw.mean() + shift_abinitio
    shift_anchor = EXP_CENTER_CM - mean_raw.mean()
    print(f"  TDDFT band centre     {mean_raw.mean():.0f} cm^-1 "
          f"({1e7/mean_raw.mean():.1f} nm)")
    print(f"  + STEOM shift         {centre_abinitio:.0f} cm^-1 "
          f"({1e7/centre_abinitio:.1f} nm)   [ab initio]")
    print(f"  experimental centre   {EXP_CENTER_CM:.0f} cm^-1 "
          f"({1e7/EXP_CENTER_CM:.1f} nm)")
    print(f"  residual origin error {centre_abinitio - EXP_CENTER_CM:+.0f} cm^-1")

    cases = {}
    for name, kd, km in (("primary_slow", k_delta, k_mean),
                         ("raw_detuning", 1.0, 1.0),
                         ("degenerate", 0.0, 0.0)):
        mean_c = mean_raw.mean() + km * (mean_raw - mean_raw.mean()) + shift_anchor
        half = 0.5 * kd * delta_raw
        cases[name] = (mean_c + half, mean_c - half)

    results, spectra = {}, {}
    for name, (ea, eb) in cases.items():
        d_eps, eps = build(ea, eb, j, triple, cos_alpha, mu2_cgs, HWHM_CM)
        theta = ELLIPTICITY_PER_DELTA_EPSILON * d_eps
        m = couplet_metrics(GRID, theta)
        om = np.hypot(ea - eb, 2 * j)
        m.update(
            eps_max=float(eps.max()),
            eps_max_ratio_to_experiment=float(eps.max() / EPS_TD_EXP),
            g_abs_peak=float(np.abs(d_eps).max() / eps.max()),
            eps_fwhm_cm=float(fwhm(GRID, eps)),
            dipole_strength_integral=float(np.trapezoid(eps / GRID, GRID)),
            mean_suppression=float(np.mean(2 * j / om)),
            median_suppression=float(np.median(2 * j / om)),
            mean_abs_delta=float(np.abs(ea - eb).mean()),
        )
        results[name] = m
        spectra[name] = (theta, eps)

    exp_theta = sum(lorentzian(GRID, *b) for b in EXP_BANDS)
    exp_full = exp_theta + lorentzian(GRID, *EXP_THIRD)
    m_exp = couplet_metrics(GRID, exp_theta)
    m_exp["eps_max"] = EPS_TD_EXP
    m_exp["g_abs_peak"] = float(
        np.abs(exp_theta).max() / ELLIPTICITY_PER_DELTA_EPSILON / EPS_TD_EXP)
    results["experiment"] = m_exp

    print("\n  case            pk-pk theta   sep    lobe    1st moment   "
          "eps_max     g_abs")
    for name in ("primary_slow", "raw_detuning", "degenerate", "experiment"):
        m = results[name]
        print(f"  {name:<14} {m['peak_to_peak']:10.3e} {m['separation_cm']:6.0f} "
              f"{m['lobe_ratio']:6.2f} {m['first_moment']:12.3e} "
              f"{m['eps_max']:10.3e} {m['g_abs_peak']:9.2e}")

    p = results["primary_slow"]
    print(f"\n  predicted / observed   peak-to-peak "
          f"{p['peak_to_peak']/m_exp['peak_to_peak']:.2f}x, "
          f"first moment {p['first_moment']/m_exp['first_moment']:.2f}x, "
          f"g_abs {p['g_abs_peak']/m_exp['g_abs_peak']:.2f}x")
    print(f"  predicted eps_max is {p['eps_max_ratio_to_experiment']:.2f}x the "
          f"measured {EPS_TD_EXP:.0f} M^-1 cm^-1, but the predicted absorption "
          f"band is only {p['eps_fwhm_cm']:.0f} cm^-1 FWHM")
    print("  eps_max is therefore width-confounded; the |mu|-only quantity is "
          f"int(eps/nu)dnu = {p['dipole_strength_integral']:.4e}, which needs the "
          "MEASURED absorption band, not its peak, to complete.")

    payload = dict(
        inputs=dict(n_frames=int(len(j)), hwhm_cm=HWHM_CM,
                    sigma_site_slow_cm=SIGMA_SITE_SLOW, huang_rhys=HUANG_RHYS,
                    omega_vib_cm=OMEGA_VIB, steom_shift_eV=STEOM_SHIFT_EV,
                    j_mean_cm=float(j.mean()), j_sd_cm=float(j.std(ddof=1)),
                    triple_mean_cgs=float(triple.mean()),
                    cos_alpha_mean=float(cos_alpha.mean()),
                    mu_a_D=float(d["mu_a"].mean()), mu_b_D=float(d["mu_b"].mean()),
                    sd_delta_raw_cm=sd_delta_raw, sd_mean_raw_cm=sd_mean_raw,
                    rescale_delta=float(k_delta), rescale_mean=float(k_mean)),
        origin=dict(tddft_centre_cm=float(mean_raw.mean()),
                    abinitio_centre_cm=float(centre_abinitio),
                    experimental_centre_cm=EXP_CENTER_CM,
                    residual_cm=float(centre_abinitio - EXP_CENTER_CM)),
        cases=results,
    )
    (OUT / "cd_from_simulation.json").write_text(json.dumps(payload, indent=2))
    np.savez(OUT / "cd_from_simulation.npz", grid_cm=GRID,
             experiment_theta=exp_theta, experiment_theta_with_third=exp_full,
             **{f"{k}_theta": v[0] for k, v in spectra.items()},
             **{f"{k}_eps": v[1] for k, v in spectra.items()})
    print(f"\n  wrote {OUT/'cd_from_simulation.json'}")


if __name__ == "__main__":
    main()

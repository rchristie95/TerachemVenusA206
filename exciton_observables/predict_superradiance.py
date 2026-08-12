#!/usr/bin/env python3
"""Predict the superradiant lifetime shortening from the computed ensemble.

Nguyen et al. measure dVenus-TD decaying 57 +/- 4 ps faster than the
single-chromophore dVenus-TDX control on a 3.026 ns lifetime. For the two-site
Hamiltonian the eigenstate dipole strengths are

    |mu_pm|^2 = mu^2 (1 +/- sin(2 theta) cos alpha),   sin(2 theta) = 2J/Omega,
    Omega = sqrt(Delta^2 + 4 J^2)

so the *thermally averaged* radiative rate of the emitting manifold is

    k_rad / k_mono = 1 - tanh(Omega / 2 kT) * (2J/Omega) * cos(alpha)

and the measured lifetime shift follows with the quantum yield Phi:

    Delta_tau / tau = Phi * [ k_rad/k_mono - 1 ]
                    = Phi * tanh(Omega/2kT) * (2J/Omega) * (-cos alpha).

Two things the manuscript's Eq. (superradiance) does not do, both of which
matter here:

  1. It omits the thermal factor tanh(Omega/2kT). At Omega ~ 300-600 cm^-1
     against kT = 208.5 cm^-1 that factor is 0.6-0.9, not 1: the upper
     exciton state is measurably populated at 300 K and it is subradiant,
     which cancels part of the lower state's enhancement.
  2. It evaluates the expression at the ENSEMBLE MEAN detuning. sin(2theta)
     is strongly convex in Delta, and the computed detuning distribution is
     broad (SD ~700 cm^-1 with substantial weight near zero), so the mean of
     the observable is not the observable at the mean. The manuscript already
     noticed the symptom -- harmonic-mean Omega 295.9 vs arithmetic 611.2.

This script therefore averages the observable frame by frame over the joint
QM/MM ensemble, where Delta, cos(alpha) and the point-dipole coupling all come
from the SAME frame and the same transition densities, so their relative signs
are physical.

Critically, the sign matters: J > 0 with cos(alpha) < 0 (our obtuse geometry)
gives a BRIGHTER emitting state and a FASTER lifetime, matching the sign of
the measurement. The anisotropy-derived cos(alpha) = +0.660 gives the opposite
sign -- subradiance, a SLOWER lifetime -- contradicting the same experiment it
is derived from.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.63
KT_CM_300K = 208.509  # k_B T at 300 K in cm^-1

# Experiment (Nguyen et al. 2025, dVenus)
TAU_TDX_NS = 3.026
DTAU_PS = 57.0
DTAU_SD_PS = 4.0
PHI_VENUS = 0.57
COS_ALPHA_EXPERIMENTAL = 0.660  # from their complete-transfer anisotropy analysis

# Ensemble coupling calibration: full transition-density coupling exceeds the
# point-dipole value by this factor of the ensemble means (production README).
TDC_OVER_PDA = 1.1872


def point_dipole_coupling_cm(mu_a, mu_b, r_a_ang, r_b_ang, epsilon):
    """Signed point-dipole coupling, distances converted to bohr correctly.

    A Coulomb sum whose distances are in Angstrom converts with
    BOHR_TO_ANGSTROM, not ANGSTROM_TO_BOHR. Here the separation is converted
    to bohr up front instead, which is the same statement.
    """
    r_vec = (r_b_ang - r_a_ang) * ANGSTROM_TO_BOHR
    r = np.linalg.norm(r_vec, axis=1)
    r_hat = r_vec / r[:, None]
    jdd = np.sum(mu_a * mu_b, axis=1) - 3.0 * np.sum(mu_a * r_hat, axis=1) * np.sum(mu_b * r_hat, axis=1)
    return jdd / (r ** 3 * epsilon) * HARTREE_TO_CM


def observable(j_cm, delta_cm, cos_alpha, kt_cm=KT_CM_300K, thermal=True):
    """Fractional radiative-rate enhancement, per frame."""
    omega = np.sqrt(delta_cm ** 2 + 4.0 * j_cm ** 2)
    mixing = 2.0 * j_cm / omega
    thermal_factor = np.tanh(omega / (2.0 * kt_cm)) if thermal else np.ones_like(omega)
    return -thermal_factor * mixing * cos_alpha, omega, mixing, thermal_factor


def block_sem(x, n_blocks=5):
    blocks = np.array_split(np.asarray(x), n_blocks)
    return float(np.std([b.mean() for b in blocks], ddof=1) / np.sqrt(n_blocks))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ensembles", nargs="+", type=Path, default=[
        Path("terachem_site_energy_cd/results/ensembles/ens_v2_all.npz"),
        Path("terachem_site_energy_cd/results/ensembles/ens_v2.npz"),
        Path("terachem_site_energy_cd/results/ensembles/ens_final.npz"),
    ])
    ap.add_argument("--epsilon", type=float, default=1.77)
    ap.add_argument("--out", type=Path, default=Path("exciton_observables/superradiance_prediction.json"))
    args = ap.parse_args()

    target = DTAU_PS / (TAU_TDX_NS * 1000.0)
    target_sd = DTAU_SD_PS / (TAU_TDX_NS * 1000.0)
    required = target / PHI_VENUS
    print(f"Experiment: Delta_tau/tau = {target:.5f} +/- {target_sd:.5f} "
          f"({DTAU_PS} +/- {DTAU_SD_PS} ps on {TAU_TDX_NS} ns)")
    print(f"  => required |(2J/Omega) cos alpha| = (1/Phi)(Dtau/tau) = {required:.4f}\n")

    results = {"experiment": {"dtau_ps": DTAU_PS, "dtau_sd_ps": DTAU_SD_PS,
                              "tau_ns": TAU_TDX_NS, "phi": PHI_VENUS,
                              "dtau_over_tau": target,
                              "required_mixing_times_cos": required},
               "ensembles": {}}

    for path in args.ensembles:
        if not path.exists():
            print(f"[skip] {path} not found")
            continue
        d = np.load(path)
        delta = d["e_a_cm"] - d["e_b_cm"]
        mu_a, mu_b = d["mu_a_au"], d["mu_b_au"]
        cos_alpha = np.sum(mu_a * mu_b, axis=1) / (
            np.linalg.norm(mu_a, axis=1) * np.linalg.norm(mu_b, axis=1))
        j_pda = point_dipole_coupling_cm(mu_a, mu_b, d["r_a_ang"], d["r_b_ang"], args.epsilon)
        j_tdc = j_pda * TDC_OVER_PDA

        enh, omega, mixing, thermal = observable(j_tdc, delta, cos_alpha)
        dtau_pred_ps = PHI_VENUS * enh * TAU_TDX_NS * 1000.0

        # Mean-field ("manuscript style") comparison: evaluate at mean |Delta|,
        # mean cos alpha, mean J, and without the thermal factor.
        j_bar, d_bar, c_bar = j_tdc.mean(), np.abs(delta).mean(), cos_alpha.mean()
        mf_no_thermal, _, mf_mix, _ = observable(j_bar, d_bar, c_bar, thermal=False)
        mf_thermal, _, _, _ = observable(j_bar, d_bar, c_bar, thermal=True)

        n = len(delta)
        name = path.stem
        print(f"=== {name}  (n = {n}) ===")
        print(f"  <|Delta|>            {np.abs(delta).mean():8.1f} cm^-1   "
              f"median {np.median(np.abs(delta)):7.1f}   SD {delta.std(ddof=1):7.1f}")
        print(f"  <cos alpha>          {cos_alpha.mean():+8.4f}          "
              f"(|cos| {np.abs(cos_alpha).mean():.4f}, SD {cos_alpha.std(ddof=1):.4f})")
        print(f"  <J_TDC>              {j_tdc.mean():8.2f} cm^-1   "
              f"(PDA {j_pda.mean():.2f}, SD {j_tdc.std(ddof=1):.2f})")
        print(f"  <Omega>              {omega.mean():8.1f} cm^-1   "
              f"harmonic {n/np.sum(1.0/omega):7.1f}")
        print(f"  <thermal factor>     {thermal.mean():8.4f}")
        print(f"  ensemble <enhancement> {enh.mean():+.5f} +/- {block_sem(enh):.5f}")
        print(f"  PREDICTED Delta_tau  {dtau_pred_ps.mean():7.2f} +/- "
              f"{PHI_VENUS*block_sem(enh)*TAU_TDX_NS*1000.0:.2f} ps"
              f"      [measured {DTAU_PS} +/- {DTAU_SD_PS} ps]")
        print(f"    mean-field, no thermal factor (manuscript form): "
              f"{PHI_VENUS*float(mf_no_thermal)*TAU_TDX_NS*1000.0:7.2f} ps")
        print(f"    mean-field, with thermal factor:                 "
              f"{PHI_VENUS*float(mf_thermal)*TAU_TDX_NS*1000.0:7.2f} ps")

        # Sign test against the anisotropy-derived angle, same J and Delta.
        enh_exp_angle, _, _, _ = observable(j_tdc, delta, np.full_like(cos_alpha,
                                                                       COS_ALPHA_EXPERIMENTAL))
        dtau_exp_angle = PHI_VENUS * enh_exp_angle.mean() * TAU_TDX_NS * 1000.0
        print(f"    with their cos alpha = +{COS_ALPHA_EXPERIMENTAL}: "
              f"{dtau_exp_angle:+7.2f} ps  "
              f"({'SLOWER - wrong sign' if dtau_exp_angle < 0 else 'faster'})")

        # Inverse problem: what detuning would our angle+J need to hit 57 ps?
        # Solve <tanh(O/2kT) (2J/O) (-cos a)> = required, scaling all |Delta|.
        scales = np.geomspace(0.05, 5.0, 4000)
        preds = np.array([observable(j_tdc, delta * s, cos_alpha)[0].mean() for s in scales])
        idx = int(np.argmin(np.abs(preds - target / PHI_VENUS)))
        implied_mean_abs_delta = float(np.abs(delta * scales[idx]).mean()) if preds.min() <= target/PHI_VENUS <= preds.max() else float("nan")
        print(f"    detuning scale reproducing the measurement: x{scales[idx]:.3f} "
              f"=> <|Delta|> = {implied_mean_abs_delta:.1f} cm^-1")

        results["ensembles"][name] = {
            "n": int(n),
            "mean_abs_delta_cm": float(np.abs(delta).mean()),
            "median_abs_delta_cm": float(np.median(np.abs(delta))),
            "mean_cos_alpha": float(cos_alpha.mean()),
            "mean_J_tdc_cm": float(j_tdc.mean()),
            "mean_J_pda_cm": float(j_pda.mean()),
            "mean_omega_cm": float(omega.mean()),
            "harmonic_omega_cm": float(n / np.sum(1.0 / omega)),
            "mean_thermal_factor": float(thermal.mean()),
            "predicted_enhancement": float(enh.mean()),
            "predicted_enhancement_sem": block_sem(enh),
            "predicted_dtau_ps": float(dtau_pred_ps.mean()),
            "predicted_dtau_sem_ps": float(PHI_VENUS * block_sem(enh) * TAU_TDX_NS * 1000.0),
            "meanfield_no_thermal_dtau_ps": float(PHI_VENUS * mf_no_thermal * TAU_TDX_NS * 1000.0),
            "meanfield_thermal_dtau_ps": float(PHI_VENUS * mf_thermal * TAU_TDX_NS * 1000.0),
            "dtau_with_experimental_angle_ps": float(dtau_exp_angle),
            "implied_mean_abs_delta_cm": implied_mean_abs_delta,
        }
        print()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

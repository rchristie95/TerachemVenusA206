#!/usr/bin/env python3
"""Does the CD couplet width close once the disorder is split fast/slow?

WHY. The standing verdict (audit: cd-forward-model-verdict) is that our
site-energy disorder is ~5x too large for the observed couplet width, which
threatens the detuning-dominated mechanism. That verdict treated
<|Delta|> = 576 cm^-1 as STATIC. It is not: the 4 fs window gives a site-energy
correlation time of 8.5 fs, and averaging over hbar/2J = 81 fs leaves only
~170-290 cm^-1 of the 494 standing.

So the old comparison double-counted -- charging the same motion once as
homogeneous dephasing and again as static detuning, exactly the hazard
cumulant_lineshape.py was written to avoid.

THE SPLIT, chosen so it cannot double-count. The cumulant is exact for Gaussian
fluctuations at ANY timescale, so everything actually measured in the 592 fs
window belongs in g(t):

    faster than 592 fs  -> cumulant lineshape from the measured C(t),
                           site A, clean:  HWHM 340 cm^-1
    slower than 592 fs  -> static detuning, from what the 1 ns ensemble adds
                           beyond the 592 fs window:
                           sd 519 (1 ns) vs 494 (592 fs) -> sigma_site = 159,
                           sigma_gap = sqrt(2)*159 = 225 cm^-1

Both numbers are computed, not assumed. Vibronic parameters likewise
(S = 0.253, omega = 1369.7 cm^-1, TeraChem CAM-B3LYP/6-31G*).

CAVEATS carried into the output: the 592 fs and 1 ns numbers come from
DIFFERENT trajectories, so sigma_site = 159 is approximate; the cumulant
lineshape at sigma*tau_c = 0.79 is intermediate Gaussian/Lorentzian while the
forward model and Nguyen's fit both use Lorentzians. A sensitivity scan over
the full range of slow-component estimates (159-287 per site) is therefore run
alongside the primary.

Read-only apart from results/cd_width_with_computed_disorder.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from vibronic_exciton import vibronic_spectra                      # noqa: E402

ENSEMBLE = ROOT / "results" / "ensembles" / "ens_v2_all.npz"
CENTRE = 19322.0
J_CM = 32.8165
GRID = np.arange(17500.0, 23500.0, 2.0)
FITW = (18000.0, 21800.0)

# computed, this session -- see tddft_ensemble_decoherence/audit/
HWHM_COMPUTED = 340.0            # cumulant, site A full 592 fs C(t)
SIGMA_SITE_SLOW = 159.0          # what 1 ns adds beyond the 592 fs window
S_COMPUTED = 0.253
OMEGA_COMPUTED = 1369.7

# Nguyen 2025 Table S3
EXP = dict(w1=293.51, w2=439.25, delta=130.79,
           asym=1.982e5 / 0.8466e5, third_offset=20774.0 - CENTRE)


def model(nu, a1, w1, a2, w2, delta, a3, nu3, w3):
    return (a1 / (1 + ((nu - (CENTRE - delta)) / w1) ** 2)
            + a2 / (1 + ((nu - (CENTRE + delta)) / w2) ** 2)
            + a3 / (1 + ((nu - nu3) / w3) ** 2))


def robust_fit(nu, y):
    """Nguyen's constrained 3-Lorentzian, scale-normalised and multi-start."""
    y = y / np.abs(y).max()
    best = None
    for d0 in (60., 130., 200., 275., 350., 450., 600.):
        for s1 in (1., -1.):
            try:
                p, _ = curve_fit(
                    model, nu, y,
                    p0=[s1, 293.51, -0.6 * s1, 439.25, d0, -0.1 * s1, 20774., 565.79],
                    bounds=([-10, 50, -10, 50, 5, -10, 20000, 50],
                            [10, 2000, 10, 2000, 2000, 10, 21600, 2000]),
                    maxfev=400000)
            except Exception:
                continue
            residual = y - model(nu, *p)
            r2 = 1 - np.sum(residual ** 2) / np.sum((y - y.mean()) ** 2)
            if best is None or r2 > best[-1]:
                best = (p[4], p[1], p[3], p[5], p[6], r2)
    return best


def build_frames(data: dict, sigma_gap_cm: float) -> list:
    """Ensemble with the detuning rescaled to a target static spread.

    Geometry (triple product, separation) and the SHAPE of the detuning
    distribution are kept; only its scale is set, because the 1 ns ensemble's
    sd(Delta) = 686 is dominated by motion that is NOT static on the exciton
    timescale.
    """
    e_a, e_b = data["e_a_cm"], data["e_b_cm"]
    triple = np.einsum("ij,ij->i", data["r_b_ang"] - data["r_a_ang"],
                       np.cross(data["mu_a_au"], data["mu_b_au"]))
    gap = e_a - e_b
    scaled = gap * (sigma_gap_cm / gap.std())
    return list(zip(CENTRE + 0.5 * scaled, CENTRE - 0.5 * scaled,
                    np.full(len(gap), J_CM), triple)), scaled


def run(frames, hwhm, S, omega):
    _, cd = vibronic_spectra(GRID, frames, hwhm, omega_cm=omega,
                             huang_rhys=S, n_max=5)
    mask = (GRID >= FITW[0]) & (GRID <= FITW[1])
    y = cd[mask]
    asym = abs(y.min()) / y.max() if y.max() > 0 else np.nan
    delta, w1, w2, a3, nu3, r2 = robust_fit(GRID[mask], y)
    return dict(delta=delta, w1=w1, w2=w2, asym=float(asym),
                third_offset=nu3 - CENTRE, r2=r2)


def main() -> None:
    data = dict(np.load(ENSEMBLE))
    print(f"target (Nguyen Table S3): widths {EXP['w1']:.1f} / {EXP['w2']:.1f} cm^-1"
          f"  (mean {0.5*(EXP['w1']+EXP['w2']):.1f}), lobe asym {EXP['asym']:.2f}\n")

    cases = [
        ("OLD  static 576, HWHM 366", 686.5, 366.4, 0.35, 1450.0),
        ("old disorder, computed HWHM", 686.5, HWHM_COMPUTED, S_COMPUTED, OMEGA_COMPUTED),
        ("PRIMARY  computed both", np.sqrt(2) * SIGMA_SITE_SLOW, HWHM_COMPUTED,
         S_COMPUTED, OMEGA_COMPUTED),
        ("zero static disorder", 1.0, HWHM_COMPUTED, S_COMPUTED, OMEGA_COMPUTED),
    ]
    rows = []
    print(f"{'case':<30} {'sd(D)':>7} {'<|D|>':>7} {'HWHM':>6} | "
          f"{'w1':>7} {'w2':>7} {'mean w':>7} {'asym':>6} {'delta':>7} {'R2':>7}")
    print("-" * 108)
    for label, sigma_gap, hwhm, S, omega in cases:
        frames, scaled = build_frames(data, sigma_gap)
        out = run(frames, hwhm, S, omega)
        out.update(case=label, sigma_gap=sigma_gap, hwhm=hwhm, S=S, omega=omega,
                   mean_abs_delta=float(np.abs(scaled).mean()))
        rows.append(out)
        print(f"{label:<30} {sigma_gap:>7.1f} {np.abs(scaled).mean():>7.1f} {hwhm:>6.0f} | "
              f"{out['w1']:>7.1f} {out['w2']:>7.1f} {0.5*(out['w1']+out['w2']):>7.1f} "
              f"{out['asym']:>6.2f} {out['delta']:>7.1f} {out['r2']:>7.4f}")

    print(f"\nsensitivity: slow site component 159-287 cm^-1, HWHM 250-400")
    print(f"{'sigma_site':>10} {'sd(D)':>7} {'HWHM':>6} | {'w1':>7} {'w2':>7} {'mean w':>7} {'asym':>6}")
    print("-" * 62)
    scan = []
    for sigma_site in (159.0, 200.0, 250.0, 287.0):
        for hwhm in (250.0, 300.0, 340.0, 400.0):
            frames, scaled = build_frames(data, np.sqrt(2) * sigma_site)
            out = run(frames, hwhm, S_COMPUTED, OMEGA_COMPUTED)
            out.update(sigma_site=sigma_site, hwhm=hwhm)
            scan.append(out)
            print(f"{sigma_site:>10.0f} {np.sqrt(2)*sigma_site:>7.1f} {hwhm:>6.0f} | "
                  f"{out['w1']:>7.1f} {out['w2']:>7.1f} {0.5*(out['w1']+out['w2']):>7.1f} {out['asym']:>6.2f}")

    destination = ROOT / "results" / "cd_width_with_computed_disorder.json"
    destination.write_text(json.dumps(
        {"experiment": EXP, "cases": rows, "sensitivity": scan,
         "computed_inputs": {"hwhm_cm": HWHM_COMPUTED,
                             "sigma_site_slow_cm": SIGMA_SITE_SLOW,
                             "huang_rhys": S_COMPUTED,
                             "omega_vib_cm": OMEGA_COMPUTED}}, indent=2) + "\n")
    print(f"\nwrote {destination.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

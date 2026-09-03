#!/usr/bin/env python3
"""Does the STEOM calibration help the predicted CD, or hurt it?

Three places where the eight matched STEOM/TDDFT site pairs could be used in
place of the raw TDDFT ensemble, tested separately:

  origin    the constant shift  -0.9745 eV, instead of anchoring the band to the
            experimental 19322 cm^-1;
  dipole    the STEOM ensemble mean |mu| = 8.92 D in place of the per-frame
            TDDFT 9.39 D (applied as a uniform scale, since STEOM dipoles exist
            for only 8 of the 500 frames);
  detuning  the STEOM-on-TDDFT regression slope 1.523, which would rescale every
            site-energy deviation.

Each is reported against the measured band: absorption FWHM 963 cm^-1, width at
0.90 of peak 252 cm^-1, eps_max 1.844e5 M^-1 cm^-1, couplet separation 500 cm^-1,
couplet peak-to-peak 2.829e5 deg cm^2 dmol^-1 (Nguyen Table S3, amplitudes in
units of 1e5 -- confirmed against the Fig. 4b axis label, which reads x10^-5,
and against the plotted extrema, which equal the summed Lorentzians).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402

MODES = [(667.1, 0.0818), (1219.9, 0.1262), (1731.2, 0.0699)]
SIGMA_SITE = 391.5
STEOM_SHIFT_CM = -0.9745050587499999 * 8065.543937
STEOM_SLOPE = 1.523
MU_TDDFT_D = 9.39
MU_STEOM_D = 8.92

MEAS = dict(fwhm=963.0, width90=252.0, eps_max=1.844e5,
            separation=500.0, peak_to_peak=2.829e5, lam_max=516.82)


def run(sigma_site, mu_scale, centre, j, triple, mu_a, mu_b):
    vib_e, fc = R.vibrational_basis(MODES, 3)
    nodes, wts = np.polynomial.hermite_e.hermegauss(16)
    wts = wts / wts.sum()
    base = np.full(j.size, centre)
    de = np.zeros_like(R.GRID)
    ep = np.zeros_like(R.GRID)
    sg = np.sqrt(2.0) * sigma_site
    ma, mb = mu_a * mu_scale, mu_b * mu_scale
    tr = triple * mu_scale ** 2                    # T is bilinear in the dipoles
    for i in range(j.size):
        for nd, wt in zip(nodes, wts):
            dl = sg * nd
            nu, rk, dk = R.frame_sticks(base[i] + 0.5 * dl, base[i] - 0.5 * dl,
                                        j[i], ma[i], mb[i], tr[i], vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * rk * nu * wt, R.GRID, de)
            R.accumulate(nu, R.D_TO_EPS * dk * nu * wt, R.GRID, ep)
    de /= j.size
    ep /= j.size
    sm = sigma_site / np.sqrt(2.0)
    de = R.gaussian_convolve(de, R.GRID, sm)
    ep = R.gaussian_convolve(ep, R.GRID, sm)
    th = R.ELLIPTICITY_PER_DELTA_EPSILON * de
    w = lambda f: (lambda k: float(R.GRID[k[-1]] - R.GRID[k[0]]))(
        np.where(ep >= f * ep.max())[0])
    return dict(lam_max=float(1e7 / R.GRID[int(ep.argmax())]),
                fwhm=w(0.50), width90=w(0.90), eps_max=float(ep.max()),
                separation=abs(float(R.GRID[int(th.argmax())]
                                     - R.GRID[int(th.argmin())])),
                peak_to_peak=float(th.max() - th.min()),
                g_abs=float(np.abs(de).max() / ep.max()))


def main() -> None:
    d = np.load(R.ENSEMBLE)
    j, triple = d["j"], d["triple"]
    mu_a, mu_b = d["mu_a_vec"], d["mu_b_vec"]
    tddft_centre = 0.5 * (d["e_a"] + d["e_b"]).mean()
    abinitio = tddft_centre + STEOM_SHIFT_CM

    variants = [
        ("baseline (anchored, TDDFT mu)", SIGMA_SITE, 1.0, R.EXP_CENTER_CM),
        ("STEOM origin", SIGMA_SITE, 1.0, abinitio),
        ("STEOM dipole 8.92 D", SIGMA_SITE, MU_STEOM_D / MU_TDDFT_D, R.EXP_CENTER_CM),
        ("STEOM detuning slope 1.523", SIGMA_SITE * STEOM_SLOPE, 1.0, R.EXP_CENTER_CM),
        ("all three STEOM corrections", SIGMA_SITE * STEOM_SLOPE,
         MU_STEOM_D / MU_TDDFT_D, abinitio),
    ]
    print("  variant                        lam_max   FWHM  w@0.90    eps_max"
          "   sep   pk-pk theta")
    out = {}
    for tag, sig, sc, ctr in variants:
        r = run(sig, sc, ctr, j, triple, mu_a, mu_b)
        out[tag] = r
        print(f"  {tag:<29} {r['lam_max']:7.1f} {r['fwhm']:6.0f} {r['width90']:7.0f} "
              f"{r['eps_max']:10.3e} {r['separation']:5.0f}  {r['peak_to_peak']:11.3e}")
    print(f"  {'MEASURED':<29} {MEAS['lam_max']:7.2f} {MEAS['fwhm']:6.0f} "
          f"{MEAS['width90']:7.0f} {MEAS['eps_max']:10.3e} {MEAS['separation']:5.0f} "
          f" {MEAS['peak_to_peak']:11.3e}")

    b = out["baseline (anchored, TDDFT mu)"]
    print(f"\n  |mu| implied by the absolute absorption intensity:")
    print(f"    eps_max ratio {b['eps_max']/MEAS['eps_max']:.2f} at matched band width"
          f"  ->  |mu| = {MU_TDDFT_D/np.sqrt(b['eps_max']/MEAS['eps_max']):.2f} D "
          f"against the computed {MU_TDDFT_D:.2f} D")
    print(f"    the CD couplet then carries |mu|^4: "
          f"{(b['eps_max']/MEAS['eps_max'])**2:.2f}x of the observed "
          f"{b['peak_to_peak']/MEAS['peak_to_peak']:.2f}x amplitude excess")
    (R.OUT / "cd_steom_variants.json").write_text(
        json.dumps(dict(measured=MEAS, variants=out), indent=2))
    print(f"\n  wrote {R.OUT/'cd_steom_variants.json'}")


if __name__ == "__main__":
    main()

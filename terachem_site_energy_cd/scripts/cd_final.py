#!/usr/bin/env python3
"""The CD spectrum, two ways, and the number that survives both systematics.

`cd_rigorous.py` already produces an absolute Delta_epsilon from the 500-frame
ensemble with nothing fitted. It gets the sign, the shape, the zero crossing,
the 481 nm vibronic shoulder and the absorption FWHM (984 vs 963 cm^-1). Two
things are off: the amplitude by 3.26x and the band-top width by 1.48x.

Of the amplitude excess, 2.19 is the transition dipole. Integrating the measured
absorption band in the same cgs convention gives |mu| = 7.74 D per site, against
9.40 D from TDDFT -- and independently Cusick's two-photon work gives 7.6-7.9 D.
Our own chain closes to 0.3% (integrating the PREDICTED band returns the 9.42 D
that went in), so the discrepancy is in the electronic structure, not the units.

So this runs the model twice:

  abinitio   per-frame TDDFT dipoles, nothing substituted
  mu_exp     the same, with |mu| set to the value the measured ABSORPTION band
             implies

The second is not a fit to the CD spectrum. It substitutes one independently
measured input -- a dipole strength from absorption, cross-checked against a
two-photon determination -- into a prediction of a different observable. It is
labelled everywhere it appears and both curves are always shown together.

The first moment M1 = -P pi <nubar J T> is reported alongside. It is exact for
any detuning and any lineshape, so it is immune to the width error; M1/D^2 is
immune to the dipole error as well, and is therefore the one number quotable
with no caveat at all.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402

MU_TDDFT_D = 9.40          # ensemble geometric mean, per site
MU_EXPT_D = 7.74           # from int(eps/nubar)dnubar on the measured band
DEBYE_ESU_CM = 1.0e-18
MEAS = dict(fwhm=963.0, width90=252.0, eps_max=1.844e5,
            separation=500.0, peak_to_peak=2.829e5)


def build(scale, j, triple, mu_a, mu_b, base, sigma_site, modes):
    """One spectrum with all dipoles scaled by `scale` (J and T go as scale^2)."""
    vib_e, fc = R.vibrational_basis(modes, 3)
    nodes, wts = np.polynomial.hermite_e.hermegauss(24)
    wts = wts / wts.sum()
    de = np.zeros_like(R.GRID)
    eps = np.zeros_like(R.GRID)
    sg = np.sqrt(2.0) * sigma_site
    ma, mb, tr, jj = mu_a * scale, mu_b * scale, triple * scale**2, j * scale**2
    for i in range(j.size):
        for nd, wt in zip(nodes, wts):
            dl = sg * nd
            nu, rk, dk = R.frame_sticks(base[i] + 0.5 * dl, base[i] - 0.5 * dl,
                                        jj[i], ma[i], mb[i], tr[i], vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * rk * nu * wt, R.GRID, de)
            R.accumulate(nu, R.D_TO_EPS * dk * nu * wt, R.GRID, eps)
    de /= j.size
    eps /= j.size
    sm = sigma_site / np.sqrt(2.0)
    de = R.gaussian_convolve(de, R.GRID, sm)
    eps = R.gaussian_convolve(eps, R.GRID, sm)
    return R.ELLIPTICITY_PER_DELTA_EPSILON * de, eps


def metrics(theta, eps):
    w = lambda f: (lambda k: float(R.GRID[k[-1]] - R.GRID[k[0]]))(
        np.where(eps >= f * eps.max())[0])
    hi, lo = int(theta.argmax()), int(theta.argmin())
    return dict(peak_to_peak=float(theta.max() - theta.min()),
                separation_cm=abs(float(R.GRID[hi] - R.GRID[lo])),
                lobe_ratio=float(abs(theta.max() / theta.min())),
                eps_max=float(eps.max()), fwhm_cm=w(0.50), width90_cm=w(0.90),
                lambda_max_nm=float(1e7 / R.GRID[int(eps.argmax())]))


def main() -> None:
    d = np.load(R.ENSEMBLE)
    j, triple = d["j"], d["triple"]
    mu_a, mu_b = d["mu_a_vec"], d["mu_b_vec"]
    mean_raw = 0.5 * (d["e_a"] + d["e_b"])
    base = np.full(j.size, R.EXP_CENTER_CM)
    sigma_total = float(np.sqrt(0.5 * (d["e_a"] - d["e_b"]).var(ddof=1)))
    sigma_site, modes, _ = R.disorder_and_modes(R.site_series(), sigma_total)

    out = {"measured": MEAS, "mu_tddft_D": MU_TDDFT_D, "mu_expt_D": MU_EXPT_D,
           "cases": {}}
    curves = {}
    for tag, scale in (("abinitio", 1.0), ("mu_exp", MU_EXPT_D / MU_TDDFT_D)):
        th, ep = build(scale, j, triple, mu_a, mu_b, base, sigma_site, modes)
        m = metrics(th, ep)
        m["mu_scale"] = scale
        # first moment: exact for any Delta and any lineshape
        s = scale ** 2
        m["M1"] = float(-R.R_TO_DELTA_EPS * np.pi
                        * np.mean(R.EXP_CENTER_CM * (j * s) * (triple * s)))
        dstr = float(np.mean(np.linalg.norm(mu_a * scale, axis=1)
                             * np.linalg.norm(mu_b * scale, axis=1))
                     * DEBYE_ESU_CM ** 2)
        m["D_cgs"] = dstr
        m["M1_over_D2"] = m["M1"] / dstr ** 2
        out["cases"][tag] = m
        curves[tag] = (th, ep)

    np.savez(R.OUT / "cd_final.npz", grid_cm=R.GRID,
             **{f"{k}_theta": v[0] for k, v in curves.items()},
             **{f"{k}_eps": v[1] for k, v in curves.items()})
    (R.OUT / "cd_final.json").write_text(json.dumps(out, indent=2) + "\n")

    print(f"  |mu|: TDDFT {MU_TDDFT_D:.2f} D   measured absorption {MU_EXPT_D:.2f} D"
          f"   (Cusick two-photon 7.6-7.9)")
    print("  case       pk-pk theta   ratio   sep   eps_max    ratio  FWHM  w@0.90")
    for tag, m in out["cases"].items():
        print(f"  {tag:<10} {m['peak_to_peak']:10.3e} {m['peak_to_peak']/MEAS['peak_to_peak']:6.2f}"
              f" {m['separation_cm']:6.0f} {m['eps_max']:9.3e} {m['eps_max']/MEAS['eps_max']:7.2f}"
              f" {m['fwhm_cm']:6.0f} {m['width90_cm']:6.0f}")
    print(f"  MEASURED   {MEAS['peak_to_peak']:10.3e}   1.00 {MEAS['separation']:6.0f}"
          f" {MEAS['eps_max']:9.3e}    1.00 {MEAS['fwhm']:6.0f} {MEAS['width90']:6.0f}")
    print(f"\n  M1/D^2 (dipole- and lineshape-independent): "
          f"{out['cases']['abinitio']['M1_over_D2']:.4e}"
          f"  [identical by construction: "
          f"{out['cases']['mu_exp']['M1_over_D2']:.4e}]")
    print(f"  wrote {R.OUT/'cd_final.json'}")


if __name__ == "__main__":
    main()

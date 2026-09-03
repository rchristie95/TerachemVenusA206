#!/usr/bin/env python3
"""Determine the lineshape from the ABSORPTION band, then predict the CD.

The standing problem: our static disorder and our vibronic coupling are both
wrong, in opposite directions. sigma_site = 391 cm^-1 where the measured band
top permits ~275; the MD's reorganisation energy gives S_total = 0.278 where the
measured blue shoulder implies ~0.56. They nearly cancel in the FWHM (984 vs 963)
and do not cancel at the band top (374 vs 252), and the CD couplet separation
inherits the band-top error.

The fix is not to guess the partition but to MEASURE it, from two independent
numbers on the absorption band:

    width at 0.90 of peak   252 cm^-1   (Nguyen Fig. 3b, 76 px/nm)
    FWHM                    963 cm^-1   (Nguyen Fig. 2e)

The band top is set almost entirely by the static Gaussian and the half width
almost entirely by the vibronic sidebands, so the two observables separate the
two unknowns. Two measurements, two unknowns: determined, not fitted.

This is legitimate because ABSORPTION IS A DIFFERENT OBSERVABLE FROM CD. Nothing
here is tuned against the couplet. Once sigma_static and the mode displacement
are fixed by the absorption band, the couplet separation, its amplitude and its
first moment are predictions, and the couplet separation in particular is the
one that was 1.26x off.

The mode FREQUENCIES stay where the measured spectral density put them
(667/1220/1731 cm^-1); only their overall displacement is scaled, by a single
factor, so the shape of the progression is not a free function.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402

MU_SCALE = 7.74 / 9.40          # measured absorption dipole, as in cd_final.py
TARGET_W90 = 252.0
TARGET_FWHM = 963.0
MODES0 = [(667.1, 0.0818), (1219.9, 0.1262), (1731.2, 0.0699)]
LAMBDA_MD = sum(w * s for w, s in MODES0)      # 329.5 cm^-1 from the classical MD


def spectrum(sigma_site, s_scale, data, n_frames, n_hermite):
    j, triple, mu_a, mu_b, base = data
    sl = slice(0, n_frames)
    j, triple, mu_a, mu_b, base = (j[sl], triple[sl], mu_a[sl], mu_b[sl], base[sl])
    modes = [(w, s * s_scale) for w, s in MODES0]
    vib_e, fc = R.vibrational_basis(modes, 3)
    nodes, wts = np.polynomial.hermite_e.hermegauss(n_hermite)
    wts = wts / wts.sum()
    de = np.zeros_like(R.GRID)
    eps = np.zeros_like(R.GRID)
    sg = np.sqrt(2.0) * sigma_site
    ma, mb = mu_a * MU_SCALE, mu_b * MU_SCALE
    tr, jj = triple * MU_SCALE**2, j * MU_SCALE**2
    for i in range(j.size):
        for nd, wt in zip(nodes, wts):
            dl = sg * nd
            nu, rk, dk = R.frame_sticks(base[i] + .5*dl, base[i] - .5*dl,
                                        jj[i], ma[i], mb[i], tr[i], vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * rk * nu * wt, R.GRID, de)
            R.accumulate(nu, R.D_TO_EPS * dk * nu * wt, R.GRID, eps)
    de /= j.size
    eps /= j.size
    sm = sigma_site / np.sqrt(2.0)
    de = R.gaussian_convolve(de, R.GRID, sm)
    eps = R.gaussian_convolve(eps, R.GRID, sm)
    return R.ELLIPTICITY_PER_DELTA_EPSILON * de, eps


def widths(eps):
    w = lambda f: (lambda k: float(R.GRID[k[-1]] - R.GRID[k[0]]))(
        np.where(eps >= f * eps.max())[0])
    return w(0.90), w(0.50)


def main() -> None:
    d = np.load(R.ENSEMBLE)
    base = np.full(d["j"].size, R.EXP_CENTER_CM)
    data = (d["j"], d["triple"], d["mu_a_vec"], d["mu_b_vec"], base)

    print(f"  targets: width@0.90 {TARGET_W90:.0f}, FWHM {TARGET_FWHM:.0f} cm^-1")
    print(f"  MD gives lambda_site = {LAMBDA_MD:.0f} cm^-1 (S_total = "
          f"{sum(s for _, s in MODES0):.3f})\n")
    print("  sigma  Sscale  lambda   w@0.90    FWHM")
    grid = []
    for sigma in (200., 250., 275., 300., 350.):
        for sc in (1.0, 1.5, 2.0, 2.5):
            th, ep = spectrum(sigma, sc, data, 60, 12)
            w90, fw = widths(ep)
            grid.append((sigma, sc, w90, fw))
            print(f"  {sigma:5.0f}  {sc:5.2f}  {sc*LAMBDA_MD:6.0f}  "
                  f"{w90:7.0f} {fw:7.0f}")

    # solve the 2x2 by local linear inversion around the closest grid point
    a = np.array(grid)
    err = np.hypot((a[:, 2] - TARGET_W90) / TARGET_W90,
                   (a[:, 3] - TARGET_FWHM) / TARGET_FWHM)
    s0, c0 = a[err.argmin(), :2]
    for _ in range(6):
        h_s, h_c = 15.0, 0.15
        f0 = np.array(widths(spectrum(s0, c0, data, 60, 12)[1]))
        js = (np.array(widths(spectrum(s0 + h_s, c0, data, 60, 12)[1])) - f0) / h_s
        jc = (np.array(widths(spectrum(s0, c0 + h_c, data, 60, 12)[1])) - f0) / h_c
        J = np.column_stack([js, jc])
        step = np.linalg.lstsq(J, np.array([TARGET_W90, TARGET_FWHM]) - f0,
                               rcond=None)[0]
        s0 = float(np.clip(s0 + step[0], 50, 500))
        c0 = float(np.clip(c0 + step[1], 0.2, 4.0))
    print(f"\n  solved: sigma_site = {s0:.0f} cm^-1, mode scale = {c0:.2f} "
          f"-> lambda_site = {c0*LAMBDA_MD:.0f} cm^-1 (S_total = "
          f"{c0*sum(s for _,s in MODES0):.3f})")
    print(f"  the classical MD under-counts the reorganisation energy by "
          f"{c0:.2f}x and over-counts the static disorder by {391.5/s0:.2f}x")

    # full ensemble at the solution: the CD is now a PREDICTION
    th, ep = spectrum(s0, c0, data, d["j"].size, 24)
    w90, fw = widths(ep)
    hi, lo = int(th.argmax()), int(th.argmin())
    sep = abs(float(R.GRID[hi] - R.GRID[lo]))
    C = R.EXP_CENTER_CM
    L = lambda A, c, w: A / (1 + ((R.GRID - c) / w) ** 2)
    obs = sum(L(-A*1e5, c, w) for A, c, w in ((3.59, C-130.79, 293.51),
              (-2.27, C+130.79, 439.25), (-0.34, 20774.0, 565.79)))
    mom = lambda y: np.trapezoid(y * (R.GRID - C) / R.GRID, R.GRID)
    out = dict(sigma_site_cm=s0, mode_scale=c0, lambda_site_cm=c0*LAMBDA_MD,
               width90_cm=w90, fwhm_cm=fw, separation_cm=sep,
               peak_to_peak=float(th.max()-th.min()),
               lobe_ratio=float(abs(th.max()/th.min())),
               eps_max=float(ep.max()), M1_ratio=float(mom(th)/mom(obs)))
    print(f"\n  FITTED to absorption:  w@0.90 {w90:.0f} (252)   FWHM {fw:.0f} (963)"
          f"   eps_max {ep.max():.3e} (1.844e5)")
    print(f"  PREDICTED for the CD:  separation {sep:.0f} cm^-1 (measured 500)"
          f"   pk-pk {out['peak_to_peak']:.3e} (2.829e5)")
    print(f"                         lobe ratio {out['lobe_ratio']:.2f} (0.43)"
          f"   M1 ratio {out['M1_ratio']:.2f}")
    np.savez(R.OUT / "cd_solved.npz", grid_cm=R.GRID, theta=th, eps=ep, observed=obs)
    (R.OUT / "cd_solved.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {R.OUT/'cd_solved.json'}")


if __name__ == "__main__":
    main()

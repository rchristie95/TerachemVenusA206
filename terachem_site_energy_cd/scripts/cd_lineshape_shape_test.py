#!/usr/bin/env python3
"""What lineshape does the measured Venus band actually demand?

cd_rigorous.py asserts there is NO homogeneous width: every classical mode is
either below kT (quasi-static, so inhomogeneous) or quantum-frozen (so
Franck-Condon).  That assertion is testable against three independent numbers
read off Nguyen's own absorption data:

    width at 0.90 of peak   252 cm^-1   (Fig. 3b, digitised, 76 px/nm)
    FWHM                    963 cm^-1   (Fig. 2e, digitised)
    CD couplet separation   500 cm^-1   (Table S3 fit)

A purely Gaussian static disorder cannot satisfy the first two at once: the
band top demands sigma ~ 275 cm^-1, the half width demands ~392.  A sharp core
with broad wings is the signature of a Lorentzian component, i.e. of exactly the
homogeneous broadening the model sets to zero.  This scans (sigma_static, gamma)
to find whether any combination satisfies all three, and what gamma that implies.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402

MEAS = dict(width90_cm=252.0, fwhm_cm=963.0, separation_cm=500.0)
MODES = [(667.1, 0.0818), (1219.9, 0.1262), (1731.2, 0.0699)]


def lorentz_convolve(spec, grid, gamma):
    if gamma <= 0.0:
        return spec
    dx = grid[1] - grid[0]
    half = min(int(np.ceil(40.0 * gamma / dx)), grid.size // 2 - 1)
    x = np.arange(-half, half + 1) * dx
    k = gamma / np.pi / (x ** 2 + gamma ** 2)
    k /= k.sum()
    return np.convolve(spec, k, mode="same")


def run(sigma_site, gamma, j, triple, mu_a, mu_b, base, n_hermite=16):
    vib_e, fc = R.vibrational_basis(MODES, 3)
    nodes, wts = np.polynomial.hermite_e.hermegauss(n_hermite)
    wts = wts / wts.sum()
    if sigma_site <= 0.0:
        nodes, wts = np.array([0.0]), np.array([1.0])
    de = np.zeros_like(R.GRID)
    ep = np.zeros_like(R.GRID)
    sg = np.sqrt(2.0) * sigma_site
    for i in range(j.size):
        for nd, wt in zip(nodes, wts):
            dl = sg * nd
            nu, rk, dk = R.frame_sticks(base[i] + 0.5 * dl, base[i] - 0.5 * dl,
                                        j[i], mu_a[i], mu_b[i], triple[i],
                                        vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * rk * nu * wt, R.GRID, de)
            R.accumulate(nu, R.D_TO_EPS * dk * nu * wt, R.GRID, ep)
    de /= j.size
    ep /= j.size
    sm = sigma_site / np.sqrt(2.0)
    if sm > 0:
        de = R.gaussian_convolve(de, R.GRID, sm)
        ep = R.gaussian_convolve(ep, R.GRID, sm)
    de = lorentz_convolve(de, R.GRID, gamma)
    ep = lorentz_convolve(ep, R.GRID, gamma)
    th = R.ELLIPTICITY_PER_DELTA_EPSILON * de
    w = lambda f: (lambda k: float(R.GRID[k[-1]] - R.GRID[k[0]]))(
        np.where(ep >= f * ep.max())[0])
    return dict(sigma=sigma_site, gamma=gamma,
                width90_cm=w(0.90), fwhm_cm=w(0.50),
                separation_cm=abs(float(R.GRID[int(th.argmax())]
                                        - R.GRID[int(th.argmin())])),
                peak_to_peak=float(th.max() - th.min()),
                eps_max=float(ep.max()),
                lobe_ratio=float(abs(th.max() / th.min())))


def main() -> None:
    d = np.load(R.ENSEMBLE)
    j, triple = d["j"], d["triple"]
    mu_a, mu_b = d["mu_a_vec"], d["mu_b_vec"]
    base = np.full(j.size, R.EXP_CENTER_CM)

    print("  measured:  width@0.90 252   FWHM 963   couplet sep 500   "
          "eps_max 1.844e+05   lobe 0.43\n")
    print("  sigma  gamma   width@0.90    FWHM   couplet sep    eps_max  lobe")
    rows = []
    for sigma in (0.0, 100.0, 150.0, 200.0, 275.0, 391.5):
        for gamma in (0.0, 100.0, 200.0, 300.0):
            r = run(sigma, gamma, j, triple, mu_a, mu_b, base)
            rows.append(r)
            flag = ""
            if abs(r["width90_cm"] - MEAS["width90_cm"]) < 30 and \
               abs(r["fwhm_cm"] - MEAS["fwhm_cm"]) < 90:
                flag = "  <== matches both widths"
            print(f"  {sigma:5.0f}  {gamma:5.0f}   {r['width90_cm']:9.0f} "
                  f"{r['fwhm_cm']:7.0f} {r['separation_cm']:12.0f}  "
                  f"{r['eps_max']:9.3e} {r['lobe_ratio']:5.2f}{flag}")

    (R.OUT / "cd_lineshape_shape_test.json").write_text(
        json.dumps(dict(measured=MEAS, scan=rows), indent=2))
    print(f"\n  wrote {R.OUT/'cd_lineshape_shape_test.json'}")


if __name__ == "__main__":
    main()

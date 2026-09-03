#!/usr/bin/env python3
"""Why is the predicted CD couplet 2x too wide?

The couplet separation of a bisignate band is set by the curvature of the
underlying absorption band at its maximum, not by the exciton gap.  Nguyen
Fig. 3b measures exactly that curvature for dVenus-TD: the top 10% of the
normalised absorption band, giving a width at 0.90 of peak of 6.745 nm =
252 cm^-1 (digitised, calibrated on their own stated peak of 516.82 nm).

That is the like-for-like quantity to compare, and it is model free -- no
Gaussian or Lorentzian extrapolation to a half width is involved.

This decomposes our predicted band into its three broadening terms, compares
each against that measurement, and then inverts: what static disorder would the
measured band permit, and what couplet separation would follow.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402

MEAS_WIDTH90_CM = 252.0            # dVenus-TD, Nguyen Fig. 3b, digitised
MEAS_PEAK_NM = 516.82


def width_at(grid, spec, frac):
    idx = np.where(spec >= frac * spec.max())[0]
    return float(grid[idx[-1]] - grid[idx[0]])


def build(sigma_site, s_scale, modes, j, triple, mu_a, mu_b, base,
          n_hermite=16):
    scaled = [(w, s * s_scale) for w, s in modes]
    vib_e, fc = R.vibrational_basis(scaled, R.N_QUANTA)
    sig_gap = np.sqrt(2.0) * sigma_site
    nodes, wts = np.polynomial.hermite_e.hermegauss(n_hermite)
    wts = wts / wts.sum()
    if sigma_site == 0.0:
        nodes, wts = np.array([0.0]), np.array([1.0])
    d_eps = np.zeros_like(R.GRID)
    eps = np.zeros_like(R.GRID)
    for i in range(j.size):
        for nd, wt in zip(nodes, wts):
            delta = sig_gap * nd
            nu, r_k, d_k = R.frame_sticks(base[i] + 0.5 * delta,
                                          base[i] - 0.5 * delta, j[i],
                                          mu_a[i], mu_b[i], triple[i], vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * r_k * nu * wt, R.GRID, d_eps)
            R.accumulate(nu, R.D_TO_EPS * d_k * nu * wt, R.GRID, eps)
    d_eps /= j.size
    eps /= j.size
    if sigma_site > 0.0:
        s_mean = sigma_site / np.sqrt(2.0)
        d_eps = R.gaussian_convolve(d_eps, R.GRID, s_mean)
        eps = R.gaussian_convolve(eps, R.GRID, s_mean)
    theta = R.ELLIPTICITY_PER_DELTA_EPSILON * d_eps
    return theta, eps


def report(tag, theta, eps):
    hi, lo = int(np.argmax(theta)), int(np.argmin(theta))
    return dict(tag=tag,
                width90_cm=width_at(R.GRID, eps, 0.90),
                fwhm_cm=width_at(R.GRID, eps, 0.50),
                separation_cm=abs(float(R.GRID[hi] - R.GRID[lo])),
                peak_to_peak=float(theta.max() - theta.min()))


def main() -> None:
    d = np.load(R.ENSEMBLE)
    e_a, e_b, j, triple = d["e_a"], d["e_b"], d["j"], d["triple"]
    mu_a, mu_b = d["mu_a_vec"], d["mu_b_vec"]
    mean_raw = 0.5 * (e_a + e_b)
    base = mean_raw + (R.EXP_CENTER_CM - mean_raw.mean())
    flat = np.full_like(base, R.EXP_CENTER_CM)     # removes ensemble mean spread

    series = R.site_series()
    sigma_total = float(np.sqrt(0.5 * (e_a - e_b).var(ddof=1)))
    sigma_static, modes, diag = R.disorder_and_modes(series, sigma_total)

    print("=== what the measured band allows ===")
    print(f"  dVenus-TD width at 0.90 of peak, Nguyen Fig. 3b:  "
          f"{MEAS_WIDTH90_CM:.0f} cm^-1")
    print(f"  measured couplet separation:                      500 cm^-1\n")

    print("=== decomposition of our predicted band ===")
    rows = []
    for tag, sig, sc, bs in (
            ("vibronic only",                     0.0, 1.0, flat),
            ("+ ensemble mean spread",            0.0, 1.0, base),
            ("+ static disorder (full model)",    sigma_static, 1.0, base),
            ("static disorder only, no vibronic", sigma_static, 0.0, base)):
        th, ep = build(sig, sc, modes, j, triple, mu_a, mu_b, bs)
        r = report(tag, th, ep)
        rows.append(r)
        print(f"  {tag:<34} width@0.90 {r['width90_cm']:6.0f}  FWHM {r['fwhm_cm']:6.0f}"
              f"  couplet sep {r['separation_cm']:6.0f}")

    print(f"\n  each term alone, against the measured {MEAS_WIDTH90_CM:.0f} cm^-1:")
    print(f"    vibronic progression only        "
          f"{rows[0]['width90_cm']/MEAS_WIDTH90_CM:5.2f}x")
    print(f"    static disorder only             "
          f"{rows[3]['width90_cm']/MEAS_WIDTH90_CM:5.2f}x")
    print(f"    full model                       "
          f"{rows[2]['width90_cm']/MEAS_WIDTH90_CM:5.2f}x")

    print("\n=== inversion: what disorder does the measured band permit? ===")
    print("  sigma_site  width@0.90   FWHM   couplet sep   pk-pk theta")
    scan = []
    for sig in (0.0, 100.0, 150.0, 200.0, 250.0, 300.0, sigma_static):
        th, ep = build(sig, 1.0, modes, j, triple, mu_a, mu_b, base)
        r = report(f"sigma={sig:.0f}", th, ep)
        r["sigma_site"] = sig
        scan.append(r)
        print(f"  {sig:9.0f}  {r['width90_cm']:9.0f} {r['fwhm_cm']:7.0f} "
              f"{r['separation_cm']:12.0f}  {r['peak_to_peak']:12.3e}")

    print("\n=== inversion: what vibronic coupling does it permit? ===")
    print("  S scale  lambda_site  width@0.90   FWHM   couplet sep")
    for sc in (0.0, 0.25, 0.5, 0.75, 1.0):
        th, ep = build(sigma_static, sc, modes, j, triple, mu_a, mu_b, base)
        r = report(f"S x{sc}", th, ep)
        lam = sc * sum(w * s for w, s in modes)
        print(f"  {sc:7.2f}  {lam:11.0f}  {r['width90_cm']:9.0f} {r['fwhm_cm']:7.0f} "
              f"{r['separation_cm']:12.0f}")

    (R.OUT / "cd_width_diagnosis.json").write_text(json.dumps(
        dict(measured_width90_cm=MEAS_WIDTH90_CM,
             measured_peak_nm=MEAS_PEAK_NM,
             sigma_site_static_cm=sigma_static,
             decomposition=rows, disorder_scan=scan), indent=2))
    print(f"\n  wrote {R.OUT/'cd_width_diagnosis.json'}")


if __name__ == "__main__":
    main()

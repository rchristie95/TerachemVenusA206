#!/usr/bin/env python3
"""Sensitivity of the rigorous CD to the single remaining discretionary choice.

cd_rigorous.py partitions the measured spectral density at 2kT: bands below are
quasi-static disorder on the absorption timescale, bands above are Franck-Condon
modes.  Everything else in that model is measured.  This scans the boundary over
kT to 4kT (the manuscript's own quantum-frozen threshold) and, separately, the
number of vibrational quanta, so the residual model dependence is stated rather
than assumed small.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cd_rigorous as R                                            # noqa: E402


def run(cutoff_factor, n_quanta, n_hermite):
    R.KT_CM = 208.5
    saved = R.disorder_and_modes.__globals__["KT_CM"]
    d = np.load(R.ENSEMBLE)
    e_a, e_b, j, triple = d["e_a"], d["e_b"], d["j"], d["triple"]
    mu_a_vec, mu_b_vec = d["mu_a_vec"], d["mu_b_vec"]
    mean_raw = 0.5 * (e_a + e_b)
    sigma_site_total = float(np.sqrt(0.5 * (e_a - e_b).var(ddof=1)))

    series = R.site_series()
    n = series[0].size
    win = np.hanning(n)
    freqs = np.fft.rfftfreq(n, 1e-15) / R.C_CM_S
    power = np.zeros(freqs.size)
    for s in series:
        power += np.abs(np.fft.rfft((s - s.mean()) * win)) ** 2
    power /= power.sum()
    var_fast = float(np.mean([s.var(ddof=1) for s in series]))
    var_slow = sigma_site_total ** 2 - var_fast

    modes, quasi = [], 0.0
    edges = (0.0,) + R.BAND_EDGES + (float(freqs[-1]) + 1.0,)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (freqs >= lo) & (freqs < hi)
        if not m.any():
            continue
        frac = float(power[m].sum())
        if frac < 1e-3:
            continue
        omega = float((freqs[m] * power[m]).sum() / power[m].sum())
        var_i = frac * var_fast
        s_i = (var_i / (2.0 * saved)) / omega
        if omega <= cutoff_factor * saved or s_i < R.S_MIN:
            quasi += var_i
        else:
            modes.append((omega, s_i))
    sigma_static = float(np.sqrt(var_slow + quasi))

    vib_e, fc = R.vibrational_basis(modes, n_quanta)
    nodes, weights = np.polynomial.hermite_e.hermegauss(n_hermite)
    weights = weights / weights.sum()
    base = mean_raw + (R.EXP_CENTER_CM - mean_raw.mean())
    sig_gap = np.sqrt(2.0) * sigma_static

    d_eps = np.zeros_like(R.GRID)
    eps = np.zeros_like(R.GRID)
    for i in range(j.size):
        for node, wt in zip(nodes, weights):
            delta = sig_gap * node
            nu, r_k, d_k = R.frame_sticks(base[i] + 0.5 * delta,
                                          base[i] - 0.5 * delta, j[i],
                                          mu_a_vec[i], mu_b_vec[i], triple[i],
                                          vib_e, fc)
            R.accumulate(nu, R.R_TO_DELTA_EPS * r_k * nu * wt, R.GRID, d_eps)
            R.accumulate(nu, R.D_TO_EPS * d_k * nu * wt, R.GRID, eps)
    d_eps = R.gaussian_convolve(d_eps / j.size, R.GRID, sigma_static / np.sqrt(2))
    eps = R.gaussian_convolve(eps / j.size, R.GRID, sigma_static / np.sqrt(2))
    theta = R.ELLIPTICITY_PER_DELTA_EPSILON * d_eps
    hi, lo = int(np.argmax(theta)), int(np.argmin(theta))
    half = eps.max() / 2.0
    k = np.where(eps >= half)[0]
    return dict(cutoff_factor=cutoff_factor, n_quanta=n_quanta,
                n_modes=len(modes), sigma_site_static_cm=sigma_static,
                separation_cm=abs(float(R.GRID[hi] - R.GRID[lo])),
                peak_to_peak=float(theta.max() - theta.min()),
                eps_max=float(eps.max()),
                eps_fwhm_cm=float(R.GRID[k[-1]] - R.GRID[k[0]]),
                g_abs=float(np.abs(d_eps).max() / eps.max()))


def main() -> None:
    rows = []
    print("  cutoff  quanta  modes  sigma_stat  couplet sep  pk-pk theta"
          "   eps_max   abs FWHM    g_abs")
    for cf, nq in ((1.0, 3), (2.0, 3), (4.0, 3), (2.0, 2), (2.0, 4)):
        r = run(cf, nq, 12)
        rows.append(r)
        print(f"  {cf:4.1f}kT  {nq:5d}  {r['n_modes']:5d}  "
              f"{r['sigma_site_static_cm']:9.1f}  {r['separation_cm']:10.0f}  "
              f"{r['peak_to_peak']:11.3e}  {r['eps_max']:9.3e}  "
              f"{r['eps_fwhm_cm']:8.0f}  {r['g_abs']:8.2e}")
    (R.OUT / "cd_rigorous_sensitivity.json").write_text(json.dumps(rows, indent=2))
    sep = [r["separation_cm"] for r in rows[:3]]
    g = [r["g_abs"] for r in rows[:3]]
    print(f"\n  across the kT-4kT boundary: separation {min(sep):.0f}-{max(sep):.0f}"
          f" cm^-1, g_abs {min(g):.2e}-{max(g):.2e}")
    print(f"  measured: separation 500 cm^-1, g_abs 3.26e-04")


if __name__ == "__main__":
    main()

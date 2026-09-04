#!/usr/bin/env python3
"""Cumulant lineshape from the 1 ps dense segment, and the disorder split it fixes.

The three v6t2 windows are 100 fs and resolve only to 334 cm^-1 -- coarser than
kT = 208 -- so the fast/slow partition behind the CD band shape has until now
been made at a boundary the data could not see. v8dense is 1 ps at 2 fs, which
resolves to 33 cm^-1, and this is the first trajectory that can place the
boundary rather than assume it.

Reports, per site and for the gap:
  * the correlation function C(t) and its correlation time;
  * the Kubo parameter sigma*tau_c, which says whether the fast part is in the
    motional-narrowing, intermediate or static limit;
  * the second-order cumulant lineshape and its width;
  * the variance below and above kT, which IS the fast/slow split, measured;
  * how that compares with the split cd_rigorous.py currently uses
    (sigma_site_static = 391.5 cm^-1, S_total = 0.278).

Read-only apart from results/v8dense_cumulant.json.
"""
from __future__ import annotations

import csv, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "terachem_site_energy_cd"))
from cumulant_lineshape import cumulant_g, lineshape_from_g, fwhm_of   # noqa: E402

KT_CM, C_CM_S, DT_FS = 208.5, 29979245800.0, 2.0
CURRENT_SIGMA_STATIC, CURRENT_S_TOTAL = 391.5, 0.278


def series():
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/v8dense_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))
    a = np.array([float(r["site_A_energy_cm-1"]) for r in rows])
    b = np.array([float(r["site_B_energy_cm-1"]) for r in rows])
    return a, b, len(rows)


def acf(x):
    x = x - x.mean()
    n = x.size
    f = np.fft.rfft(x, 2 * n)
    c = np.fft.irfft(f * np.conj(f), 2 * n)[:n].real
    return c / np.arange(n, 0, -1)          # unbiased


def main() -> None:
    a, b, n = series()
    lag = np.arange(n) * DT_FS
    out = {"n_frames": n, "dt_fs": DT_FS, "span_ps": n * DT_FS / 1000.0,
           "resolution_cm": 1.0 / (n * DT_FS * 1e-15 * C_CM_S)}
    print(f"  {n} frames at {DT_FS} fs = {out['span_ps']:.2f} ps, "
          f"resolving to {out['resolution_cm']:.0f} cm^-1 "
          f"(v6t2 reached only 334)")

    for tag, x in (("site_A", a), ("site_B", b), ("gap", a - b)):
        c = acf(x)
        sd = float(x.std(ddof=1))
        tau = float(np.trapezoid(c[:c.size // 2], lag[:c.size // 2]) / c[0]) if c[0] else 0.0
        # cumulant_g takes C in cm^-2 and applies TWO_PI_C**2 itself; passing a
        # pre-converted C squares the conversion, makes g ~1e-11, and returns the
        # 1/T window function (about 20 cm^-1 here) instead of a lineshape.
        g = cumulant_g(lag, c)
        grid = np.arange(-4000.0, 4000.0, 2.0)
        L = lineshape_from_g(lag, g, grid, 0.0)
        fw = fwhm_of(grid, L)
        # variance below / above kT from the power spectrum
        freqs = np.fft.rfftfreq(x.size, DT_FS * 1e-15) / C_CM_S
        p = np.abs(np.fft.rfft((x - x.mean()) * np.hanning(x.size))) ** 2
        p /= p.sum()
        below = float(p[freqs <= KT_CM].sum())
        rec = dict(sd_cm=sd, tau_c_fs=tau, kubo_sigma_tau=float(
            sd * 2 * np.pi * C_CM_S * 1e-15 * tau), cumulant_fwhm_cm=fw,
            var_fraction_below_kT=below,
            sd_below_kT_cm=float(sd * np.sqrt(below)),
            sd_above_kT_cm=float(sd * np.sqrt(1 - below)))
        out[tag] = rec
        print(f"  {tag:7s} sd {sd:6.1f}  tau_c {tau:5.1f} fs  sigma*tau {rec['kubo_sigma_tau']:.2f}"
              f"  cumulant FWHM {fw:6.0f} cm^-1"
              f"  |  below kT {100*below:4.1f}% -> {rec['sd_below_kT_cm']:5.1f},"
              f" above -> {rec['sd_above_kT_cm']:5.1f}")

    print(f"\n  the split cd_rigorous.py uses: sigma_static {CURRENT_SIGMA_STATIC:.1f} cm^-1, "
          f"S_total {CURRENT_S_TOTAL:.3f}")
    print(f"  measured here (site rms, below kT = quasi-static):        "
          f"{np.hypot(out['site_A']['sd_below_kT_cm'], out['site_B']['sd_below_kT_cm'])/np.sqrt(2):.1f} cm^-1")
    lam = (np.hypot(out['site_A']['sd_above_kT_cm'],
                    out['site_B']['sd_above_kT_cm']) / np.sqrt(2)) ** 2 / (2 * KT_CM)
    print(f"  reorganisation energy from the above-kT part:             {lam:.0f} cm^-1")
    out["comparison"] = {"cd_rigorous_sigma_static_cm": CURRENT_SIGMA_STATIC,
                         "cd_rigorous_S_total": CURRENT_S_TOTAL,
                         "measured_lambda_site_cm": float(lam)}
    (ROOT / "results" / "v8dense_cumulant.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {ROOT/'results'/'v8dense_cumulant.json'}")


if __name__ == "__main__":
    main()

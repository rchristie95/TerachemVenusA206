#!/usr/bin/env python3
"""Homogeneous linewidth from the REAL QM/MM gap autocorrelation.

Item 1 of the CD diagnosis. Everything about whether our site-energy disorder is
too large hinges on the homogeneous linewidth, and the two available estimates
disagree 7x:

  proxy (first-order Delta_q . V electrostatic gap, 4 fs/8 ps)  601 cm^-1 FWHM
  Nguyen Table S4 coherence time, dVenus-TD 128 fs               83 cm^-1 FWHM

The proxy's own docstring warns it underestimates the amplitude while the
correlation TIME is transferable, so it cannot settle this. This script uses the
actual QM/MM TDDFT site energies from the 8 ps / 4 fs high-cadence campaign
(tddft_ensemble_decoherence/manifests/short_frames.csv), of which 193 frames
completed -- 114 of them contiguous from t=0, i.e. 456 fs. That is far longer
than either candidate decay time, so it discriminates cleanly.

If the real linewidth is ~600 cm^-1, the observed couplet width is entirely
homogeneous, leaving no room for static disorder and forcing <|Delta|> <~ 115
against our computed 576 -- which would undermine the detuning-dominated
mechanism. If it is ~83 cm^-1, the couplet width is mostly inhomogeneous and
our 576 is broadly consistent.

Read-only apart from results/real_qmmm_linewidth.json. No network.
"""
import csv, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cumulant_lineshape import cumulant_g, lineshape_from_g, fwhm_of  # noqa: E402

MAN = ROOT.parent / "tddft_ensemble_decoherence" / "manifests" / "short_frames.csv"
DT_FS = 4.0


def contiguous(rows):
    idx = np.array(sorted(int(r["frame_index"]) for r in rows))
    best = cur = 1; bs = cs = idx[0]
    for a, b in zip(idx, idx[1:]):
        if b - a == 1:
            cur += 1
        else:
            if cur > best: best, bs = cur, cs
            cur = 1; cs = b
    if cur > best: best, bs = cur, cs
    return bs, best


def acf(x, nlag):
    x = x - x.mean(); n = len(x)
    return np.array([np.dot(x[:n-k], x[k:]) / (n - k) for k in range(nlag)])


def main():
    rows = [r for r in csv.DictReader(open(MAN)) if r["convergence_status"] == "complete"]
    start, run = contiguous(rows)
    keep = {start + i for i in range(run)}
    sel = sorted((r for r in rows if int(r["frame_index"]) in keep),
                 key=lambda r: int(r["frame_index"]))
    ea = np.array([float(r["site_A_energy_cm-1"]) for r in sel])
    eb = np.array([float(r["site_B_energy_cm-1"]) for r in sel])
    gap = ea - eb
    print(f"REAL QM/MM TDDFT, {len(sel)} contiguous frames from index {start}"
          f"  ({len(sel)*DT_FS:.0f} fs at {DT_FS:.0f} fs)\n")
    print(f"  site A  {ea.mean():9.1f} +- {ea.std(ddof=1):6.1f} cm^-1")
    print(f"  site B  {eb.mean():9.1f} +- {eb.std(ddof=1):6.1f} cm^-1")
    print(f"  gap     {gap.mean():+9.1f} +- {gap.std(ddof=1):6.1f} cm^-1"
          f"   <|gap|> = {np.abs(gap).mean():.1f}")
    r = np.corrcoef(ea, eb)[0, 1]
    print(f"  site-site correlation r = {r:+.3f}"
          f"   ({'correlated -> gap narrower than sqrt(2)*sigma' if r > 0.2 else 'uncorrelated'})")

    nlag = len(gap) // 2
    lag = np.arange(nlag) * DT_FS
    grid = np.arange(-3000., 3000., 1.0)
    print(f"\n{'channel':>10} {'sigma':>8} {'FWHM':>9} {'HWHM':>9}")
    out = {}
    for key, series in (("site A", ea), ("site B", eb), ("gap", gap)):
        C = acf(series, nlag)
        g = cumulant_g(lag, C, amplitude_scale=1.0)
        y = lineshape_from_g(lag, g, grid, 0.0)
        f = fwhm_of(grid, y)
        out[key] = {"sigma_cm": float(series.std(ddof=1)), "fwhm_cm": float(f)}
        print(f"{key:>10} {series.std(ddof=1):8.1f} {f:9.1f} {f/2:9.1f}  cm^-1")

    site_hwhm = 0.5 * (out["site A"]["fwhm_cm"] + out["site B"]["fwhm_cm"]) / 2
    print(f"\n  REAL QM/MM homogeneous HWHM (mean of sites) = {site_hwhm:.1f} cm^-1")
    print(f"  proxy (Delta_q . V)                          = 300.5 cm^-1")
    print(f"  Nguyen Table S4 (T2 = 128 fs)                =  41.5 cm^-1")
    json.dump({"n_frames": len(sel), "dt_fs": DT_FS, "channels": out,
               "real_site_hwhm_cm": site_hwhm, "proxy_hwhm_cm": 300.5,
               "nguyen_T2_hwhm_cm": 41.5},
              open(ROOT / "results" / "real_qmmm_linewidth.json", "w"), indent=2)
    print("\nwrote results/real_qmmm_linewidth.json")


if __name__ == "__main__":
    main()

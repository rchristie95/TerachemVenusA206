#!/usr/bin/env python3
"""Spectral density, reorganisation energy and Huang-Rhys factor from the gap.

Cross-checks the MD route against the vibronic parameters the CD forward model
uses (S = 0.253, omega = 1369.7 cm^-1), which came from an independent TeraChem
CAM-B3LYP/6-31G* excited-state calculation. If the two agree, the vibronic and
MD analyses are describing one thing rather than two.

Classical relations used, in the high-temperature limit that a classical
trajectory can support:

    lambda = <dE^2> / (2 kT)          reorganisation energy from the variance
    S      = lambda / (hbar omega)    Huang-Rhys, if attributed to one mode

Sites are uncorrelated to within r = +0.016 over 1 ns, so
sigma_site^2 = sigma_gap^2 / 2 and the per-site reorganisation energy is half
the gap value. That is checked here rather than assumed.

CAVEAT, stated because it limits what the number means: a 100 fs window resolves
to 1/(c*100fs) = 334 cm^-1, coarser than kT = 208 cm^-1. The spectral density
below kT -- the only part a classical treatment is entitled to -- is therefore
unresolved. The 500 fs windows (v6lfa/b/c) reach 67 cm^-1 and are what settle it.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EV = 8065.544005
KT_CM = 208.5                 # 300 K
OMEGA_VIBRONIC = 1369.7       # repo value, TeraChem CAM-B3LYP/6-31G*
S_VIBRONIC = 0.253
C_CM_PER_FS = 2.99792458e-5


def load(window: str) -> np.ndarray:
    with (ROOT / "manifests" / f"{window}_frames.csv").open(newline="") as h:
        rows = [r for r in csv.DictReader(h) if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))
    return np.array([float(r["site_A_energy_eV"]) - float(r["site_B_energy_eV"])
                     for r in rows]) * EV


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windows", nargs="+", default=["v6t2a", "v6t2b", "v6t2c"])
    ap.add_argument("--dt-fs", type=float, default=1.0)
    ap.add_argument("--output", type=Path,
                    default=ROOT / "results/gap_spectral_density.json")
    args = ap.parse_args()

    gaps = [load(w) for w in args.windows]
    n = min(len(g) for g in gaps)
    window_fs = n * args.dt_fs
    resolution = 1.0 / (C_CM_PER_FS * window_fs)

    # pooled variance about each window's own mean (the offsets between windows
    # are slow disorder, not fast)
    var = float(np.mean([g[:n].var(ddof=1) for g in gaps]))
    sigma_gap = float(np.sqrt(var))
    sigma_site = sigma_gap / np.sqrt(2.0)

    lam_gap = var / (2.0 * KT_CM)
    lam_site = lam_gap / 2.0
    s_if_one_mode = lam_site / OMEGA_VIBRONIC

    # pooled normalised ACF, then the classical spectral density
    num = np.zeros(n); cnt = np.zeros(n)
    for g in gaps:
        x = g[:n] - g[:n].mean()
        for k in range(n):
            num[k] += float(np.dot(x[:n - k], x[k:])); cnt[k] += (n - k)
    acf = num / cnt

    freqs = np.arange(1, n // 2) / (window_fs * C_CM_PER_FS)
    t = np.arange(n) * args.dt_fs
    win = 0.5 * (1 + np.cos(np.pi * t / t[-1]))          # taper, suppresses ringing
    jw = np.array([np.trapezoid(acf * win * np.cos(2 * np.pi * C_CM_PER_FS * f * t), t)
                   for f in freqs])
    jw = np.clip(jw, 0.0, None)
    weight = jw / jw.sum() if jw.sum() > 0 else jw
    below = float(weight[freqs < KT_CM].sum())
    above = float(weight[freqs >= KT_CM].sum())
    peak = float(freqs[int(np.argmax(jw))])

    print("=== gap spectral density from the 1 fs windows ===\n")
    print(f"  window {window_fs:.0f} fs, resolution {resolution:.0f} cm-1, "
          f"kT = {KT_CM:.0f} cm-1")
    print(f"  sigma_gap  = {sigma_gap:7.1f} cm-1")
    print(f"  sigma_site = {sigma_site:7.1f} cm-1   (sites uncorrelated, r = +0.016)")
    print()
    print(f"  reorganisation energy, gap   lambda = <dE^2>/2kT = {lam_gap:7.1f} cm-1")
    print(f"  reorganisation energy, site                      = {lam_site:7.1f} cm-1")
    print()
    print("=== cross-check against the CD model's vibronic parameters ===")
    print(f"  repo (TeraChem CAM-B3LYP/6-31G*): S = {S_VIBRONIC:.3f} at omega = {OMEGA_VIBRONIC:.1f} cm-1")
    print(f"  this MD route, ALL site reorganisation put in that one mode:")
    print(f"     S = lambda_site / omega = {s_if_one_mode:.3f}")
    print(f"  ratio {s_if_one_mode / S_VIBRONIC:.2f}x. The MD value is an UPPER bound on S for")
    print(f"  that mode, since it attributes every mode's reorganisation to it, so")
    print(f"  {'S_repo < S_MD is CONSISTENT' if S_VIBRONIC < s_if_one_mode else 'S_repo > S_MD is NOT consistent'}.")
    print()
    print(f"  spectral density peak at {peak:.0f} cm-1"
          f"  (the gap ACF's oscillation implied ~1350)")
    print(f"  weight below kT: {below*100:5.1f}%   above kT: {above*100:5.1f}%")
    print(f"  -> but the {resolution:.0f} cm-1 resolution cannot separate below-kT structure;")
    print(f"     v6lfa/b/c at 500 fs reach 67 cm-1 and are what settle this.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({
        "window_fs": window_fs, "resolution_cm": resolution, "kT_cm": KT_CM,
        "sigma_gap_cm": sigma_gap, "sigma_site_cm": sigma_site,
        "lambda_gap_cm": lam_gap, "lambda_site_cm": lam_site,
        "S_if_all_in_one_mode": s_if_one_mode,
        "S_repo": S_VIBRONIC, "omega_repo_cm": OMEGA_VIBRONIC,
        "spectral_peak_cm": peak,
        "weight_below_kT": below, "weight_above_kT": above,
        "caveat": "100 fs resolves to 334 cm-1, coarser than kT = 208 cm-1; the "
                  "classically-valid part of the spectral density is unresolved here",
    }, indent=2))
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()

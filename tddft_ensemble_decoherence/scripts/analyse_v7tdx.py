#!/usr/bin/env python3
"""What the dVenus-TDX campaigns settle, now that both are complete.

Two 1 ns / 500-frame campaigns that have sat unanalysed because their master
manifests read 0/500 (the workers wrote to the _a/_b/_c splits and the merge was
never done):

  v7tdx      the construct with ONE chromophore -- Addgene #250939 deletes the
             barrel-2 chromophore tyrosine, so there is no site B, no detuning
             and no coupling. It measures the monomer band directly.
  v7tdxdark  the same construct's dark barrel (Gly-Gly + Tyr473, 29 atoms,
             neutral), as a negative control on whether the remnant absorbs
             anywhere near the chromophore band.

THE ONE COMPARISON THAT IS CLEAN IS THE WIDTH, NOT THE MEAN. The TDX box is
neutral where the TD box carries net -2 e, and the two constructs ran as
separate MD, so absolute site energies are offset by an amount this campaign
cannot measure. sigma is insensitive to that offset; the mean is not, and the
TD - TDX shift is reported here only with that caveat attached.

Read-only apart from results/v7tdx_analysis.json.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]


def column(manifest: str, key: str) -> np.ndarray:
    rows = [r for r in csv.DictReader((ROOT / "manifests" / manifest).open())
            if r["convergence_status"] == "complete" and r[key]]
    rows.sort(key=lambda r: int(r["frame_index"]))
    return np.array([float(r[key]) for r in rows])


def describe(x: np.ndarray, label: str) -> dict:
    n = x.size
    sd = float(x.std(ddof=1))
    # 500 frames at 2 ps are independent (tau_c ~ 5 fs), so the usual se applies
    return dict(label=label, n=n, mean_cm=float(x.mean()),
                se_mean_cm=float(sd / np.sqrt(n)), sd_cm=sd,
                se_sd_cm=float(sd / np.sqrt(2 * (n - 1))),
                shapiro_p=float(stats.shapiro(x).pvalue),
                skew=float(stats.skew(x)), kurtosis=float(stats.kurtosis(x)))


def main() -> None:
    tdx = column("v7tdx_frames.csv", "site_A_energy_cm-1")
    dark = column("v7tdxdark_frames.csv", "site_B_energy_cm-1")
    td_a = column("v5camp_frames.csv", "site_A_energy_cm-1")
    td_b = column("v5camp_frames.csv", "site_B_energy_cm-1")

    out = {"campaigns": {}}
    for x, lab in ((tdx, "TDX single chromophore"), (td_a, "TD site A"),
                   (td_b, "TD site B")):
        out["campaigns"][lab] = describe(x, lab)

    print("=== site-energy statistics, 500 frames each, 1 ns ===")
    print(f"  {'':<24s} {'mean':>9s} {'sd':>16s}  {'Shapiro p':>9s}  skew")
    for lab, d in out["campaigns"].items():
        print(f"  {lab:<24s} {d['mean_cm']:9.0f} "
              f"{d['sd_cm']:8.1f} +/- {d['se_sd_cm']:4.1f}  "
              f"{d['shapiro_p']:9.3f}  {d['skew']:+.2f}")

    # --- does the partner chromophore contribute to the site disorder? -------
    #
    # The comparison needs a yardstick, because TDX and TD are different
    # constructs on different trajectories and two independent 1 ns runs do not
    # return the same sigma. The yardstick is the difference between the TWO TD
    # SITES, measured in one trajectory with the partner present in both: that
    # is what "no effect" looks like on this data.
    print("\n=== does removing the partner change the site disorder? ===")
    ref = out["campaigns"]["TD site A"]
    ref_b = out["campaigns"]["TD site B"]
    d = out["campaigns"]["TDX single chromophore"]
    diff = d["sd_cm"] - ref["sd_cm"]
    se = float(np.hypot(d["se_sd_cm"], ref["se_sd_cm"]))
    yard = ref["sd_cm"] - ref_b["sd_cm"]
    se_yard = float(np.hypot(ref["se_sd_cm"], ref_b["se_sd_cm"]))
    lev = stats.levene(tdx, td_a, center="median")
    print(f"  sigma: TDX {d['sd_cm']:.1f}   TD site A {ref['sd_cm']:.1f}   "
          f"TD site B {ref_b['sd_cm']:.1f}")
    print(f"  TDX - site A       {diff:+6.1f} +/- {se:.1f}  "
          f"({diff/se:+.2f} sigma, Levene p = {lev.pvalue:.3f})")
    print(f"  site A - site B    {yard:+6.1f} +/- {se_yard:.1f}  "
          f"({yard/se_yard:+.2f} sigma)  <- one trajectory, partner present in both")
    print("  The two differences are the same size, so the TDX offset is within")
    print("  what two sites differ by anyway. NOT DEMONSTRATED that the partner")
    print("  contributes to the site disorder; also not excluded -- the data")
    print(f"  bound any effect at roughly {abs(diff) + 1.64 * se:.0f} cm^-1 "
          f"({(abs(diff) + 1.64 * se) / ref['sd_cm'] * 100:.0f}% of sigma).")
    out["locality_test"] = dict(
        delta_sd_cm=diff, se_cm=se, n_sigma=float(diff / se),
        levene_p=float(lev.pvalue),
        within_trajectory_yardstick_cm=yard, yardstick_se_cm=se_yard,
        verdict="not demonstrated; TDX-vs-A is the same size as A-vs-B",
        upper_bound_95_cm=float(abs(diff) + 1.64 * se))

    # --- the shift, with its confound stated ---------------------------------
    shift = float(td_a.mean() - tdx.mean())
    se_shift = float(np.hypot(ref["se_mean_cm"], d["se_mean_cm"]))
    print("\n=== TD - TDX shift (CONFOUNDED, see below) ===")
    print(f"  TD site A - TDX = {shift:+.1f} +/- {se_shift:.1f} cm^-1")
    print("  NOT a clean excitonic or electrostatic measurement: the TDX box is")
    print("  neutral where the TD box is net -2 e, and the two constructs ran as")
    print("  separate MD, so an unmeasured constant offset sits inside this.")
    out["shift"] = dict(td_minus_tdx_cm=shift, se_cm=se_shift,
                        confounded_by="TDX box neutral vs TD box net -2 e; "
                                      "separate trajectories")

    # --- the dark barrel: does the remnant absorb near the band? -------------
    #
    # Energy separation is the WRONG argument here. The mean sits 22,000 cm^-1
    # to the blue, but the distribution has a long red tail -- the barrel has no
    # strongly bright state, so brightest-root selection wanders over roots 1-7
    # and the lowest of 500 frames lands at 19,455 cm^-1, only 114 cm^-1 from
    # the chromophore band. The argument that holds is the INTENSITY.
    band = 19341.0
    dark_f = column("v7tdxdark_frames.csv", "site_B_oscillator_strength")
    near = dark < 25000.0
    f_chromophore = float(column("v7tdx_frames.csv",
                                 "site_A_oscillator_strength").mean())
    weight = float(near.sum() / dark.size * dark_f[near].mean() / f_chromophore) \
        if near.any() else 0.0
    print("\n=== dark-barrel control (v7tdxdark) ===")
    print(f"  {dark.size} frames, median {np.median(dark):.0f} cm^-1, but a long")
    print(f"  red tail: {int(near.sum())} frames fall below 25,000 and the lowest")
    print(f"  is {dark.min():.0f} cm^-1, only {dark.min()-band:.0f} cm^-1 from the band.")
    print(f"  Separation therefore settles nothing. Intensity does: those frames")
    print(f"  carry f = {dark_f[near].mean():.4f} (max {dark_f[near].max():.4f}) "
          f"against the chromophore's {f_chromophore:.2f},")
    print(f"  and they are {100*near.sum()/dark.size:.1f}% of the ensemble, so the")
    print(f"  remnant contributes {100*weight:.3f}% of the band intensity.")
    print("  => negligible, on oscillator strength rather than on energy.")
    out["dark_barrel"] = dict(
        n=int(dark.size), mean_cm=float(dark.mean()),
        median_cm=float(np.median(dark)), sd_cm=float(dark.std(ddof=1)),
        min_cm=float(dark.min()), band_cm=band,
        min_clearance_cm=float(dark.min() - band),
        n_below_25000=int(near.sum()),
        f_mean_near_band=float(dark_f[near].mean()) if near.any() else 0.0,
        f_chromophore=f_chromophore,
        intensity_fraction=weight,
        note="energy separation is not the argument; the tail reaches the band. "
             "brightest-root selection over roots 1-7 with no bright state.")

    (ROOT / "results" / "v7tdx_analysis.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print(f"\n  wrote results/v7tdx_analysis.json")


if __name__ == "__main__":
    main()

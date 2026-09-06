#!/usr/bin/env python3
"""The gap spectral density below v8dense's 33 cm^-1 floor, from the 8 ps record.

The quantum cumulant left one thing open. Its two residuals -- the static width
the band needs (161 cm^-1) and the under-predicted peak-to-peak Stokes shift
(94 against 402) -- both point below the lowest frequency a 1 ps window can
resolve. sigma cannot settle it: variance genuinely at 4 cm^-1 is reported at
33 by a 1 ps record, giving the SAME sigma but a very different lineshape,
because 33 cm^-1 still motionally narrows while 4 cm^-1 is static on the
absorption timescale.

This reads the stride-20 campaign (100 frames at 80 fs over 7.92 ps, resolution
4.2 cm^-1, Nyquist 208) and answers it directly.

The 33-208 cm^-1 OVERLAP with v8dense is the control. The two records come from
different trajectories (tandem_v4_8ps_4fs and tandem_v8_1ps_2fs), so their
spectral densities are only spliceable if they agree where both can see, and
that is checked here rather than assumed.

Read-only apart from results/lf_spectral_density.json.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v8dense_quantum_cumulant import (                                # noqa: E402
    power_spectrum, cumulant_g, lineshape, band_metrics, solve_static,
    emission_lineshape, first_moment,
    TWO_PI_C, KT_CM, EXP_PEAK_CM, EXP_FWHM_CM, EXP_ASYM, EXP_STOKES_CM)

LF_DT_FS = 80.0
OVERLAP = (33.0, 208.0)


def series(files, key):
    rows = [r for f in files for r in csv.DictReader((ROOT / "manifests" / f).open())
            if r["convergence_status"] == "complete" and r[key]]
    rows.sort(key=lambda r: int(r["frame_index"]))
    return np.array([float(r[key]) for r in rows])


def band_variance(omega, s, lo, hi):
    """Variance in [lo, hi] cm^-1, one-sided, in cm^-2."""
    nu = omega / TWO_PI_C
    w = s.copy()
    w[1:-1] *= 2.0
    m = (nu >= lo) & (nu < hi)
    if not m.any():
        return 0.0
    return float(np.trapezoid(w[m], omega[m]) / (2.0 * np.pi) / TWO_PI_C ** 2)


def main() -> None:
    out = {"overlap_cm": list(OVERLAP)}
    print("=== the band below 33 cm^-1, measured ===")
    print(f"  8 ps at {LF_DT_FS:.0f} fs -> resolution 4.2 cm^-1, Nyquist 208")
    print(f"  1 ps at 2 fs  -> resolution 33  cm^-1, Nyquist 8339\n")

    for site in ("A", "B"):
        key = f"site_{site}_energy_cm-1"
        lf = series([f"v4full_lf{w}_frames.csv" for w in "abc"], key)
        d8 = series([f"v8dense_{w}_frames.csv" for w in "abc"], key)

        om_lf, s_lf, var_lf = power_spectrum(lf, LF_DT_FS)
        om_d8, s_d8, var_d8 = power_spectrum(d8, 2.0)

        # --- the control: do they agree where both can see? -----------------
        v_lf = band_variance(om_lf, s_lf, *OVERLAP)
        v_d8 = band_variance(om_d8, s_d8, *OVERLAP)
        ratio = v_lf / v_d8 if v_d8 else float("nan")

        # --- the new information: 4.2 to 33 cm^-1 ---------------------------
        res = 1.0 / (lf.size * LF_DT_FS * 1e-15 * 2.99792458e10)
        v_new = band_variance(om_lf, s_lf, res, OVERLAP[0])
        sd_new = float(np.sqrt(max(v_new, 0.0)))
        sd_tot_lf = float(np.sqrt(var_lf) / TWO_PI_C)

        rec = dict(
            n_lf=int(lf.size), sd_lf_cm=sd_tot_lf,
            sd_dense_cm=float(np.sqrt(var_d8) / TWO_PI_C),
            overlap_sd_lf_cm=float(np.sqrt(v_lf)),
            overlap_sd_dense_cm=float(np.sqrt(v_d8)),
            overlap_variance_ratio=ratio,
            resolution_cm=res,
            sd_below_33_cm=sd_new,
            var_fraction_below_33=float(v_new / (sd_tot_lf ** 2)))
        out[f"site_{site}"] = rec

        print(f"  --- site {site} ---")
        print(f"    OVERLAP CONTROL {OVERLAP[0]:.0f}-{OVERLAP[1]:.0f} cm^-1:  "
              f"8 ps gives sd {np.sqrt(v_lf):.0f},  1 ps gives {np.sqrt(v_d8):.0f}"
              f"   variance ratio {ratio:.2f}")
        print(f"    NEW BAND {res:.1f}-{OVERLAP[0]:.0f} cm^-1: sd {sd_new:.0f} cm^-1"
              f"  = {100*rec['var_fraction_below_33']:.1f}% of the total variance")
        print(f"    (the cumulant needed a static 161 cm^-1 for site A, "
              f"276 for site B)\n")

    (ROOT / "results" / "lf_spectral_density.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print("  wrote results/lf_spectral_density.json")


if __name__ == "__main__":
    main()

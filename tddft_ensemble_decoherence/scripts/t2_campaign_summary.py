#!/usr/bin/env python3
"""Combine the three 1 fs T2 windows into one answer.

Each window is 100 frames at 1 fs, seeded from a v5 frame chosen to span the
detuning at near-constant coupling. The point of running three rather than one
contiguous 300 fs block is that T2 then gets a *configurational* error bar, and
the question "does T2 depend on the instantaneous detuning?" becomes measurable
instead of assumed.

The autocorrelation is built from the 100 new frames ONLY. The v5 seed frame is
a legitimate t = 0 for each window and attaching it would add a 101st point, but
it was computed in a different run, and lag 1 is the single most load-bearing
lag for a ~12 fs T2 -- a cross-run inconsistency there would propagate straight
into the answer. It is reported below as a sensitivity check, never as the
headline.

What the reused STEOM references do and do not buy here: they fix the absolute
site energies and the |mu| scale, and they let this campaign stay consistent
with the v5 calibration without 17 h of new CPU. They do NOT affect T2. The
cumulant sees only the centred gap, and a constant STEOM shift cancels exactly
in E_A - E_B -- summarize_100fs50.py asserts precisely this invariance.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from summarize_100fs50 import (
    acf_unbiased,
    cumulant,
    direct_coherence,
    first_crossing,
    moving_block_bootstrap_t2,
    positive_correlation_time,
)

ROOT = Path(__file__).resolve().parents[1]
EV_TO_CM = 8065.544005
HBAR_EV_FS = 0.6582119569

# window -> (v5 seed frame, |Delta| at that frame cm^-1, ensemble percentile)
WINDOWS = {
    "v6t2a": (900, 826.6, 71.2),
    "v6t2b": (1400, 58.6, 4.4),
    "v6t2c": (1900, 514.2, 51.4),
}


def guarded_crossing(time_fs: np.ndarray, coherence: np.ndarray,
                     floor: float) -> float | None:
    r"""1/e crossing, refusing to report one taken from the sampling noise floor.

    direct_coherence averages exp(-i \int Delta dt / hbar) over time origins, and
    at lag L only (n - L) origins remain. The magnitude of a mean of m random
    phases does not decay to zero; it fluctuates about ~1/sqrt(m). Once the curve
    reaches that level, everything after it is noise, and a 1/e crossing found
    there is an artefact.

    On the 4 fs data behind the manuscript's numbers this is not hypothetical.
    With n = 50 and tau_c ~ 12 fs there are ~17 independent samples, and the
    gap-only curve decays to 0.064 by 36 fs, then RECURS to 0.466 at 64 fs --
    crossing 1/e three times. Which crossing you report changes T2 by 5x. This
    function reports only a crossing reached before the curve first touches the
    floor, and returns None rather than a number it cannot defend.
    """
    target = float(np.exp(-1.0))
    # Only trust the curve up to the first time it touches the noise floor.
    touched = np.flatnonzero(coherence <= floor)
    horizon = int(touched[0]) if len(touched) else len(coherence)
    if horizon < 2:
        return None
    window = coherence[: horizon + 1]
    if window.min() > target:
        return None
    return first_crossing(time_fs[: horizon + 1], window)


def load(window: str) -> dict:
    path = ROOT / "manifests" / f"{window}_frames.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    complete = [r for r in rows if r["convergence_status"] == "complete"]
    failed = [r for r in rows if r["convergence_status"] == "failed"]
    if failed:
        raise RuntimeError(
            f"{window}: {len(failed)} failed frames. The cumulant needs the "
            f"contiguous measured series; gaps are not interpolated: "
            + ",".join(r["frame_index"] for r in failed)
        )
    if len(complete) != len(rows):
        raise RuntimeError(f"{window}: only {len(complete)}/{len(rows)} complete")
    t = np.asarray([float(r["time_fs"]) for r in complete])
    dt = float(np.median(np.diff(t)))
    if not np.allclose(np.diff(t), dt) or dt <= 0.0:
        raise RuntimeError(f"{window}: non-contiguous cadence")
    return {
        "dt_fs": dt,
        "ea": np.asarray([float(r["site_A_energy_eV"]) for r in complete]),
        "eb": np.asarray([float(r["site_B_energy_eV"]) for r in complete]),
        "jscr": np.asarray([float(r["J_screened_cm-1"]) for r in complete]),
    }


def analyse(window: str, data: dict) -> dict:
    dt = data["dt_fs"]
    gap = data["ea"] - data["eb"]
    result = cumulant(gap, dt)
    tau_c = positive_correlation_time(result["correlation_eV2"], dt)
    sigma_ev = float(np.std(gap, ddof=1))
    # Kubo parameter sigma*tau_c with sigma as an angular frequency (rad/fs).
    kubo = (sigma_ev / HBAR_EV_FS) * tau_c
    boot = moving_block_bootstrap_t2(gap, dt, max(3, int(np.ceil(2 * tau_c / dt))))
    direct_full = direct_coherence(data["ea"], data["eb"], data["jscr"], dt)
    direct_gap = direct_coherence(data["ea"], data["eb"], np.zeros_like(data["jscr"]), dt)
    time_fs = np.arange(len(gap), dtype=float) * dt
    n_indep = len(gap) * dt / tau_c if tau_c > 0 else float("nan")
    floor = 1.0 / np.sqrt(n_indep) if n_indep > 0 else 1.0

    # The stride test run UPWARD. On the old 4 fs data this could only be run
    # downward -- striding a series that is already under-sampled cannot detect
    # under-sampling, so its reassuring result carried no information. At 1 fs
    # it is a real convergence check: if T2* moves when the cadence is coarsened
    # to 2 and 4 fs, then the published 4 fs numbers were resolution-limited.
    stride = {}
    for s in (1, 2, 4):
        if len(gap) // s >= 8:
            stride[s] = cumulant(gap[::s], dt * s)["T2star_fs"]

    return {
        "n_frames": int(len(gap)),
        "dt_fs": dt,
        "gap_mean_cm": float(np.mean(gap) * EV_TO_CM),
        "gap_abs_mean_cm": float(np.mean(np.abs(gap)) * EV_TO_CM),
        "gap_sd_cm": float(sigma_ev * EV_TO_CM),
        "J_mean_cm": float(np.mean(data["jscr"])),
        "J_sd_cm": float(np.std(data["jscr"], ddof=1)),
        "T2star_cumulant_fs": result["T2star_fs"],
        "T2star_bootstrap": boot,
        "bath_correlation_time_fs": tau_c,
        "kubo_sigma_tau": kubo,
        "n_independent_samples": float(n_indep),
        "direct_noise_floor": float(floor),
        "T2_direct_full_signed_J_fs": guarded_crossing(time_fs, direct_full, floor),
        "T2_direct_gap_only_fs": guarded_crossing(time_fs, direct_gap, floor),
        "T2_direct_full_unguarded_fs": first_crossing(time_fs, direct_full),
        "T2_direct_gap_only_unguarded_fs": first_crossing(time_fs, direct_gap),
        "stride_sensitivity_fs": stride,
        "gap_acf_normalized": (
            result["correlation_eV2"] / result["correlation_eV2"][0]
        ).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/t2_campaign_summary.json")
    args = parser.parse_args()

    per_window = {}
    for window in WINDOWS:
        per_window[window] = analyse(window, load(window))

    print("=== T2 campaign: three 100 fs windows at 1 fs ===")
    print()
    header = (f"{'window':8s} {'seed':>5s} {'|D| seed':>9s} {'|D| window':>11s} "
              f"{'sd(D)':>8s} {'J':>7s} {'T2* cum':>9s} {'tau_c':>7s} {'Kubo':>6s}")
    print(header)
    print("-" * len(header))
    for window, (seed, det, pct) in WINDOWS.items():
        r = per_window[window]
        t2 = r["T2star_cumulant_fs"]
        print(f"{window:8s} {seed:5d} {det:9.1f} {r['gap_abs_mean_cm']:11.1f} "
              f"{r['gap_sd_cm']:8.1f} {r['J_mean_cm']:7.2f} "
              f"{(f'{t2:.2f}' if t2 is not None else 'none'):>9s} "
              f"{r['bath_correlation_time_fs']:7.2f} {r['kubo_sigma_tau']:6.2f}")

    t2s = [per_window[w]["T2star_cumulant_fs"] for w in WINDOWS]
    if all(v is not None for v in t2s):
        t2a = np.asarray(t2s, dtype=float)
        dets = np.asarray([per_window[w]["gap_abs_mean_cm"] for w in WINDOWS])
        print()
        print(f"  T2* across windows: mean {t2a.mean():.2f} fs, "
              f"sd {t2a.std(ddof=1):.2f} fs, range {t2a.min():.2f}-{t2a.max():.2f}")
        # Three points: a correlation coefficient here is descriptive only.
        r = float(np.corrcoef(dets, t2a)[0, 1])
        print(f"  T2* vs window <|Delta|>: r = {r:+.3f} over 3 points "
              f"(descriptive only, not a test)")
        fwhm = [1.0 / (np.pi * 2.99792458e-5 * v) for v in t2a]
        print(f"  homogeneous FWHM 1/(pi c T2*): "
              f"{min(fwhm):.0f}-{max(fwhm):.0f} cm^-1")

    print()
    print("  direct two-site propagation (e^-1 crossing, noise-floor guarded):")
    for window in WINDOWS:
        r = per_window[window]
        print(f"    {window}: with signed J {_f(r['T2_direct_full_signed_J_fs']):>8s} fs   "
              f"gap only {_f(r['T2_direct_gap_only_fs']):>8s} fs   "
              f"[{r['n_independent_samples']:.0f} indep samples, "
              f"floor {r['direct_noise_floor']:.2f}]")
        print(f"      unguarded (for comparison only): "
              f"{_f(r['T2_direct_full_unguarded_fs']):>8s} / "
              f"{_f(r['T2_direct_gap_only_unguarded_fs']):>8s} fs")
    print("    The cumulant is the estimator to quote. It integrates C(tau) only")
    print("    out to t ~ T2, where the ACF is well sampled; the direct curve is")
    print("    read at long lag, where it is not. On the old 4 fs data the")
    print("    gap-only direct curve recurs from 0.064 at 36 fs to 0.466 at 64 fs,")
    print("    crossing 1/e three times -- which is where a 60 fs number can come")
    print("    from without any physics behind it.")

    print()
    print("  stride test, run upward from 1 fs (T2* in fs):")
    for window in WINDOWS:
        s = per_window[window]["stride_sensitivity_fs"]
        cells = "  ".join(f"{k} fs: {_f(v):>7s}" for k, v in sorted(s.items()))
        print(f"    {window}: {cells}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({
        "windows": {w: {**per_window[w], "seed_frame": WINDOWS[w][0],
                        "seed_detuning_cm": WINDOWS[w][1],
                        "seed_ensemble_percentile": WINDOWS[w][2]}
                    for w in WINDOWS},
        "provenance": {
            "construct": "511-aa Addgene 250938 with AF3 C-terminal tails",
            "embedding": "full-system point charges, link-only exclusion",
            "profile": "def2-SVP, 7-root density-overlap tracking",
            "steom_references": "reused from the v5 sweep; constant shift, "
                                "cancels in the gap and does not affect T2",
            "supersedes": "rapid_100fs50_fixedmm (495-aa construct, 12 A cutoff, 4 fs)",
        },
    }, indent=2))
    print()
    print(f"wrote {args.output}")


def _f(v) -> str:
    return "none" if v is None else f"{v:.2f}"


if __name__ == "__main__":
    main()

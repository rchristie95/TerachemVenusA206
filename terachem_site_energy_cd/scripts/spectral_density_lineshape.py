#!/usr/bin/env python3
"""Lineshape from the full gap spectral density, not one effective HWHM.

WHY. The couplet-width result (342.7 computed vs 366.4 observed) currently
depends on a fast/slow split placed at the 592 fs window length: everything
faster went into the cumulant, everything slower was added as static disorder.
That cut is defensible but arbitrary, and it is the one free choice left in a
result advertised as parameter-free.

The cumulant is exact for Gaussian fluctuations at ANY timescale, so the honest
object is the measured correlation function itself plus an explicit statement of
what the window cannot see:

    resolved by the 592 fs window : omega > ~56 cm^-1
    sampled by the 1 ns ensemble  : omega < ~7 cm^-1 (effectively static)
    UNCONSTRAINED                 : 7 - 56 cm^-1

This script computes C(t) and its spectral density from the 4 fs campaign
(site A, bright-moment estimator), builds the cumulant lineshape from it
directly, and then reports how the couplet width moves as the unconstrained
low-frequency variance is swept. Sensitivity is stated rather than chosen.

Read-only apart from results/spectral_density_lineshape.json.
"""
from __future__ import annotations
import json, glob, os, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cumulant_lineshape import cumulant_g, lineshape_from_g, fwhm_of   # noqa: E402
sys.path.insert(0, str(ROOT.parent / "tddft_ensemble_decoherence/scripts"))
from run_ensemble import _bright_moment                                # noqa: E402

EV_TO_CM = 8065.544
DT_FS = 4.0
TWO_PI_C = 2.0 * np.pi * 2.99792458e-5     # cm^-1 -> rad/fs
BASE = ROOT.parent / "tddft_ensemble_decoherence/frames/short"


def site_series(site: str) -> np.ndarray:
    frames = sorted(int(p.name.split("_")[-1]) for p in BASE.glob("frame_*"))
    out = []
    for f in frames:
        p = BASE / f"frame_{f:04d}" / f"site_{site}" / "density_summary.json"
        if not p.exists():
            continue
        roots = json.loads(p.read_text())["roots"]
        out.append(_bright_moment({r["root"]: r for r in roots}) * EV_TO_CM)
    return np.array(out)


def main() -> None:
    a, b = site_series("A"), site_series("B")
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    print(f"4 fs window: {n} frames = {(n-1)*DT_FS:.0f} fs")
    for lab, x in (("site A", a), ("site B", b), ("gap A-B", a - b)):
        print(f"   {lab:9s} sd {x.std():7.1f} cm^-1")

    # correlation function of the SITE energy (site A: the clean one)
    x = a - a.mean()
    nl = n // 2
    C = np.array([np.dot(x[:n - l], x[l:]) / (n - l) for l in range(nl)])
    lag = np.arange(nl) * DT_FS

    # spectral density: cosine transform of C(t)
    w = np.linspace(5.0, 2000.0, 800)                 # cm^-1
    Jw = np.array([np.trapezoid(C * np.cos(TWO_PI_C * wi * lag), lag) for wi in w])
    lo, hi = 1.0 / ((n - 1) * DT_FS * 2.99792458e-5), 1.0 / (2 * DT_FS * 2.99792458e-5)
    print(f"\nwindow resolves {lo:.0f} - {hi:.0f} cm^-1 "
          f"(Nyquist {hi:.0f}); below {lo:.0f} cm^-1 is unconstrained")
    peaks = w[1:-1][(Jw[1:-1] > Jw[:-2]) & (Jw[1:-1] > Jw[2:]) & (Jw[1:-1] > 0.15 * Jw.max())]
    print(f"spectral density peaks at: {', '.join(f'{p:.0f}' for p in peaks[:8])} cm^-1")

    grid = np.arange(17500.0, 23500.0, 2.0)
    centre = 19322.0
    g = cumulant_g(lag, C)
    base = fwhm_of(grid, lineshape_from_g(lag, g, grid, centre)) / 2.0
    print(f"\ncumulant HWHM from the measured C(t) alone: {base:.1f} cm^-1")

    # sweep the UNCONSTRAINED low-frequency variance as extra static disorder
    print(f"\n{'extra static sd':>16} {'total site sd':>14} {'HWHM':>8} {'couplet-equivalent':>19}")
    rows = []
    for extra in (0.0, 100.0, 150.0, 200.0, 250.0, 300.0):
        tot = float(np.sqrt(a.var() + extra**2))
        # a static Gaussian of sd `extra` convolves the homogeneous lineshape
        y = lineshape_from_g(lag, g, grid, centre)
        if extra > 0:
            s = extra
            k = np.exp(-0.5 * ((grid - centre) / s) ** 2); k /= k.sum()
            y = np.convolve(y, k, mode="same")
        hw = fwhm_of(grid, y) / 2.0
        rows.append(dict(extra_static_sd=extra, total_site_sd=tot, hwhm=hw))
        print(f"{extra:>16.0f} {tot:>14.1f} {hw:>8.1f} {hw:>19.1f}")
    print(f"\nobserved couplet widths (Nguyen Table S3): 293.5 and 439.3, mean 366.4")
    out = ROOT / "results" / "spectral_density_lineshape.json"
    out.write_text(json.dumps(dict(
        n_frames=int(n), window_fs=float((n - 1) * DT_FS),
        site_A_sd=float(a.std()), site_B_sd=float(b.std()), gap_sd=float((a - b).std()),
        resolved_cm=[float(lo), float(hi)],
        spectral_peaks_cm=[float(p) for p in peaks[:8]],
        hwhm_from_Ct=float(base), sweep=rows), indent=2) + "\n")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Item 2: confront the Holstein vibronic exciton CD with Nguyen 2025 Table S3.

WHY. The purely electronic two-state model fails two independent ways
(scripts/nguyen_fit_on_simulation.py):
  * run through THEIR constrained fit our spectrum returns delta = 275.9 cm^-1
    against their 130.79, and no (Delta, linewidth) combination gets below
    ~276 -- delta rises as the detuning falls;
  * our couplet is near-conservative (lobe ratio 1.06) where theirs is strongly
    asymmetric (1.982e5 / 0.847e5 = 2.34), and a conservative couplet is a
    theorem of the two-state model, not an accident of parameters.
Both point at intensity the electronic model does not contain. Nguyen's own
third component sits ~1452 cm^-1 to the blue, carries ~18% of the couplet area
with the SAME sign as the negative lobe, and is never assigned -- a C=C/C=N
vibronic sideband on the natural reading.

WHAT. Build the ensemble CD from the same 95-frame QM/MM ensemble using the
one-particle Holstein model in vibronic_exciton.py, scan (S, omega_vib), and
score each against four independent Table S3 observables:
    delta (their constrained fit re-run on our spectrum), the two fitted
    widths, the lobe asymmetry, and the third-band offset.

The fit is scale-normalised and multi-start: the constrained 3-Lorentzian has
multiple local minima and curve_fit lands in different ones depending on the
data's absolute scale (an unnormalised first pass returned 228 where the global
optimum is 276).

Read-only apart from results/vibronic_cd_vs_experiment.json. No network.
"""
from pathlib import Path
import json, sys
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from vibronic_exciton import vibronic_spectra          # noqa: E402

ENSEMBLE = ROOT / "results" / "ensembles" / "ens_v2_all.npz"
CENTRE = 19322.0
J_CM = 32.8165
HWHM = 0.5 * (293.51 + 439.25)
GRID = np.arange(17500.0, 23500.0, 2.0)
FITW = (18000.0, 21800.0)
EXP = dict(delta=130.79, w1=293.51, w2=439.25, asym=1.982e5 / 0.8466e5,
           third_offset=20774.0 - CENTRE)


def model(nu, a1, w1, a2, w2, delta, a3, nu3, w3):
    return (a1 / (1 + ((nu - (CENTRE - delta)) / w1) ** 2)
            + a2 / (1 + ((nu - (CENTRE + delta)) / w2) ** 2)
            + a3 / (1 + ((nu - nu3) / w3) ** 2))


def robust_fit(nu, y):
    y = y / np.abs(y).max()
    best = None
    for d0 in (60., 130., 200., 275., 350., 450., 600.):
        for s1 in (1., -1.):
            try:
                p, _ = curve_fit(model, nu, y,
                                 p0=[s1, 293.51, -0.6 * s1, 439.25, d0, -0.1 * s1, 20774., 565.79],
                                 bounds=([-10, 50, -10, 50, 5, -10, 20000, 50],
                                         [10, 2000, 10, 2000, 2000, 10, 21600, 2000]),
                                 maxfev=400000)
            except Exception:
                continue
            r = y - model(nu, *p)
            r2 = 1 - np.sum(r ** 2) / np.sum((y - y.mean()) ** 2)
            if best is None or r2 > best[-1]:
                best = (p[4], p[1], p[3], p[5], p[6], r2)
    return best


def main():
    d = np.load(ENSEMBLE)
    ea, eb = d["e_a_cm"], d["e_b_cm"]
    T = np.einsum("ij,ij->i", d["r_b_ang"] - d["r_a_ang"],
                  np.cross(d["mu_a_au"], d["mu_b_au"]))
    shift = CENTRE - 0.5 * (ea + eb).mean()
    frames = list(zip(ea + shift, eb + shift, np.full(len(ea), J_CM), T))
    m = (GRID >= FITW[0]) & (GRID <= FITW[1])

    print(f"ensemble n={len(frames)}  <|D|>={np.abs(ea-eb).mean():.1f}  J={J_CM:.2f}")
    print(f"targets: delta {EXP['delta']:.2f} | widths {EXP['w1']:.0f}/{EXP['w2']:.0f}"
          f" | lobe asym {EXP['asym']:.2f} | third band +{EXP['third_offset']:.0f} cm^-1\n")
    print(f"{'S':>5} {'w_vib':>6} {'delta':>8} {'w1':>7} {'w2':>7} {'asym':>6} {'3rd off':>8} {'R2':>8}")

    rows = []
    for S in (0.0, 0.2, 0.35, 0.5, 0.7, 1.0):
        for wv in (1450.0,):
            _, cd = vibronic_spectra(GRID, frames, HWHM, omega_cm=wv,
                                     huang_rhys=S, n_max=5)
            y = cd[m]
            asym = abs(y.min()) / y.max() if y.max() > 0 else np.nan
            dl, w1, w2, a3, nu3, r2 = robust_fit(GRID[m], y)
            rows.append(dict(S=S, omega=wv, delta=dl, w1=w1, w2=w2,
                             asym=float(asym), third_offset=nu3 - CENTRE, r2=r2))
            print(f"{S:5.2f} {wv:6.0f} {dl:8.2f} {w1:7.1f} {w2:7.1f} {asym:6.2f}"
                  f" {nu3-CENTRE:8.0f} {r2:8.5f}")

    out = ROOT / "results" / "vibronic_cd_vs_experiment.json"
    out.write_text(json.dumps({"experiment": EXP, "scan": rows}, indent=2))
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Test B: apply Nguyen's own constrained 3-Lorentzian fit to OUR simulated CD.

The dispute over the dVenus couplet is an INVERSION dispute: their fit returns an
apparent half-splitting of 130.79 cm^-1, our ensemble has <|Delta|> = 576 and
J = 32.8. Comparing 130.79 to J, or to <Omega> = 584, or to the CD-weighted
harmonic <Omega>_CD = 304, all compare different quantities -- and in the
overlapped regime (HWHM 366 > splitting) the derivation note says none of those
mappings holds anyway.

The like-for-like test is to run THEIR reduction on OUR spectrum and read off the
same number they read off. That is what this does.

Protocol, from Nguyen 2025 Table S3: three peak-normalised Lorentzians,
    f(nu) = sum_k A_k / (1 + ((nu - nu_k)/w_k)^2),
with the two couplet components constrained SYMMETRIC about a centre fixed at
the TDX absorption maximum (19322 cm^-1) and separated by 2*delta. Free
parameters: A1, w1, A2, w2, delta, and the third band (A3, nu3, w3).
Their published solution: A=(+3.59,-2.27,-0.34)e5, delta=130.79,
w=(293.51,439.25,565.79), nu3=20774.

CONTROL: the identical fit is first run on the experimental spectrum. If it
recovers delta = 130.79 the protocol is faithful, and whatever it returns on the
simulated spectrum is an apples-to-apples answer.

Read-only apart from results/nguyen_fit_on_simulation.json. No network.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent
SPEC = ROOT / "results" / "absolute_cd_spectrum.npz"
CENTRE_CM = 19322.0          # TDX absorption maximum -- FIXED, as in their fit
PUBLISHED_DELTA = 130.79
FIT_WINDOW = (18000.0, 21800.0)


def model(nu, a1, w1, a2, w2, delta, a3, nu3, w3):
    """Constrained 3-Lorentzian: couplet symmetric about the fixed centre."""
    l1 = a1 / (1.0 + ((nu - (CENTRE_CM - delta)) / w1) ** 2)
    l2 = a2 / (1.0 + ((nu - (CENTRE_CM + delta)) / w2) ** 2)
    l3 = a3 / (1.0 + ((nu - nu3) / w3) ** 2)
    return l1 + l2 + l3


def fit(nu, y, label):
    scale = np.abs(y).max()
    p0 = [scale, 293.51, -0.6 * scale, 439.25, PUBLISHED_DELTA, -0.1 * scale, 20774.0, 565.79]
    lo = [-10 * scale, 50.0, -10 * scale, 50.0, 5.0, -10 * scale, 20000.0, 50.0]
    hi = [10 * scale, 2000.0, 10 * scale, 2000.0, 2000.0, 10 * scale, 21600.0, 2000.0]
    p, _ = curve_fit(model, nu, y, p0=p0, bounds=(lo, hi), maxfev=200000)
    resid = y - model(nu, *p)
    r2 = 1.0 - np.sum(resid**2) / np.sum((y - y.mean()) ** 2)
    print(f"  {label}")
    print(f"    delta (apparent half-splitting) = {p[4]:8.2f} cm^-1"
          f"   -> apparent 2|J| = {2*p[4]:.1f}")
    print(f"    A1 {p[0]:+.3e} w1 {p[1]:6.1f} | A2 {p[2]:+.3e} w2 {p[3]:6.1f}"
          f" | A3 {p[5]:+.3e} nu3 {p[6]:.0f} w3 {p[7]:.0f}")
    print(f"    R^2 = {r2:.5f}")
    return p, r2


def main():
    d = np.load(SPEC)
    nu, pred, exp = d["grid_cm"], d["predicted_theta"], d["experimental_theta"]
    m = (nu >= FIT_WINDOW[0]) & (nu <= FIT_WINDOW[1])
    nu, pred, exp = nu[m], pred[m], exp[m]

    print(f"Nguyen constrained 3-Lorentzian fit, centre FIXED at {CENTRE_CM:.0f} cm^-1")
    print(f"window {FIT_WINDOW[0]:.0f}-{FIT_WINDOW[1]:.0f} cm^-1, {len(nu)} points\n")
    print("CONTROL -- their fit re-run on the experimental spectrum")
    pe, r2e = fit(nu, exp, "experimental")
    ok = abs(pe[4] - PUBLISHED_DELTA) / PUBLISHED_DELTA < 0.10
    print(f"    published delta = {PUBLISHED_DELTA:.2f}  -> protocol "
          f"{'REPRODUCES it' if ok else 'DOES NOT reproduce it'}"
          f" ({100*(pe[4]-PUBLISHED_DELTA)/PUBLISHED_DELTA:+.1f}%)\n")

    print("TEST -- the same fit applied to our simulated spectrum")
    pp, r2p = fit(nu, pred, "simulated (J = 32.82, <|Delta|> = 576)")

    print(f"\n{'='*66}")
    print(f"  their reduction on their data : delta = {pe[4]:7.2f} cm^-1")
    print(f"  their reduction on our  data : delta = {pp[4]:7.2f} cm^-1")
    print(f"  ratio (ours / theirs)        = {pp[4]/pe[4]:.2f}")
    print(f"  our actual J                 = 32.82 cm^-1")
    print(f"  inflation their fit imposes on our ensemble: {pp[4]/32.82:.2f}x")
    print(f"{'='*66}")

    out = ROOT / "results" / "nguyen_fit_on_simulation.json"
    out.write_text(json.dumps({
        "centre_fixed_cm": CENTRE_CM, "window_cm": list(FIT_WINDOW),
        "published_delta_cm": PUBLISHED_DELTA,
        "control_experimental": {"delta_cm": pe[4], "r2": r2e,
                                 "params": list(map(float, pe))},
        "test_simulated": {"delta_cm": pp[4], "r2": r2p,
                           "params": list(map(float, pp))},
        "ratio_sim_over_exp": pp[4] / pe[4],
        "j_true_cm": 32.8165, "inflation_over_true_J": pp[4] / 32.8165,
    }, indent=2))
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

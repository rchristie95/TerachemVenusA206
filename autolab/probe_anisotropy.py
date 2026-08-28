#!/usr/bin/env python3
"""How much of the anisotropy constraint survives its own stated assumption?

joint_constraints.py inverts the R0 drop with

    R_TD = (R_TDX/2) * [1 + (3 cos^2 alpha - 1)/2]

whose factor of 1/2 is an assumption, flagged in that script's own docstring:
complete sub-IRF transfer, with the excitation ending on either chromophore
with equal probability. The general form is

    R_TD = R_TDX * [ (1 - f) + f * d2 ],    d2 = (3 cos^2 alpha - 1)/2

with f the fraction of excitations that end up on the *partner*. f = 1/2
reproduces the published formula exactly (checked below, and the script asserts
it).

That assumption is worth probing because the dimer is detuning-dominated
(|Delta| ~ 253-1310 >> 2|J| = 65.6), so transfer is one-way downhill, not a
symmetric shuffle. After full thermal equilibration between two sites split by
Delta, with w the probability that the *upper* site is the one initially
excited:

    f(Delta, w) = w * p_low + (1 - w) * p_up,   p from Boltzmann at 300 K

Note w = 1/2 gives f = 1/2 for ANY Delta -- the published inversion is exactly
robust to detuning as long as the two sites are excited with equal probability.
So the free parameter that can actually move |cos alpha| is w, i.e. unequal
two-photon cross-sections for the two chromophores, not the detuning.

This scans w and reports the implied |cos alpha| and alpha. It answers: how
badly does w have to be skewed before the anisotropy stops demanding 131.3 deg
and starts agreeing with the crystal register at 110.4 deg?
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "exciton_observables"))
import joint_constraints as jc  # noqa: E402

DELTAS_CM = {"CD splitting": 253.2, "computed QM/MM": 576.4, "third handle": 1310.0}
CRYSTAL_ALPHA, TARGET_ALPHA = 110.44, 131.3


def transferred_fraction(delta_cm: float, w_upper: float) -> float:
    """Equilibrium fraction on the partner site, averaged over which site is excited."""
    boltz = np.exp(-delta_cm / jc.KT)
    p_low, p_up = 1.0 / (1.0 + boltz), boltz / (1.0 + boltz)
    return w_upper * p_low + (1.0 - w_upper) * p_up


def cos_alpha_from(f: float, r_td: float = jc.R_TD, r_tdx: float = jc.R_TDX) -> float:
    """Invert R_TD = R_TDX [(1-f) + f d2] for |cos alpha|; nan where unphysical."""
    if f <= 0.0:
        return float("nan")
    d2 = (r_td / r_tdx - (1.0 - f)) / f
    c2 = (2.0 * d2 + 1.0) / 3.0
    return float(np.sqrt(c2)) if 0.0 <= c2 <= 1.0 else float("nan")


def main() -> None:
    # The generalisation must reproduce the published number at f = 1/2.
    published, _ = jc.anisotropy_cos()
    assert abs(cos_alpha_from(0.5) - published) < 1e-9, "generalisation broke the f=1/2 limit"
    print(f"sanity: f = 1/2 reproduces joint_constraints |cos alpha| = {published:.4f}  OK\n")

    print("f is independent of Delta at w = 1/2 -- detuning alone cannot rescue the fit:")
    for name, d in DELTAS_CM.items():
        f = transferred_fraction(d, 0.5)
        print(f"  {name:<16} Delta = {d:>7.1f} cm^-1   f = {f:.6f}   "
              f"|cos a| = {cos_alpha_from(f):.4f}")

    print("\nScan over w, the probability the UPPER site is the one excited")
    print("(Delta = 576.4 cm^-1; at 300 K the downhill site takes essentially everything):\n")
    print(f"{'w_upper':>9}{'f':>10}{'|cos a|':>10}{'alpha_deg':>11}   note")
    rows = []
    target_hits = {}
    for w in np.linspace(0.0, 1.0, 21):
        f = transferred_fraction(DELTAS_CM["computed QM/MM"], float(w))
        c = cos_alpha_from(f)
        alpha = float(np.degrees(np.arccos(-c))) if np.isfinite(c) else float("nan")
        note = ""
        if np.isfinite(alpha):
            if abs(alpha - CRYSTAL_ALPHA) < 1.5:
                note = "<-- agrees with the crystal register"
                target_hits["crystal"] = float(w)
            elif abs(alpha - TARGET_ALPHA) < 1.5:
                note = "<-- the published 131.3 deg"
                target_hits["published"] = float(w)
        elif f > 0:
            note = "unphysical (no real alpha)"
        rows.append({"w_upper": float(w), "f": float(f),
                     "abs_cos_alpha": c, "alpha_deg": alpha})
        print(f"{w:>9.2f}{f:>10.4f}{c:>10.4f}{alpha:>11.1f}   {note}")

    finite = [r for r in rows if np.isfinite(r["alpha_deg"])]
    reach = any(abs(r["alpha_deg"] - CRYSTAL_ALPHA) < 1.5 for r in finite)
    print(f"\ncrystal register reachable by skewing w alone: {reach}")
    if not reach:
        span = (min(r["alpha_deg"] for r in finite), max(r["alpha_deg"] for r in finite))
        print(f"  alpha spans only {span[0]:.1f}-{span[1]:.1f} deg over all w in [0,1];")
        print("  the anisotropy constraint is NOT rescued by the transfer model, so the")
        print("  disagreement is geometric or lies in R0 itself.")

    out = {"published_cos_alpha_at_f_half": float(published),
           "delta_independence_at_w_half": {k: transferred_fraction(v, 0.5)
                                            for k, v in DELTAS_CM.items()},
           "scan_delta_cm": DELTAS_CM["computed QM/MM"],
           "scan": rows,
           "crystal_reachable": bool(reach),
           "w_at_targets": target_hits}
    dest = Path(__file__).resolve().parent / "results" / "probe_anisotropy.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nWrote {dest}")


if __name__ == "__main__":
    main()

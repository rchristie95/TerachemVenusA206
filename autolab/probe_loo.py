#!/usr/bin/env python3
"""Leave-one-out on the four-observable joint fit: which observable is the liar?

exciton_observables/joint_constraints.py gets chi^2 = 45.95 for 2 dof, so the
two-state model plus the four measurements are mutually inconsistent. That is a
falsification, but it does not say *what* is wrong. Dropping each observable in
turn and refitting the remaining three (1 dof) does: if one removal collapses
chi^2 to O(1) and the others do not, that observable carries the inconsistency
on its own and is the thing to attack.

All constants, formulae and grids are imported from joint_constraints.py -- no
physics is restated here, so this cannot drift from the published fit.

Costs seconds. Run it before spending a GPU-hour on anything geometric.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "exciton_observables"))
import joint_constraints as jc  # noqa: E402


def chi_surfaces():
    """Rebuild the four chi^2 surfaces exactly as joint_constraints.main() does."""
    cos_exp, cos_exp_sd = jc.anisotropy_cos()
    shift_cos = jc.SHIFT / jc.J0
    shift_cos_sd = abs(shift_cos) * np.hypot(jc.SHIFT_SD / jc.SHIFT, jc.J_SD / jc.J0)
    target = jc.DTAU / jc.TAU_PS
    target_sd = target * np.hypot(jc.DTAU_SD / jc.DTAU, jc.PHI_SD / jc.PHI)

    cos_grid = np.linspace(0.0, 1.0, 1201)
    del_grid = np.linspace(0.0, 2500.0, 2501)
    C, D = np.meshgrid(cos_grid, del_grid, indexing="ij")

    chi = {
        "superradiance": ((jc.superradiance(C, D) - target) / target_sd) ** 2,
        "anisotropy": ((C - cos_exp) / cos_exp_sd) ** 2 * np.ones_like(D),
        "CD splitting": ((D - jc.CD_DELTA) / jc.CD_DELTA_SD) ** 2 * np.ones_like(C),
        "red shift": ((C - shift_cos) / shift_cos_sd) ** 2 * np.ones_like(D),
    }
    return chi, cos_grid, del_grid


def best(total, cos_grid, del_grid):
    a, b = np.unravel_index(np.argmin(total), total.shape)
    return float(cos_grid[a]), float(del_grid[b]), float(total[a, b]), (a, b)


def main() -> None:
    chi, cos_grid, del_grid = chi_surfaces()
    names = list(chi)

    c_all, d_all, x_all, _ = best(sum(chi.values()), cos_grid, del_grid)
    print(f"ALL FOUR       |cos a| = {c_all:.3f}  |Delta| = {d_all:.0f}  "
          f"chi2 = {x_all:.2f} / 2 dof\n")

    print(f"{'dropped':<16}{'|cos a|':>9}{'|Delta|':>10}{'chi2':>9}{'alpha_deg':>11}  verdict")
    rows = {}
    for n in names:
        total = sum(v for k, v in chi.items() if k != n)
        c, d, x, _ = best(total, cos_grid, del_grid)
        alpha = float(np.degrees(np.arccos(-c)))  # superradiance forces cos alpha < 0
        rows[n] = {"abs_cos_alpha": c, "abs_delta_cm": d, "chi2": x,
                   "dof": 1, "alpha_deg": alpha}
        tag = "EXONERATES the rest" if x < 4.0 else ("helps" if x < x_all / 3 else "")
        print(f"{n:<16}{c:>9.3f}{d:>10.0f}{x:>9.2f}{alpha:>11.1f}  {tag}")

    worst = min(rows, key=lambda n: rows[n]["chi2"])
    drop = x_all - rows[worst]["chi2"]
    exonerating = sorted((n for n in names if rows[n]["chi2"] < 4.0),
                         key=lambda n: rows[n]["chi2"])
    print(f"\nlargest single culprit: '{worst}' "
          f"(chi2 {x_all:.1f} -> {rows[worst]['chi2']:.1f}, drop {drop:.1f})")

    if len(exonerating) > 1:
        # A tie is not an ambiguity to be broken -- it localises the tension.
        alphas = [rows[n]["alpha_deg"] for n in exonerating]
        deltas = [rows[n]["abs_delta_cm"] for n in exonerating]
        print(f"  -> {len(exonerating)} removals each give a consistent fit: "
              f"{', '.join(exonerating)}")
        print(f"     they agree on alpha ({min(alphas):.1f}-{max(alphas):.1f} deg) and "
              f"disagree on |Delta| ({min(deltas):.0f} vs {max(deltas):.0f} cm^-1).")
        print("     The geometry is NOT the problem. Two observables demand "
              "incompatible detunings; that pair is the target.")
    elif exonerating:
        print("  -> the remaining three ARE mutually consistent. Attack this observable.")
    else:
        print("  -> no single observable explains the tension; the two-state MODEL is "
              "the suspect, not any one measurement.")

    out = {"all_four": {"abs_cos_alpha": c_all, "abs_delta_cm": d_all,
                        "chi2": x_all, "dof": 2},
           "leave_one_out": rows,
           "largest_culprit": worst,
           "exonerating_removals": exonerating,
           "chi2_drop": drop,
           "single_observable_explains": bool(rows[worst]["chi2"] < 4.0)}
    dest = Path(__file__).resolve().parent / "results" / "probe_loo.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nWrote {dest}")


if __name__ == "__main__":
    main()

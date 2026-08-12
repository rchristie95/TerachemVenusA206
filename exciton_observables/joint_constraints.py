#!/usr/bin/env python3
"""Are the experimental observables mutually consistent in a two-state model?

Four measurements constrain only two unknowns, (|cos alpha|, |Delta|), with J
fixed at 32.82 +/- 1.55 cm^-1 by the transition-density calculation. The system
is over-determined, which means it can be falsified -- and nobody has tested it.

Observables, each reduced to a constraint on (|cos alpha|, |Delta|):

  1. SUPERRADIANCE. dVenus-TD decays 57 +/- 4 ps faster than TDX on 3.026 ns.
       Phi * tanh(Omega/2kT) * (2J/Omega) * |cos alpha| = Dtau/tau
     Couples both unknowns. Omega = sqrt(Delta^2 + 4J^2).

  2. LIMITING ANISOTROPY. R0 falls 0.52 (TDX) -> 0.30 (TD) under two-photon
     excitation. With complete sub-IRF transfer -- verified: the Marcus rate
     from our own J and Delta gives 0.9-1.7 ps against a 26.5 ps IRF --
       R_TD = (R_TDX/2) * [1 + (3 cos^2 alpha - 1)/2]
     Constrains |cos alpha| alone.

  3. CD COUPLET SPLITTING. Nguyen's constrained three-Lorentzian fit gives an
     apparent splitting of 261.58 cm^-1. The apparent-vs-true inversion in
     terachem_site_energy_cd/SPLITTING_PROFILE_INVESTIGATION.md maps that onto
       <|Delta|> = 253.2 cm^-1   (the unconstrained 370.78 fit gives 364.9)
     Constrains |Delta| alone.

  4. ABSORPTION RED SHIFT. The dipole-strength-weighted first moment shifts by
     exactly J cos(alpha), independent of Delta, so it survives ensemble
     averaging. Observed TDX->TD shift 35.3 cm^-1, of which 15.6 +/- 4.0 is
     electrostatic (TDX deletes Tyr67 and with it the partner's -1 e), leaving
       J |cos alpha| = 19.7 +/- 4.0 cm^-1
     Constrains |cos alpha| alone.

Uncertainties: (1) and (4) are propagated from the quoted experimental errors
and the J spread. (2) has no published error bar -- each trace is an average of
4 replicates -- so R0 is assigned +/- 0.02, which is generous for a TRA
intercept and is the dominant assumption in this script. (3)'s formal inversion
error is 0.81 cm^-1, which is plainly an underestimate of the lineshape-model
systematic; the gap between the constrained and unconstrained fits (253 vs 365)
is used as its real uncertainty instead.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

KT = 208.509
J0, J_SD = 32.82, 1.554
PHI, PHI_SD = 0.57, 0.03
TAU_PS = 3026.0
DTAU, DTAU_SD = 57.0, 4.0
R_TDX, R_TD, R_SD = 0.52, 0.30, 0.02
CD_DELTA, CD_DELTA_SD = 253.2, 55.9          # 253.2 vs 364.9 -> spread/2 as sigma
SHIFT, SHIFT_SD = 19.7, 4.0                  # cm^-1, electrostatics removed

OUR = {"crystal register": 0.349, "MD ensemble 67.6 ns": 0.213}
OUR_DELTA = (576.4, 71.0)


def superradiance(cos_a, delta, j=J0, phi=PHI):
    omega = np.sqrt(delta ** 2 + 4.0 * j ** 2)
    return phi * np.tanh(omega / (2.0 * KT)) * (2.0 * j / omega) * cos_a


def anisotropy_cos():
    """|cos alpha| and 1-sigma from the R0 drop."""
    def invert(r_td, r_tdx):
        d2 = 2.0 * r_td / r_tdx - 1.0
        return np.sqrt(max((2.0 * d2 + 1.0) / 3.0, 0.0))
    c = invert(R_TD, R_TDX)
    hi = invert(R_TD + R_SD, R_TDX - R_SD)
    lo = invert(R_TD - R_SD, R_TDX + R_SD)
    return c, (hi - lo) / 2.0


def main():
    cos_exp, cos_exp_sd = anisotropy_cos()
    shift_cos, shift_cos_sd = SHIFT / J0, abs(SHIFT / J0) * np.hypot(SHIFT_SD / SHIFT, J_SD / J0)
    target = DTAU / TAU_PS
    target_sd = target * np.hypot(DTAU_SD / DTAU, PHI_SD / PHI)

    print("Constraints on (|cos alpha|, |Delta|), two-state model, "
          f"J = {J0} +/- {J_SD} cm^-1\n")
    print(f"  1 superradiance   Dtau/tau = {target:.5f} +/- {target_sd:.5f}   (couples both)")
    print(f"  2 anisotropy      |cos alpha| = {cos_exp:.3f} +/- {cos_exp_sd:.3f}")
    print(f"  3 CD splitting    |Delta|     = {CD_DELTA:.1f} +/- {CD_DELTA_SD:.1f} cm^-1")
    print(f"  4 red shift       |cos alpha| = {shift_cos:.3f} +/- {shift_cos_sd:.3f}")
    print(f"\n  computed         |cos alpha| = {OUR['crystal register']:.3f} (crystal) / "
          f"{OUR['MD ensemble 67.6 ns']:.3f} (MD),  |Delta| = {OUR_DELTA[0]:.0f} +/- {OUR_DELTA[1]:.0f}\n")

    cos_grid = np.linspace(0.0, 1.0, 1201)
    del_grid = np.linspace(0.0, 2500.0, 2501)
    C, D = np.meshgrid(cos_grid, del_grid, indexing="ij")

    chi = {
        "superradiance": ((superradiance(C, D) - target) / target_sd) ** 2,
        "anisotropy": ((C - cos_exp) / cos_exp_sd) ** 2 * np.ones_like(D),
        "CD splitting": ((D - CD_DELTA) / CD_DELTA_SD) ** 2 * np.ones_like(C),
        "red shift": ((C - shift_cos) / shift_cos_sd) ** 2 * np.ones_like(D),
    }
    names = list(chi)

    print("Pairwise best fits (chi^2 at the joint optimum of each pair):\n")
    print(f"{'pair':<32}{'|cos a|':>9}{'|Delta|':>10}{'chi^2':>9}")
    for i in range(len(names)):
        for k in range(i + 1, len(names)):
            tot = chi[names[i]] + chi[names[k]]
            a, b = np.unravel_index(np.argmin(tot), tot.shape)
            print(f"{names[i]+' + '+names[k]:<32}{cos_grid[a]:>9.3f}{del_grid[b]:>10.0f}{tot[a,b]:>9.2f}")

    total = sum(chi.values())
    a, b = np.unravel_index(np.argmin(total), total.shape)
    dof = len(names) - 2
    print(f"\nGLOBAL fit of all four: |cos alpha| = {cos_grid[a]:.3f} "
          f"(alpha = {np.degrees(np.arccos(-cos_grid[a])):.1f} deg), "
          f"|Delta| = {del_grid[b]:.0f} cm^-1")
    print(f"  chi^2 = {total[a,b]:.1f} for {dof} dof")
    print("  per-observable residuals at that optimum (in sigma):")
    for n in names:
        print(f"    {n:<16}{np.sqrt(chi[n][a,b]):+7.2f}")

    # How badly does each single observable disagree with our computed values?
    print(f"\nOur computed pair (|cos a| = {OUR['crystal register']:.3f}, "
          f"|Delta| = {OUR_DELTA[0]:.0f}) tested against each observable:")
    ia = int(np.argmin(np.abs(cos_grid - OUR["crystal register"])))
    ib = int(np.argmin(np.abs(del_grid - OUR_DELTA[0])))
    for n in names:
        print(f"    {n:<16}{np.sqrt(chi[n][ia,ib]):+7.2f} sigma")

    out = {
        "constraints": {
            "superradiance_dtau_over_tau": [target, target_sd],
            "anisotropy_abs_cos_alpha": [float(cos_exp), float(cos_exp_sd)],
            "cd_splitting_abs_delta_cm": [CD_DELTA, CD_DELTA_SD],
            "red_shift_abs_cos_alpha": [float(shift_cos), float(shift_cos_sd)],
        },
        "global_fit": {"abs_cos_alpha": float(cos_grid[a]),
                       "abs_delta_cm": float(del_grid[b]),
                       "chi2": float(total[a, b]), "dof": dof},
        "computed": {"abs_cos_alpha_crystal": OUR["crystal register"],
                     "abs_cos_alpha_md": OUR["MD ensemble 67.6 ns"],
                     "abs_delta_cm": OUR_DELTA},
    }
    Path("exciton_observables/joint_constraints.json").write_text(json.dumps(out, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    colors = {"superradiance": "tab:red", "anisotropy": "tab:blue",
              "CD splitting": "tab:green", "red shift": "tab:orange"}
    for n in names:
        ax.contourf(D, C, chi[n], levels=[0, 1], colors=[colors[n]], alpha=0.28)
        ax.contour(D, C, chi[n], levels=[1], colors=[colors[n]], linewidths=1.4)
        ax.plot([], [], color=colors[n], lw=6, alpha=0.5, label=f"{n} (1$\\sigma$)")
    ax.errorbar([OUR_DELTA[0]], [OUR["crystal register"]], xerr=[OUR_DELTA[1]],
                fmt="k*", ms=16, capsize=4, label="computed (crystal)")
    ax.errorbar([OUR_DELTA[0]], [OUR["MD ensemble 67.6 ns"]], xerr=[OUR_DELTA[1]],
                fmt="ks", ms=8, capsize=4, mfc="white", label="computed (MD)")
    ax.set_xlabel(r"$|\Delta|$  (cm$^{-1}$)")
    ax.set_ylabel(r"$|\cos\alpha|$")
    ax.set_xlim(0, 1800); ax.set_ylim(0, 0.9)
    ax.set_title("Two-state model: four observables, two unknowns")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig("exciton_observables/joint_constraints.png", dpi=150)
    print("\nWrote exciton_observables/joint_constraints.{json,png}")


if __name__ == "__main__":
    main()

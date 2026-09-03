#!/usr/bin/env python3
"""Extract J from the supermolecular calculation and compare it to J_TDC.

Two independent estimators, both applied to the same pair of roots:

  splitting   J = 0.5 * sqrt(Omega^2 - Delta^2)
              Needs Delta, and at frame 1832 Delta is tiny, so this is the
              estimator that dies first if the numerical floor is high.

  dipole      mu_pm = c_A^pm mu_A + c_B^pm mu_B, solved for the mixing angle
              theta, then J = (Omega/2) sin(2 theta). Needs no Delta at all,
              so it survives a poorly resolved Delta and is the one to trust
              if the two disagree.

THE COMPARISON IS AGAINST J_vacuum, NOT J_screened. A supermolecule in vacuum
has no dielectric screening; production divides the Coulomb integral by
eps_opt = 1.77. At frame 1832: J_screened = 24.549, so J_vacuum = 43.45 cm^-1.
Comparing to the screened value would manufacture a 1.77x discrepancy.

Identifying the exciton pair: a charge-transfer state between fragments 24 A
apart carries |mu| ~ 100+ D against the monomers' ~9 D, so contaminants are
rejected on transition-dipole magnitude. Every root is logged so that rejection
is auditable, and the dipole sum rule |mu_+|^2 + |mu_-|^2 = |mu_A|^2 + |mu_B|^2
is checked rather than assumed.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EV_CM = 8065.543937
DEBYE_PER_AU = 2.5417464519

# frame 1832, production def2-SVP values
J_SCREENED_CM = 24.549
EPS_OPT = 1.77
J_VACUUM_CM = J_SCREENED_CM * EPS_OPT
DELTA_PROD_CM = 0.17
CT_DIPOLE_DEBYE = 30.0          # anything above this is inter-fragment CT


sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_terachem_tddft import parse as parse_tc                 # noqa: E402


def roots_of(path: Path):
    """(energies_eV, oscillator strengths, transition dipoles in a.u.).

    Reuses the production parser rather than a second regex, so the
    supermolecule and the ensemble are read out of TeraChem identically.
    """
    parsed = parse_tc(path)
    if parsed["status"] != "converged":
        raise RuntimeError(f"{path.parent.name}: {parsed['status']}")
    rows = sorted(parsed["roots"], key=lambda r: r["root"])
    e = np.array([r["energy_eV"] for r in rows])
    f = np.array([r["oscillator_strength"] for r in rows])
    mu = np.array([r["transition_dipole_au"] for r in rows], dtype=float)
    return e, f, mu


def load_summary(path: Path) -> dict:
    return json.loads(path.read_text())


def exciton_pair(roots: np.ndarray, dipoles: np.ndarray, f: np.ndarray):
    """The two locally-excited bright roots: largest f among non-CT states."""
    mu = np.linalg.norm(dipoles, axis=1) * DEBYE_PER_AU
    local = mu < CT_DIPOLE_DEBYE
    if local.sum() < 2:
        raise RuntimeError(
            f"fewer than two non-CT roots (|mu| < {CT_DIPOLE_DEBYE} D); "
            "the window is CT-contaminated")
    order = np.argsort(np.where(local, f, -np.inf))[::-1]
    pair = sorted(order[:2])
    return pair, mu, local


def j_from_dipoles(mu_plus, mu_minus, mu_a, mu_b, omega_cm):
    """Mixing angle by least squares on mu_pm = c_A mu_A + c_B mu_B."""
    basis = np.column_stack([mu_a, mu_b])
    c_plus, *_ = np.linalg.lstsq(basis, mu_plus, rcond=None)
    c_minus, *_ = np.linalg.lstsq(basis, mu_minus, rcond=None)
    # normalise each exciton vector, then theta from the upper state
    c_plus = c_plus / np.linalg.norm(c_plus)
    c_minus = c_minus / np.linalg.norm(c_minus)
    sin2theta = float(abs(2.0 * c_plus[0] * c_plus[1]))
    return 0.5 * omega_cm * sin2theta, c_plus, c_minus


def main() -> None:
    base = ROOT / "super_test" / "frame_1832"
    prep = load_summary(base / "preparation.json")
    out = {
        "frame": 1832,
        "preparation": prep["jobs"],
        "centroid_separation_A": prep["centroid_separation_A"],
        "reference": {
            "J_screened_cm": J_SCREENED_CM, "eps_opt": EPS_OPT,
            "J_vacuum_cm": J_VACUUM_CM, "Delta_production_cm": DELTA_PROD_CM,
            "predicted_Omega_screened_cm": float(
                np.hypot(DELTA_PROD_CM, 2 * J_SCREENED_CM)),
            "note": ("a vacuum supermolecule is UNSCREENED, so compare against "
                     "J_vacuum = J_screened * eps_opt"),
        },
        "jobs": {},
    }

    def analyse(tag, mu_ref=None):
        path = base / tag / "tddft.out"
        if not path.exists():
            return None
        try:
            e, f, mu = roots_of(path)
        except Exception as error:                              # noqa: BLE001
            out["jobs"][tag] = {"error": repr(error)}
            return None
        pair, mu_mag, local = exciton_pair(e, mu, f)
        lo, hi = pair
        omega_cm = float((e[hi] - e[lo]) * EV_CM)
        rec = {
            "n_roots": int(e.size),
            "exciton_pair_roots": [int(lo) + 1, int(hi) + 1],
            "E_lower_eV": float(e[lo]), "E_upper_eV": float(e[hi]),
            "Omega_cm": omega_cm,
            "f_pair": [float(f[lo]), float(f[hi])],
            "mu_pair_D": [float(mu_mag[lo]), float(mu_mag[hi])],
            "n_ct_rejected": int((~local).sum()),
            "all_roots": [
                {"root": i + 1, "E_eV": float(e[i]), "f": float(f[i]),
                 "mu_D": float(mu_mag[i]), "local": bool(local[i])}
                for i in range(e.size)],
        }
        if mu_ref is not None:
            mu_a, mu_b = mu_ref
            j_dip, c_p, c_m = j_from_dipoles(mu[hi], mu[lo], mu_a, mu_b, omega_cm)
            rec["J_dipole_cm"] = float(j_dip)
            rec["mixing_upper"] = [float(x) for x in c_p]
            rec["mixing_lower"] = [float(x) for x in c_m]
            # sum rule: the pair must carry the monomers' total dipole strength
            rec["dipole_sum_rule"] = float(
                (mu_mag[lo] ** 2 + mu_mag[hi] ** 2)
                / ((np.linalg.norm(mu_a) ** 2 + np.linalg.norm(mu_b) ** 2)
                   * DEBYE_PER_AU ** 2))
        out["jobs"][tag] = rec
        return rec, mu, pair

    # --- the 200 A reference: isolated monomers through the identical path ---
    far = analyse("super_vac_far")
    if far is not None:
        rec, mu, (lo, hi) = far
        mu_a, mu_b = mu[lo], mu[hi]
        delta_vac_cm = rec["Omega_cm"]
        # The SAME dipole estimator applied here must return zero: the fragments
        # are 200 A apart. Whatever it returns IS the numerical floor.
        j_floor, _, _ = j_from_dipoles(mu[hi], mu[lo], mu_a, mu_b,
                                       rec["Omega_cm"])
        out["jobs"]["super_vac_far"]["Delta_vac_cm"] = delta_vac_cm
        out["jobs"]["super_vac_far"]["J_floor_cm"] = float(j_floor)

        true = analyse("super_vac_true", mu_ref=(mu_a, mu_b))
        if true is not None:
            rec_t = true[0]
            om = rec_t["Omega_cm"]
            inside = om ** 2 - delta_vac_cm ** 2
            rec_t["Delta_vac_cm"] = delta_vac_cm
            rec_t["J_splitting_cm"] = (
                float(0.5 * np.sqrt(inside)) if inside > 0 else None)
            rec_t["J_splitting_note"] = (
                None if inside > 0 else
                "Omega < Delta: the splitting is below the resolvable floor")
            for key, val in (("J_splitting_cm", rec_t["J_splitting_cm"]),
                             ("J_dipole_cm", rec_t.get("J_dipole_cm"))):
                if val:
                    rec_t[key.replace("_cm", "_over_J_vacuum")] = float(
                        val / J_VACUUM_CM)

    (ROOT / "results" / "supermolecule_validation.json").write_text(
        json.dumps(out, indent=2) + "\n")

    print(f"  frame 1832, def2-SVPD, 88 atoms, centroid separation "
          f"{out['centroid_separation_A']:.2f} A")
    print(f"  production: J_screened {J_SCREENED_CM:.3f}  ->  J_vacuum "
          f"{J_VACUUM_CM:.2f} cm^-1 (the comparison target)")
    for tag, rec in out["jobs"].items():
        if "error" in rec:
            print(f"  {tag:<15} ERROR {rec['error']}")
            continue
        print(f"  {tag:<15} roots {rec['exciton_pair_roots']} of {rec['n_roots']}"
              f"  Omega {rec['Omega_cm']:9.2f} cm^-1"
              f"  |mu| {rec['mu_pair_D'][0]:.2f}/{rec['mu_pair_D'][1]:.2f} D"
              f"  CT rejected {rec['n_ct_rejected']}")
        for k in ("J_floor_cm", "J_splitting_cm", "J_dipole_cm",
                  "dipole_sum_rule"):
            if rec.get(k) is not None:
                print(f"      {k:<18} {rec[k]:.4f}")
        if rec.get("J_splitting_note"):
            print(f"      NOTE {rec['J_splitting_note']}")
    print(f"\n  wrote {ROOT/'results'/'supermolecule_validation.json'}")


if __name__ == "__main__":
    main()

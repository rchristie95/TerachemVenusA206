#!/usr/bin/env python3
"""Derived-parameter comparison against Cusick et al. 2026.

We have never run two-photon absorption. Omega is Cusick's instrument, not their
result: what comes out of their analysis is a set of ordinary molecular
parameters, and every one of those is something this pipeline produces natively.
So the comparison happens at the level of derived parameters.

Three quantities are computed here, all from geometry that is already on disk:

  delta   the interchromophore angle. Cusick define x_hat || (mu_a - mu_b) and
          y_hat || (mu_a + mu_b), with delta the angle between each monomer's mu
          and x_hat. That reduces to cos(theta_ab) = -cos(2 delta), hence
          delta = 90 - theta_ab/2, where theta_ab is the plain interdipole angle.
          Their values: 31 deg (vdW/crystal dimer), 14-20 deg (tandem).

  kappa   the true 3D Kasha orientation factor. Their eq 13,
          nu_bar = mu^2 (1 + cos^2 delta)/(h c R^3), only follows if R_hat is
          parallel to x_hat -- an assumption of their quasi-2D model, not a
          measurement. We have the real 3D structure, so we can test it: compare
          the true kappa against -(1 + cos^2 delta) evaluated at the SAME delta.

  d(nu)   the 1PA absorption shift, Kasha:
          delta_nu = J (tan^2 delta - 1)/(tan^2 delta + 1) = -J cos(2 delta)
                   = J cos(theta_ab).
          Cusick measure -35.3 +/- 0.4 cm^-1, attributing <=12 cm^-1 to a Stark
          contribution, leaving about -23 cm^-1 excitonic.

Three of these are dimensionless angles, so unlike the coupling itself they are
immune to the epsilon = 1.77 screening-convention fork (see screening_table.py).

SIGN AND PHASE. The overall sign of a transition dipole is an arbitrary phase:
mu_A -> -mu_A leaves all physics unchanged but flips the sign of cos(theta_ab),
of kappa, and of the 1PA shift. We therefore report signed values AND the
phase-robust folded forms, and compare magnitudes. delta is folded to [0, 90].
Nothing here silently takes an absolute value.

UNITS TRAP. coupling_geometry.npz stores mu in atomic units (e bohr) and r in
Angstrom, and records no unit metadata. That split is set in coupling_dcd_steom.py
(:183 applies ANGSTROM_TO_BOHR to the dipole; :217-218 leave the origins in
Angstrom). The closure test in main() re-derives the stored J_pda_cm column from
the raw vectors and will fail loudly if that convention ever changes.

Inputs are read-only. Nothing here overwrites an existing result file.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PRODUCTION = ROOT / "coupling_nvt_production_cr2_1000_20260721"
GEOMETRY_NPZ = PRODUCTION / "coupling_geometry.npz"
SAMPLES_CSV = PRODUCTION / "coupling_samples.csv"
STATIC_JSON = ROOT / "coupling_paper_steom_thermal" / "dipole_geometry.json"
OUT_DIR = Path(__file__).resolve().parent / "results"

# Matching coupling_dcd_steom.py exactly -- see the Bohr/Angstrom note in
# coupling_core.py: a Coulomb sum evaluated in Angstrom converts with
# BOHR_TO_ANGSTROM, not its reciprocal.
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.6313702

# Cusick et al. 2026, Table S2 and main text.
CUSICK = {
    "delta_vdw_deg": 31.0,
    "delta_tandem_deg": (14.0, 20.0),
    "shift_1pa_total_cm": -35.3,
    "shift_1pa_total_err_cm": 0.4,
    "shift_1pa_stark_max_cm": 12.0,
    "gamma0_deg": 22.0,
}
# Their measured shift minus the largest Stark contribution they allow.
CUSICK_EXCITONIC_SHIFT_CM = (
    CUSICK["shift_1pa_total_cm"] + CUSICK["shift_1pa_stark_max_cm"]
)


def derived(mu_a, mu_b, r_a, r_b, epsilon):
    """theta_ab, delta, kappa and the PDA coupling from raw vectors.

    mu in atomic units, r in Angstrom. Arrays are (n, 3); scalars work too via
    atleast_2d. Returns a dict of (n,) arrays.
    """
    mu_a = np.atleast_2d(np.asarray(mu_a, dtype=float))
    mu_b = np.atleast_2d(np.asarray(mu_b, dtype=float))
    r_a = np.atleast_2d(np.asarray(r_a, dtype=float))
    r_b = np.atleast_2d(np.asarray(r_b, dtype=float))

    mag_a = np.linalg.norm(mu_a, axis=1)
    mag_b = np.linalg.norm(mu_b, axis=1)

    r_vec_ang = r_b - r_a                      # A -> B, matching coupling_dcd_steom.py:220
    separation = np.linalg.norm(r_vec_ang, axis=1)
    r_vec = r_vec_ang * ANGSTROM_TO_BOHR
    r = np.linalg.norm(r_vec, axis=1)
    r_hat = r_vec / r[:, None]

    dot = np.einsum("ij,ij->i", mu_a, mu_b)
    cos_theta = dot / (mag_a * mag_b)
    theta = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

    # Unnormalised orientation factor, exactly as coupling_dcd_steom.py:228.
    jdd = dot - 3.0 * np.einsum("ij,ij->i", mu_a, r_hat) * np.einsum(
        "ij,ij->i", mu_b, r_hat
    )
    kappa = jdd / (mag_a * mag_b)
    j_pda_cm = jdd / (r**3 * epsilon) * HARTREE_TO_CM

    # delta = 90 - theta/2 lands in [0, 90] for theta in [0, 180] already, but
    # fold explicitly so a future signed-theta convention cannot leak through.
    delta = 90.0 - theta / 2.0
    delta = np.abs(((delta + 90.0) % 180.0) - 90.0)

    # Their quasi-2D form, evaluated at OUR delta. This is the self-consistency
    # test of eq 13: does the real 3D kappa match the 2D one at the same angle?
    kappa_2d = -(1.0 + np.cos(np.radians(delta)) ** 2)

    return {
        "theta_ab_deg": theta,
        "delta_deg": delta,
        "kappa": kappa,
        "kappa_2d": kappa_2d,
        "kappa_ratio_abs": np.abs(kappa) / np.abs(kappa_2d),
        "separation_A": separation,
        "J_pda_cm": j_pda_cm,
        "mu_mag_au": mag_a,
    }


def shift_1pa_cm(j_cm, delta_deg):
    """Kasha 1PA shift: J (tan^2 d - 1)/(tan^2 d + 1), identically -J cos(2d)."""
    return -j_cm * np.cos(2.0 * np.radians(delta_deg))


def stats(x):
    x = np.asarray(x, dtype=float)
    return {
        "mean": float(x.mean()),
        "std": float(x.std()),
        "min": float(x.min()),
        "p05": float(np.percentile(x, 5)),
        "p50": float(np.percentile(x, 50)),
        "p95": float(np.percentile(x, 95)),
        "max": float(x.max()),
    }


def closure_test(d, samples_csv):
    """Re-derive the stored columns from raw vectors. Gate on everything else."""
    rows = list(csv.DictReader(open(samples_csv)))
    stored = {
        k: np.array([float(row[k]) for row in rows])
        for k in ("J_pda_cm", "angle_deg", "separation_A")
    }
    checks = {
        "J_pda_cm": float(np.abs(d["J_pda_cm"] - stored["J_pda_cm"]).max()),
        "angle_deg": float(np.abs(d["theta_ab_deg"] - stored["angle_deg"]).max()),
        "separation_A": float(np.abs(d["separation_A"] - stored["separation_A"]).max()),
    }
    # J_pda is stored at CSV precision, so 1e-6 is the honest tolerance there.
    tolerances = {"J_pda_cm": 1e-5, "angle_deg": 1e-9, "separation_A": 1e-9}
    failures = {k: v for k, v in checks.items() if v > tolerances[k]}
    if failures:
        raise SystemExit(
            f"CLOSURE TEST FAILED -- units or frame convention has changed: {failures}"
        )
    return {"max_abs_diff": checks, "tolerances": tolerances, "n_frames": len(rows)}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- production tandem ensemble, 1000 frames ---------------------------
    npz = np.load(GEOMETRY_NPZ)
    epsilon = float(npz["epsilon"])
    prod = derived(npz["mu_A"], npz["mu_B"], npz["r_A"], npz["r_B"], epsilon)
    closure = closure_test(prod, SAMPLES_CSV)
    j_cm = np.asarray(npz["J_cm"], dtype=float)
    prod["J_cm"] = j_cm
    prod["shift_1pa_cm"] = shift_1pa_cm(j_cm, prod["delta_deg"])
    frames = np.asarray(npz["frame"], dtype=int)

    # ---- static crystal / vdW dimer ----------------------------------------
    static_raw = json.loads(STATIC_JSON.read_text())
    stat = derived(
        static_raw["mu_A"], static_raw["mu_B"],
        static_raw["r_A"], static_raw["r_B"], epsilon,
    )
    # No TDC coupling is stored alongside this geometry, so the static 1PA shift
    # uses its own PDA coupling rather than a TDC value it does not have.
    stat["shift_1pa_cm"] = shift_1pa_cm(stat["J_pda_cm"], stat["delta_deg"])

    summary = {
        "provenance": {
            "production_npz": str(GEOMETRY_NPZ.relative_to(ROOT)),
            "production_samples": str(SAMPLES_CSV.relative_to(ROOT)),
            "static_json": str(STATIC_JSON.relative_to(ROOT)),
            "static_dimer": static_raw.get("dimer"),
            "static_source_density": static_raw.get("source_density"),
            "epsilon": epsilon,
            "note": (
                "mu in atomic units, r in Angstrom; J_cm and J_pda_cm are ALREADY "
                "screened by epsilon. delta/theta/kappa are dimensionless and "
                "screening-independent."
            ),
        },
        "closure_test": closure,
        "production_tandem": {
            "n_frames": int(len(frames)),
            "theta_ab_deg": stats(prod["theta_ab_deg"]),
            "delta_deg": stats(prod["delta_deg"]),
            "kappa": stats(prod["kappa"]),
            "kappa_2d": stats(prod["kappa_2d"]),
            "kappa_ratio_abs": stats(prod["kappa_ratio_abs"]),
            "separation_A": stats(prod["separation_A"]),
            "J_cm": stats(j_cm),
            "J_pda_cm": stats(prod["J_pda_cm"]),
            "shift_1pa_cm": stats(prod["shift_1pa_cm"]),
        },
        "static_crystal_vdw": {
            k: float(v[0])
            for k, v in stat.items()
        },
        "cusick_reference": CUSICK,
        "cusick_excitonic_shift_cm": CUSICK_EXCITONIC_SHIFT_CM,
    }

    (OUT_DIR / "derived_parameters.json").write_text(json.dumps(summary, indent=2))

    with open(OUT_DIR / "derived_parameters.csv", "w", newline="") as fh:
        cols = [
            "frame", "theta_ab_deg", "delta_deg", "kappa", "kappa_2d",
            "kappa_ratio_abs", "separation_A", "J_cm", "J_pda_cm", "shift_1pa_cm",
        ]
        w = csv.writer(fh)
        w.writerow(cols)
        for i in range(len(frames)):
            w.writerow(
                [frames[i]] + [f"{prod[c][i]:.10g}" for c in cols[1:]]
            )

    # ---- report -------------------------------------------------------------
    d = summary["production_tandem"]
    s = summary["static_crystal_vdw"]
    lo, hi = CUSICK["delta_tandem_deg"]
    print(f"closure test passed  {closure['max_abs_diff']}")
    print()
    print("PRODUCTION TANDEM ENSEMBLE (n=1000)")
    print(f"  theta_ab   {d['theta_ab_deg']['mean']:8.3f} +/- {d['theta_ab_deg']['std']:.3f} deg")
    print(f"  delta      {d['delta_deg']['mean']:8.3f} +/- {d['delta_deg']['std']:.3f} deg"
          f"      Cusick tandem {lo}-{hi} deg")
    print(f"  kappa      {d['kappa']['mean']:8.4f} +/- {d['kappa']['std']:.4f}")
    print(f"  kappa_2d   {d['kappa_2d']['mean']:8.4f}   (their -(1+cos^2 delta) at our delta)")
    print(f"  |ratio|    {d['kappa_ratio_abs']['mean']:8.4f} +/- {d['kappa_ratio_abs']['std']:.4f}"
          f"   -> 2D assumption off by {abs(1-d['kappa_ratio_abs']['mean'])*100:.1f}%")
    print(f"  R          {d['separation_A']['mean']:8.3f} +/- {d['separation_A']['std']:.3f} A")
    print(f"  J (TDC)    {d['J_cm']['mean']:8.3f} +/- {d['J_cm']['std']:.3f} cm^-1  (screened)")
    print(f"  J (PDA)    {d['J_pda_cm']['mean']:8.3f} +/- {d['J_pda_cm']['std']:.3f} cm^-1  (screened)")
    print(f"  1PA shift  {d['shift_1pa_cm']['mean']:8.3f} +/- {d['shift_1pa_cm']['std']:.3f} cm^-1"
          f"   Cusick excitonic ~{CUSICK_EXCITONIC_SHIFT_CM:.1f}")
    print()
    print("STATIC CRYSTAL / vdW DIMER")
    print(f"  theta_ab   {s['theta_ab_deg']:8.3f} deg")
    print(f"  delta      {s['delta_deg']:8.3f} deg      Cusick vdW {CUSICK['delta_vdw_deg']} deg")
    print(f"  kappa      {s['kappa']:8.4f}   vs 2D {s['kappa_2d']:8.4f}"
          f"   -> off by {abs(1-s['kappa_ratio_abs'])*100:.1f}%")
    print(f"  R          {s['separation_A']:8.3f} A")
    print()
    print(f"wrote {OUT_DIR/'derived_parameters.json'}")
    print(f"wrote {OUT_DIR/'derived_parameters.csv'}")


if __name__ == "__main__":
    main()

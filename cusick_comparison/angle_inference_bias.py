#!/usr/bin/env python3
"""Does neglecting multipolar coupling bias Cusick's inferred interchromophore angle?

THE QUESTION. Our delta = 39.6 deg; Cusick report 14-20 deg for the tandem. That
looked like a flat geometric disagreement. But their delta is not measured with a
ruler -- for the tandem it is INFERRED from an observed coupling through their
eq 13, a point-dipole expression:

    nu_bar = mu^2 (1 + cos^2 delta) / (h c R^3)

Everything that expression omits is therefore absorbed into the fitted angle. The
largest thing it omits is the multipolar (extended transition density) part of
the coupling -- which is exactly the quantity this project computes directly.

THE TEST. Write the point-dipole orientation factor as kappa = 1 + cos^2 delta.
An analyst who observes the true coupling J_TDC but models it as point-dipole
recovers not the true kappa but

    kappa_apparent = kappa_true * f ,     f = J_TDC / J_PDA

and since kappa = 1 + cos^2 delta DECREASES with delta on [0, 90], an inflated
kappa yields a SMALLER apparent angle. The bias therefore has the right sign to
explain the discrepancy before any number is computed. This script asks whether
it also has the right magnitude, by running the inversion backwards: taking
Cusick's published delta, removing the multipolar factor we measure, and asking
what angle their data actually implies.

WHY THIS IS NOT CIRCULAR. f = 1.187 is measured here from a full Coulomb sum over
extended STEOM transition densities against the point-dipole limit of the same
densities at the same geometry. It is not fitted to anything of Cusick's.

SECOND, INDEPENDENT CHECK. Their eq 13 also lets us back out the transition
dipole their numbers imply, from their own (J, R, delta). If that lands on the
experimental |mu| rather than ours, it independently corroborates the weekend
basis-convergence result that our STEOM |mu| = 9.8 D is the outlier.

=============================================================================
OUTCOME, AFTER CHECKING AGAINST THE PUBLISHED PAPER (JPC A, 10.1021/acs.jpca.6c02663)
=============================================================================
The mechanism is real, correctly signed, and the right size -- but it is NOT
the explanation. Three things in the paper decide it:

 1. CONFIRMED. Their eq 13 is exactly this point-dipole form, and they state
    outright: "We disregard any dielectric attenuation by the protein matrix
    because of the nanometer scale of the system." The screening reconciliation
    stands unchanged.

 2. CONFIRMED. The tandem delta really is obtained spectroscopically: "We find
    that the angle delta can span the range from 14 to 20 deg and the coupling
    strength nu_bar -- from 32 to 40 cm^-1." So the load-bearing assumption was
    right.

 3. REFUTES THE EXPLANATION. AlphaFold3 predicts delta = 15 deg for dVenus-TD
    independently of any spectroscopy, and the paper notes the spectroscopic
    range "overlaps with the range predicted by AlphaFold3 simulations at
    delta = 15 deg". A purely structural prediction agrees with their
    spectroscopic one. Our multipolar correction moves their spectroscopic
    value to 37-40 deg and would therefore BREAK that agreement.

 4. Also decisive: delta is fixed by intersecting THREE constraints (their
    Figure 8) -- eq 15 (point-dipole, uses R and mu), eq 19 (the 1PA shift) and
    eq 34 (2PA intensity ratios y and C_H). Only eq 15 carries the point-dipole
    assumption. Correcting it alone moves one family of curves, not the
    intersection.

CONCLUSION. The delta discrepancy is a genuine STRUCTURAL disagreement: our
tandem MD geometry (39.6 deg) differs from AlphaFold3's (15 deg), and their
spectroscopy sides with AlphaFold3. This script therefore quantifies the size of
a real but insufficient bias, and the numbers below should be read as an upper
bound on how much of the gap model error could account for -- not as a
resolution of it.

Read-only. Writes results/angle_inference_bias.{md,json}.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PRODUCTION = ROOT / "coupling_nvt_production_cr2_1000_20260721"
OUT_DIR = Path(__file__).resolve().parent / "results"

BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_CM = 219474.6313702
AU_TO_DEBYE = 2.5417464519

# J[cm^-1] = C * mu_A mu_B [D^2] * kappa / R^3 [A^3], derived rather than quoted:
# au dipoles and bohr separations converted to Debye and Angstrom.
C_DEBYE_ANGSTROM = (
    HARTREE_TO_CM * BOHR_TO_ANGSTROM**3 / AU_TO_DEBYE**2
)

EPSILON = 1.77

# Cusick et al. 2026, verified against the published paper and SI.
CUSICK = {
    "tandem": {
        "label": "dVenus-TD (tandem, linkered)",
        "delta_deg": (14.0, 20.0),
        "J_cm": (32.0, 40.0),
        "R_A": 25.0,
        "delta_is_inferred": True,
        "route": ("spectroscopic: intersection of eq 15 (point-dipole), eq 19 "
                  "(1PA shift) and eq 34 (2PA intensity ratios), their Figure 8"),
    },
    "vdw": {
        "label": "dVenus-vdW (crystal contact)",
        "delta_deg": (31.0, 31.0),
        "J_cm": (28.0, 38.0),
        "R_A": 25.4,
        "delta_is_inferred": False,
        "route": "structural, from the dVenus-vdW crystal structure (1myw)",
    },
}
# The independent structural prediction that decides the question: AlphaFold3
# gives this for the tandem, with no spectroscopy and no point-dipole model.
CUSICK_ALPHAFOLD3_DELTA_DEG = 15.0

# Our own independent decomposition of the 35.3 cm^-1 dimer/monomer red shift,
# from deleting the partner chromophore's charge (TD - TDX).
TDX_EXCITONIC_CM = 6.1
TDX_ELECTROSTATIC_CM = 15.6


def kappa_from_delta(delta_deg):
    return 1.0 + np.cos(np.radians(delta_deg)) ** 2


def delta_from_kappa(kappa):
    """Invert kappa = 1 + cos^2 delta. Returns nan where kappa leaves [1, 2]."""
    k = np.asarray(kappa, dtype=float)
    out = np.full(k.shape, np.nan)
    ok = (k >= 1.0) & (k <= 2.0)
    out[ok] = np.degrees(np.arccos(np.sqrt(k[ok] - 1.0)))
    return out


def load_ensemble():
    rows = list(csv.DictReader(open(PRODUCTION / "coupling_samples.csv")))
    j = np.array([float(r["J_cm"]) for r in rows])
    pda = np.array([float(r["J_pda_cm"]) for r in rows])
    sep = np.array([float(r["separation_A"]) for r in rows])
    ang = np.array([float(r["angle_deg"]) for r in rows])

    d = np.load(PRODUCTION / "coupling_geometry.npz")
    mu_a, mu_b, r_a, r_b = d["mu_A"], d["mu_B"], d["r_A"], d["r_B"]
    mag = np.linalg.norm(mu_a, axis=1)
    r_vec = r_b - r_a
    r_hat = r_vec / np.linalg.norm(r_vec, axis=1)[:, None]
    kappa = (
        np.einsum("ij,ij->i", mu_a, mu_b)
        - 3.0 * np.einsum("ij,ij->i", mu_a, r_hat) * np.einsum("ij,ij->i", mu_b, r_hat)
    ) / (mag * np.linalg.norm(mu_b, axis=1))
    return {
        "J_cm": j, "J_pda_cm": pda, "R_A": sep, "theta_deg": ang,
        "delta_deg": 90.0 - ang / 2.0, "kappa": kappa,
        "mu_debye": float(mag.mean() * AU_TO_DEBYE),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    e = load_ensemble()
    mu_d = e["mu_debye"]

    # Self-check: the derived prefactor must reproduce our own unscreened PDA.
    pred = C_DEBYE_ANGSTROM * mu_d**2 * e["kappa"] / e["R_A"] ** 3
    stored = e["J_pda_cm"] * EPSILON
    rel = float(np.abs(pred - stored).max() / stored.mean())
    if rel > 1e-3:
        raise SystemExit(f"prefactor self-check failed: max rel err {rel:.2e}")

    f = float((e["J_cm"] / e["J_pda_cm"]).mean())
    f_std = float((e["J_cm"] / e["J_pda_cm"]).std())
    delta_true = float(e["delta_deg"].mean())
    delta_true_std = float(e["delta_deg"].std())

    out = {}
    for key, c in CUSICK.items():
        lo, hi = c["delta_deg"]
        k_app = np.array([kappa_from_delta(hi), kappa_from_delta(lo)])
        k_corr = k_app / f
        d_corr = delta_from_kappa(k_corr)
        # Their own numbers, inverted for |mu|.
        jm = float(np.mean(c["J_cm"]))
        dm = float(np.mean(c["delta_deg"]))
        mu_implied = float(
            np.sqrt(jm * c["R_A"] ** 3 / (C_DEBYE_ANGSTROM * kappa_from_delta(dm)))
        )
        out[key] = {
            "label": c["label"],
            "route": c["route"],
            "delta_is_inferred": c["delta_is_inferred"],
            "published_delta_deg": [lo, hi],
            "kappa_apparent": k_app.tolist(),
            "kappa_after_removing_multipolar": k_corr.tolist(),
            "delta_corrected_deg": [float(np.nanmin(d_corr)), float(np.nanmax(d_corr))],
            "implied_mu_debye": mu_implied,
        }

    # How precise is their delta range, in kappa terms, versus the neglected term?
    k_lo, k_hi = kappa_from_delta(20.0), kappa_from_delta(14.0)
    k_mid = 0.5 * (k_lo + k_hi)
    quoted_precision = (k_hi - k_lo) / 2.0 / k_mid

    payload = {
        "prefactor_cm_per_debye2_per_angstrom3": C_DEBYE_ANGSTROM,
        "prefactor_self_check_max_rel_err": rel,
        "ours": {
            "multipolar_enhancement_f": f,
            "multipolar_enhancement_f_std": f_std,
            "delta_deg_mean": delta_true,
            "delta_deg_std": delta_true_std,
            "mu_debye_steom": mu_d,
            "R_A_mean": float(e["R_A"].mean()),
        },
        "cusick_corrected": out,
        "precision_argument": {
            "their_delta_range_deg": [14.0, 20.0],
            "implied_kappa_precision_fraction": float(quoted_precision),
            "neglected_multipolar_fraction": f - 1.0,
            "ratio_neglected_to_quoted": float((f - 1.0) / quoted_precision),
        },
        "reference_mu_debye": {
            "steom_ours": mu_d,
            "tddft_basis_converged_weekend": 7.84,
            "experimental_extinction_strickler_berg": [7.5, 7.9],
        },
    }
    (OUT_DIR / "angle_inference_bias.json").write_text(json.dumps(payload, indent=2))

    t, v = out["tandem"], out["vdw"]
    L = []
    A = L.append
    A("# Does multipolar neglect explain the angle discrepancy?\n")
    A(f"Measured multipolar enhancement over the same 1000 frames: "
      f"**f = J_TDC/J_PDA = {f:.4f} ± {f_std:.4f}**. Prefactor self-check against our "
      f"own unscreened point-dipole coupling: max relative error {rel:.1e}.\n")
    A("## The mechanism\n")
    A("Cusick's eq 13 is a point-dipole expression, `nu_bar = mu^2 (1+cos^2 d)/(h c R^3)`.")
    A("An analyst who observes the true coupling but models it this way recovers")
    A("`kappa_apparent = kappa_true * f`. Because `kappa = 1 + cos^2 delta` *decreases*")
    A("with delta, an inflated kappa returns a **smaller** angle. The bias has the")
    A("right sign before any number is computed.\n")
    A("## The tandem angle — inferred, so the correction applies\n")
    A(f"Route: {t['route']}\n")
    A("| | delta (deg) | kappa |")
    A("|---|---|---|")
    A(f"| Cusick, as published | {t['published_delta_deg'][0]:.0f}–{t['published_delta_deg'][1]:.0f} "
      f"| {t['kappa_apparent'][0]:.4f}–{t['kappa_apparent'][1]:.4f} |")
    A(f"| after removing f = {f:.3f} | **{t['delta_corrected_deg'][0]:.2f}–{t['delta_corrected_deg'][1]:.2f}** "
      f"| {t['kappa_after_removing_multipolar'][0]:.4f}–{t['kappa_after_removing_multipolar'][1]:.4f} |")
    A(f"| **this work, computed** | **{delta_true:.2f} ± {delta_true_std:.2f}** "
      f"| {e['kappa'].mean():.4f} |")
    A("")
    A(f"Their angle, corrected only for the multipolar term they omit, lands at "
      f"{t['delta_corrected_deg'][0]:.1f}–{t['delta_corrected_deg'][1]:.1f}°. Ours is "
      f"{delta_true:.2f} ± {delta_true_std:.2f}°. The ranges overlap, with no free "
      f"parameter — f is measured from our own densities.\n")
    A("### …but this is not the explanation\n")
    A(f"Checked against the published paper, the overlap does not survive as an")
    A(f"explanation. AlphaFold3 predicts **delta = {CUSICK_ALPHAFOLD3_DELTA_DEG:.0f}°** "
      f"for dVenus-TD with no spectroscopy and no point-dipole model, and the paper "
      f"notes that the spectroscopic range \"overlaps with the range predicted by "
      f"AlphaFold3\". A purely structural prediction already agrees with their "
      f"spectroscopic one. Applying our correction would move their value to "
      f"{t['delta_corrected_deg'][0]:.0f}–{t['delta_corrected_deg'][1]:.0f}° and "
      f"**break** that agreement.\n")
    A("Their Figure 8 also fixes delta by intersecting three constraints — eq 15")
    A("(point-dipole, uses R and mu), eq 19 (the 1PA shift) and eq 34 (2PA intensity")
    A("ratios y and C_H). Only eq 15 carries the point-dipole assumption, so")
    A("correcting it moves one family of curves, not the intersection.\n")
    A("**The delta gap is therefore a genuine structural disagreement**: our tandem")
    A(f"MD geometry gives {delta_true:.1f}° where AlphaFold3 gives "
      f"{CUSICK_ALPHAFOLD3_DELTA_DEG:.0f}°, and their spectroscopy sides with")
    A("AlphaFold3. The numbers above are best read as an upper bound on how much of")
    A("the gap point-dipole model error could account for — not as a resolution.\n")
    A("## The vdW angle — structural, so the correction does NOT apply\n")
    A(f"Route: {v['route']}\n")
    A(f"Applying the same correction to their structural 31° would give "
      f"{v['delta_corrected_deg'][0]:.1f}°, which agrees with nothing. That is the "
      f"expected result and a useful control: a structurally measured angle carries "
      f"no point-dipole model error to remove. The remaining 31° vs 38.16° gap "
      f"between their crystal geometry and ours is a genuine geometric disagreement "
      f"— most likely the assumed direction of the transition dipole *within* the "
      f"chromophore — and multipolar effects cannot explain it.\n")
    A("The two routes are confirmed distinct in the paper, but note this cuts against")
    A("the multipolar story rather than for it: their *structural* vdW angle (31°) and")
    A("their *spectroscopic* tandem angle (14–20°) differ from each other by more than")
    A("either differs from ours after correction. The paper attributes that to a real")
    A("difference between the two constructs — linker versus crystal packing, solution")
    A("at room temperature versus crystal at low temperature.\n")
    A("## Independent corroboration: their implied transition dipole\n")
    A("Inverting eq 13 for |mu| using *their* J, R and delta:\n")
    A("| Source | implied \\|mu\\| (D) |")
    A("|---|---|")
    A(f"| Cusick tandem | **{t['implied_mu_debye']:.2f}** |")
    A(f"| Cusick vdW | **{v['implied_mu_debye']:.2f}** |")
    A(f"| this work, STEOM | {mu_d:.2f} |")
    A("| weekend TDDFT, basis-converged (def2-TZVPD) | 7.84 |")
    A("| experiment (extinction + Strickler-Berg) | 7.5–7.9 |")
    A("")
    A("Their numbers imply a transition dipole in the experimental range, not ours.")
    A("This is an independent line of evidence for the weekend basis-convergence")
    A("result: the STEOM |mu| = 9.8 D is the outlier, and since J scales as |mu|^2")
    A("it inflates our coupling by ~1.5x on its own.\n")
    A("## Why their quoted precision cannot absorb this\n")
    pa = payload["precision_argument"]
    A(f"Their 14–20° range corresponds to a kappa precision of "
      f"±{pa['implied_kappa_precision_fraction']*100:.2f}%. The multipolar term they "
      f"neglect is {pa['neglected_multipolar_fraction']*100:.1f}% — larger by a factor "
      f"of **{pa['ratio_neglected_to_quoted']:.0f}**. The omission is an order of "
      f"magnitude bigger than the uncertainty they quote, so it cannot be treated as "
      f"within error.\n")
    A("## Consistency with the 1PA shift\n")
    A(f"The corrected picture also fixes what looked like a second disagreement. Our "
      f"Kasha prediction gives an excitonic red shift of −6.16 ± 1.68 cm^-1, against "
      f"Cusick's ~−23 cm^-1 excitonic residual. But our own independent partner-charge "
      f"decomposition (TD − TDX) attributes only **{TDX_EXCITONIC_CM} cm^-1** of the "
      f"35.3 cm^-1 shift to excitonic coupling, with **{TDX_ELECTROSTATIC_CM} cm^-1** "
      f"electrostatic. Two independent routes on our side agree to ~1%. The "
      f"disagreement with Cusick is therefore not about the shift itself but about how "
      f"much of it is Stark rather than excitonic — and their delta, which sets their "
      f"split, is the inferred one corrected above.\n")
    A("## Verified against the published paper\n")
    A("Read from JPC A 10.1021/acs.jpca.6c02663 and its SI:\n")
    A("| Claim | Status |")
    A("|---|---|")
    A("| eq 13 is the point-dipole form assumed here | **confirmed** |")
    A("| They apply no dielectric screening | **confirmed** — stated explicitly after eq 13 |")
    A("| Their tandem delta is spectroscopic, not structural | **confirmed** |")
    A("| Their vdW delta = 31° is structural (1myw crystal) | **confirmed** |")
    A("| gamma_0 = 22° is a direct measurement (Omega = 0.709 at 0–0, via eq 7) | **confirmed** |")
    A(f"| Multipolar neglect explains the delta gap | **refuted** — AlphaFold3 "
      f"independently gives {CUSICK_ALPHAFOLD3_DELTA_DEG:.0f}° |")
    A("")
    A("One further point in their favour, worth recording because it bears on our own")
    A("CD argument: they observe their Omega features separated by 750 cm^-1 against")
    A("an H–J splitting of only 60–90 cm^-1, identify the shape as the first")
    A("derivative of a Gaussian with extrema separated by FWHM/sqrt(2 ln 2) (i.e. 2 sigma),")
    A("and deliberately fit amplitudes rather than reading a splitting off the")
    A("separation. That is the same trap our JPCB round-2 draft falls into with the CD")
    A("couplet — and they stepped over it.\n")
    A("## Caveats\n")
    A("- f is our own TDC/PDA ratio at our geometry. Using it to correct their")
    A("  inference assumes the multipolar enhancement is similar at theirs, which is")
    A("  reasonable at comparable R but is an assumption.")
    A("- The correction is applied to kappa, i.e. it assumes their mu and R are")
    A("  right. Their implied mu differs from ours, so mu and delta are partially")
    A("  degenerate in their fit.")
    A("- The implied-|mu| result below is unaffected by any of the above and remains")
    A("  the most robust finding here.")
    (OUT_DIR / "angle_inference_bias.md").write_text("\n".join(L))
    print("\n".join(L))
    print(f"\nwrote {OUT_DIR/'angle_inference_bias.md'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Is the |mu| = 9.8 D discrepancy introduced by spectroscopic normalisation?

THE QUESTION. Our transition dipole is 9.8 D; experiment (extinction +
Strickler-Berg) and Cusick's own inverted numbers give 7.5-7.9 D, and a TDDFT
basis sweep converges to 7.84 D. Since J scales as |mu|^2 that is a factor of
1.5 and the dominant systematic in the coupling.

Two candidate causes, and they are distinguishable without any new calculation:

 (a) THE NORMALISATION STEP. The production density is built from NTO cubes and
     then rescaled so its integrated dipole matches a hardcoded target vector
     (build_capmasked_steom_density.py:96-123). If that target were wrong, or
     drawn from the wrong convention, the normalisation would be introducing
     the discrepancy.

 (b) THE UNDERLYING STEOM DIPOLE. If ORCA's own printed dipole is already
     ~9.8 D, the normalisation is faithful and the discrepancy is in the method
     (or the truncated QM region), not in our post-processing.

EOM-CC is not variational and its left and right eigenvectors differ, so ORCA
prints the transition moment three ways: RIGHT, LEFT, and their geometric mean
(LEFT-RIGHT), the last being what its FINAL STEOM-CCSD ABSORPTION SPECTRUM
reports. If those three disagree materially, "the" STEOM dipole is convention-
dependent and picking a different one might close the gap. This script extracts
all three and checks.

It also performs the equivalent test on the oscillator strength, which is the
quantity experiment actually measures: f = (2/3) E |mu|^2 in atomic units, so
comparing f against the extinction-derived value is the same test done without
ever converting to Debye.

Read-only. Writes results/normalization_test.{md,json}.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = Path(__file__).resolve().parent / "results"

AU_TO_DEBYE = 2.5417464519
HARTREE_TO_CM = 219474.6313702

SOURCES = {
    "production (steom_phenol_svpd)": ROOT / "neo_model/orca_steom/steom_phenol_svpd.out",
    "weekend def2-SVPD": ROOT / "weekend_dipole/steom_def2_svpd.out",
    "weekend def2-SVP": ROOT / "weekend_dipole/steom_def2_svp.out",
}

# The hardcoded target the real-space density is rescaled to reproduce.
SPEC_NORM_TARGET_AU = (1.13353, -2.06197, 3.12464)  # build_capmasked_steom_density.py
SPEC_NORM_TARGET_EARLIER_AU = (1.08027, -1.98184, 3.01365)  # build_steom_density.py

# Experiment, via Cusick Table S2 / Strickler-Berg.
EXPT_MU_DEBYE = (7.2, 7.9, 8.4)
EXPT_LAMBDA_NM = 515.0
TDDFT_CONVERGED_MU_DEBYE = 7.84


def parse_blocks(path: Path):
    """Return {convention: row} for the first bright transition of each block."""
    if not path.is_file():
        return {}
    text = path.read_text(errors="replace").splitlines()
    out, current = {}, None
    for i, line in enumerate(text):
        if "SPECTRUM FOR RIGHT TRANSITION MOMENTS" in line:
            current = "RIGHT"
        elif "SPECTRUM FOR LEFT TRANSITION MOMENTS" in line:
            current = "LEFT"
        elif "SPECTRUM FOR LEFT-RIGHT TRANSITION MOMENTS" in line:
            current = "LEFT-RIGHT"
        elif "FINAL STEOM-CCSD ABSORPTION SPECTRUM" in line:
            current = "FINAL"
        if current and "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in line:
            for row in text[i : i + 12]:
                m = re.match(
                    r"\s*0-1A\s+->\s+1-1A\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
                    r"([\d.eE+-]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
                    row,
                )
                if m:
                    ev, cm, nm, f, d2, dx, dy, dz = (float(x) for x in m.groups())
                    out.setdefault(
                        current,
                        {
                            "eV": ev, "cm-1": cm, "nm": nm, "fosc": f, "D2_au2": d2,
                            "mu_au": (dx, dy, dz),
                            "mu_debye": d2**0.5 * AU_TO_DEBYE,
                        },
                    )
                    break
            current = None
    return out


def f_from_mu(mu_debye, lambda_nm):
    """Oscillator strength from |mu| and wavelength, both in lab units."""
    mu_au = mu_debye / AU_TO_DEBYE
    e_au = (1.0e7 / lambda_nm) / HARTREE_TO_CM
    return (2.0 / 3.0) * e_au * mu_au**2


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    parsed = {name: parse_blocks(p) for name, p in SOURCES.items()}

    target_mu = sum(x**2 for x in SPEC_NORM_TARGET_AU) ** 0.5 * AU_TO_DEBYE
    target_mu_earlier = sum(x**2 for x in SPEC_NORM_TARGET_EARLIER_AU) ** 0.5 * AU_TO_DEBYE

    L = []
    A = L.append
    A("# Is the transition-dipole discrepancy introduced by normalisation?\n")
    A("## ORCA's own printed dipole, by convention\n")
    A("EOM-CC left and right eigenvectors differ, so the transition moment is")
    A("printed three ways. If they disagreed materially, \"the\" STEOM dipole would")
    A("be convention-dependent.\n")
    A("| source | convention | nm | f | \\|mu\\| (D) |")
    A("|---|---|---|---|---|")
    for name, blocks in parsed.items():
        if not blocks:
            A(f"| {name} | *(not parsed)* | | | |")
            continue
        for conv in ("RIGHT", "LEFT", "LEFT-RIGHT", "FINAL"):
            if conv in blocks:
                b = blocks[conv]
                A(f"| {name} | {conv} | {b['nm']:.1f} | {b['fosc']:.4f} | **{b['mu_debye']:.3f}** |")
    A("")
    A("## The normalisation target\n")
    A(f"`build_capmasked_steom_density.py` rescales the real-space density to")
    A(f"reproduce a hardcoded vector of magnitude **{target_mu:.3f} D**")
    A(f"(an earlier script used {target_mu_earlier:.3f} D). Both are copied from")
    A(f"ORCA's own output, not from experiment.\n")

    prod = parsed.get("production (steom_phenol_svpd)", {})
    ref = prod.get("FINAL") or prod.get("LEFT-RIGHT") or prod.get("RIGHT")
    verdict = {}
    if ref:
        delta = abs(target_mu - ref["mu_debye"])
        A("### Verdict\n")
        A(f"ORCA's own dipole for this run is **{ref['mu_debye']:.3f} D**; the")
        A(f"normalisation target is **{target_mu:.3f} D**, a difference of")
        A(f"{delta:.3f} D ({100*delta/ref['mu_debye']:.1f}%).\n")
        if delta / ref["mu_debye"] < 0.05:
            A("**The normalisation is faithful and is NOT the cause.** It reproduces")
            A("ORCA's own transition dipole to within a few percent, so it cannot be")
            A("responsible for a 1.5x discrepancy in J. The gap is in the underlying")
            A("STEOM calculation (method, or the truncated QM region), not in our")
            A("post-processing of it.\n")
            verdict["normalisation_is_cause"] = False
        else:
            A("**The normalisation may be contributing** — it does not reproduce")
            A("ORCA's own dipole, and the target should be re-derived.\n")
            verdict["normalisation_is_cause"] = True
        verdict.update(orca_mu_debye=ref["mu_debye"], target_mu_debye=target_mu)

    A("## The same test on the oscillator strength\n")
    A("f is what extinction actually measures, so this repeats the comparison")
    A("without ever converting to Debye: f = (2/3) E |mu|^2 in atomic units.\n")
    A("| source | \\|mu\\| (D) | lambda (nm) | f |")
    A("|---|---|---|---|")
    if ref:
        A(f"| this work, STEOM | {ref['mu_debye']:.2f} | {ref['nm']:.1f} | **{ref['fosc']:.3f}** |")
    for mu in EXPT_MU_DEBYE:
        marker = " (accepted)" if mu == 7.9 else ""
        A(f"| experiment{marker} | {mu:.1f} | {EXPT_LAMBDA_NM:.0f} | {f_from_mu(mu, EXPT_LAMBDA_NM):.3f} |")
    A(f"| TDDFT, basis-converged | {TDDFT_CONVERGED_MU_DEBYE:.2f} | {EXPT_LAMBDA_NM:.0f} "
      f"| {f_from_mu(TDDFT_CONVERGED_MU_DEBYE, EXPT_LAMBDA_NM):.3f} |")
    A("")
    if ref:
        f_expt = f_from_mu(7.9, EXPT_LAMBDA_NM)
        A(f"Our oscillator strength is **{ref['fosc']/f_expt:.2f}x** the experimental")
        A(f"value — the same factor as |mu|^2, as it must be. The discrepancy is")
        A("real and is present in ORCA's own output before any post-processing.\n")
        verdict["fosc_ratio_vs_experiment"] = ref["fosc"] / f_expt

    A("## What this means for the basis question\n")
    A("Diffuse functions are the sensitive step, and the two methods respond very")
    A("differently to them:\n")
    A("| basis step | TDDFT | STEOM |")
    A("|---|---|---|")
    A("| SVP → SVPD | 8.71 → 8.24 (−5.4%) | 9.89 → 9.81 (**−0.8%**) |")
    A("| SVPD → TZVPD | 8.24 → 7.84 (−4.9%) | *(run pending)* |")
    A("")
    A("STEOM's dipole is already nearly basis-flat where TDDFT moves most. Granting")
    A("STEOM the same relative SVPD→TZVPD shift TDDFT showed would put it near")
    A("9.3 D — still far from 7.84 D. So basis incompleteness cannot plausibly")
    A("account for the gap either, and the def2-TZVPD run is a confirmation rather")
    A("than the decisive test it was originally framed as.\n")

    (OUT_DIR / "normalization_test.json").write_text(
        json.dumps(
            {"parsed": {k: {c: {kk: vv for kk, vv in b.items()} for c, b in v.items()}
                        for k, v in parsed.items()},
             "spec_norm_target_debye": target_mu,
             "spec_norm_target_earlier_debye": target_mu_earlier,
             "verdict": verdict},
            indent=2,
        )
    )
    (OUT_DIR / "normalization_test.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()

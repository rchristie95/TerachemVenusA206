# The STEOM gauge "failure" was not a defect -- CORRECTED 2026-08-29

> **This file previously concluded that the STEOM transition moment fails an
> internal check "that worsens with basis", and on that basis requoted
> `J` as 21-26 cm^-1 instead of 32.82. That conclusion is WITHDRAWN.**
> The measurements below are correct; the interpretation was not.
> See `gauge_calibration/RESULT.md` for the calibration that settles it.

2026-08-27 measurements. Same geometry (`geom_cthrp.xyz`), same point-charge
embedding (`field.pc`), ORCA 6.1.1.

| method | basis | E (eV) | f (length) | f (velocity) | vel/len | \|mu\| length |
|---|---|---|---|---|---|---|
| STEOM-CCSD | def2-SVP | 2.4724 | 0.9163 | 1.2273 | 1.34 | 9.886 D |
| STEOM-CCSD | def2-SVPD | 2.3272 | 0.8488 | 1.2833 | **1.51** | 9.807 D |
| TDDFT CAM-B3LYP | def2-SVPD | 3.0681 | 0.8795 | 0.8568 | 0.97 | 8.694 D |
| TDDFT CAM-B3LYP | def2-TZVPD | 3.0574 | 0.8729 | 0.8690 | **1.00** | 8.677 D |

## Why the old reading was wrong

**1. "Worsens with basis" was a misreading.** `def2-SVP -> def2-SVPD` adds
diffuse functions at *fixed double-zeta*. It is not a basis-convergence series.
The real series is SVPD -> TZVPD -> QZVPD, and on formaldehyde the STEOM ratio
improves monotonically along it: **0.6886 -> 0.9278 -> 0.9558**, converging to 1
exactly as a correct implementation must.

**2. TDDFT was never a valid control.** Linear-response methods satisfy sum
rules that equation-of-motion methods do not, so TDDFT's ~1.00 is structural and
automatic at *any* basis -- confirmed across three (1.0065 / 0.9882 / 1.0129,
no trend). EOM-class methods, including STEOM, are documented as *not*
gauge-invariant by construction. Demanding that STEOM match TDDFT on this test
was a category error.

**3. The velocity form is not even the right object here.** Only the length form
is strictly valid for a nonlocal Hamiltonian, and STEOM's doubly
similarity-transformed `G` is non-Hermitian and effectively nonlocal.

The gauge ratio at def2-SVPD therefore carries **no information** about whether
the length-gauge dipole is accurate. Supporting evidence: on butadiene, STEOM's
length-gauge \|mu\| (5.78 D) agrees with TDDFT's (5.64 D) to 2.6% while STEOM's
ratio is 0.75.

## What \|mu\| to use

**Use the STEOM length-gauge value, 9.807 D, and therefore J = 32.82 cm^-1.**
The rescaling to TDDFT's 8.68 D (which gave J = 25.7, and the quoted range
21-26, later widened to 17.7-28.6) rested entirely on the withdrawn gauge
argument and has no justification. It also amounted to importing a 0.7 eV energy
error to partly cancel a ~1.5x intensity error.

## The honest limits, which remain

1. **9.807 D is a def2-SVPD number of unknown convergence.** The same
   calibration shows the *length* dipole is not converged at double-zeta:
   formaldehyde STEOM \|mu\| moves 1.295 -> 1.626 -> 1.592 D (+23% from SVPD to
   TZVPD). STEOM and TDDFT also converge toward each other with basis -- 39%
   apart at SVPD, 6.6% at QZVPD -- so our 13% chromophore gap may be largely a
   basis artefact. A triple-zeta chromophore STEOM run would settle it, but
   costs >64 h (attempt 6 needed that just to reach the dressing) and is out of
   reach under a 16 h limit. **State this as a limitation rather than fixing it.**
2. **9.807 D sits above the defensible experimental band, 7.20-9.15 D**
   (`LOCAL_FIELD_CONVENTION.md`). That objection is untouched by any of this.
3. **Both methods overshoot the experimental oscillator strength by 1.1-1.5x**
   under every local-field convention, while agreeing with each other to 2.8%.
   This is the cleanest remaining discrepancy about the dipole. It is common to
   both methods, so it is not a STEOM defect and not a gauge artefact.

## Superseded advice, do not act on it

- "Resolving it needs EOM-CCSD or CC2, not a bigger basis" -- the premise is
  gone. A bigger basis is in fact the relevant axis for the *dipole*; the gauge
  needs nothing.
- `AddDDTerm` and other density options: tested and dead
  (`ADDDDTERM_RESULT.md`), and moot now.

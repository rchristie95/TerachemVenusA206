# Is the transition-dipole discrepancy introduced by normalisation?

## ORCA's own printed dipole, by convention

EOM-CC left and right eigenvectors differ, so the transition moment is
printed three ways. If they disagreed materially, "the" STEOM dipole would
be convention-dependent.

| source | convention | nm | f | \|mu\| (D) |
|---|---|---|---|---|
| production (steom_phenol_svpd) | RIGHT | 523.9 | 0.8733 | **9.864** |
| production (steom_phenol_svpd) | LEFT | 523.9 | 0.8970 | **9.997** |
| production (steom_phenol_svpd) | LEFT-RIGHT | 523.9 | 0.8850 | **9.930** |
| production (steom_phenol_svpd) | FINAL | 523.9 | 0.8850 | **9.930** |
| weekend def2-SVPD | RIGHT | 532.8 | 0.8378 | **9.744** |
| weekend def2-SVPD | LEFT | 532.8 | 0.8600 | **9.872** |
| weekend def2-SVPD | LEFT-RIGHT | 532.8 | 0.8488 | **9.807** |
| weekend def2-SVP | RIGHT | 501.5 | 0.8999 | **9.797** |
| weekend def2-SVP | LEFT | 501.5 | 0.9330 | **9.975** |
| weekend def2-SVP | LEFT-RIGHT | 501.5 | 0.9163 | **9.886** |

## The normalisation target

`build_capmasked_steom_density.py` rescales the real-space density to
reproduce a hardcoded vector of magnitude **9.942 D**
(an earlier script used 9.570 D). Both are copied from
ORCA's own output, not from experiment.

### Verdict

ORCA's own dipole for this run is **9.930 D**; the
normalisation target is **9.942 D**, a difference of
0.012 D (0.1%).

**The normalisation is faithful and is NOT the cause.** It reproduces
ORCA's own transition dipole to within a few percent, so it cannot be
responsible for a 1.5x discrepancy in J. The gap is in the underlying
STEOM calculation (method, or the truncated QM region), not in our
post-processing of it.

## The same test on the oscillator strength

f is what extinction actually measures, so this repeats the comparison
without ever converting to Debye: f = (2/3) E |mu|^2 in atomic units.

| source | \|mu\| (D) | lambda (nm) | f |
|---|---|---|---|
| this work, STEOM | 9.93 | 523.9 | **0.885** |
| experiment | 7.2 | 515 | 0.473 |
| experiment (accepted) | 7.9 | 515 | 0.570 |
| experiment | 8.4 | 515 | 0.644 |
| TDDFT, basis-converged | 7.84 | 515 | 0.561 |

Our oscillator strength is **1.55x** the experimental
value — the same factor as |mu|^2, as it must be. The discrepancy is
real and is present in ORCA's own output before any post-processing.

## What this means for the basis question

Diffuse functions are the sensitive step, and the two methods respond very
differently to them:

| basis step | TDDFT | STEOM |
|---|---|---|
| SVP → SVPD | 8.71 → 8.24 (−5.4%) | 9.89 → 9.81 (**−0.8%**) |
| SVPD → TZVPD | 8.24 → 7.84 (−4.9%) | *(run pending)* |

STEOM's dipole is already nearly basis-flat where TDDFT moves most. Granting
STEOM the same relative SVPD→TZVPD shift TDDFT showed would put it near
9.3 D — still far from 7.84 D. So basis incompleteness cannot plausibly
account for the gap either, and the def2-TZVPD run is a confirmation rather
than the decisive test it was originally framed as.

# STEOM -> TDDFT site-energy calibration, v5 campaign

Running record. The 500-frame TDDFT ensemble carries the statistics; these
DLPNO-STEOM-CCSD references carry the absolute scale. Each reference uses the
*byte-identical* geometry and 75,035-charge embedding file that TeraChem was
given for the same frame, so the difference between them is method and nothing
else.

## Why this exists

`run_ensemble.py` corrects every TDDFT site energy by one hardcoded constant,
`DELTA_STEOM_EV = -1.06778287`. It was derived on a different geometry and a
different embedding. Whether it transfers -- and whether the correct shift is
even a constant -- is what these four frames test.

## Points so far

| frame | site | TDDFT (eV) | STEOM (eV) | shift (eV) | \|mu\| (D) | f | gauge |
|---|---|---:|---:|---:|---:|---:|---:|
| 900 | A | 3.3583 | 2.4016 (516.3 nm) | **-0.9567** | 9.517 | 0.825 | 1.253 |
| 900 | B | 3.2558 | 2.2517 (550.6 nm) | **-1.0041** | 10.044 | 0.861 | 1.238 |
| 1900 | A | 3.2346 | 2.3010 (538.8 nm) | **-0.9336** | 8.763 | - | - |
| 1900 | B | 3.1708 | 2.0857 (594.4 nm) | **-1.0852** | 7.996 | - | - |
| 2400 | A | 3.2539 | 2.2774 (544.4 nm) | **-0.9765** | 9.507 | 0.781 | 1.391 |
| 1400 | A | 3.2654 | 2.3416 (529.5 nm) | **-0.9238** | 8.760 | - | - |
| 1400 | B | 3.2581 | 2.2468 (551.8 nm) | **-1.0113** | 9.466 | - | - |
| 2400 | B | 3.2402 | 2.3352 (531.0 nm) | **-0.9049** | 7.291 | - | - |

**COMPLETE: 8/8 jobs, four frames x two sites.** Mean shift **-0.9745 eV**,
sd 471 cm^-1, +752 cm^-1 from the pipeline constant.

## Final result, n = 8

| claim | statistic | verdict |
|---|---|---|
| the constant's VALUE is wrong | +752 +/- 167 cm^-1, t = 4.52 on 7 dof, **p = 0.0027** | **supported** |
| the shift VARIES beyond dressing noise | chi^2 = 16.0 on 7 dof, **p = 0.025** | **supported** |
| the variation is SITE-dependent | paired t = -1.14, **p = 0.336** | **rejected** |

`DELTA_STEOM_EV = -1.06778287` should be **-0.9745 eV** on this construct,
trajectory and embedding. And a single constant is the wrong functional form:
the residual scatter, 471 cm^-1, is larger than STEOM's own 312 cm^-1
active-space noise. What drives that scatter is not site identity
(B - A runs -382, -706, -1223, +577 cm^-1 -- three negative, one positive).

### Consequence for the paper

A constant cancels exactly in a difference, so the **value** being wrong moves
absolute site energies (~750 cm^-1, ~19 nm) and leaves the detuning, the
coupling and the design rule untouched.

The **variation** would not cancel, and it is real at p = 0.025. But it is
unattributed: not site, and four frames cannot resolve geometry dependence from
a 312 cm^-1 noise floor. Quantifying it needs repeat STEOM runs on one fixed
geometry to measure the noise directly, which is a separate experiment. Until
then the defensible statement is that per-site energies carry an additional
~470 cm^-1 of method-transfer uncertainty that the constant-shift construction
does not express.

### Read this section's history as a warning

The conclusion inverted twice before settling:

| n | "constant's value is wrong" | "shift varies" | "site-dependent" |
|---|---|---|---|
| 3 | supported, t = 5.0 | retracted (290 vs 312) | not tested |
| 4 | **retracted**, p = 0.118 | **supported**, p = 0.029 | proposed |
| 6 | **supported again**, p = 0.017 | supported, p = 0.021 | **dead**, p = 0.578 |
| **8** | **supported, p = 0.0027** | **supported, p = 0.025** | **rejected, p = 0.336** |

Three or four points at this noise level will order themselves into whatever
story is asked of them. Only the n = 8 row should be quoted.

### The ensemble-averaged dipole falls inside the experimental band

\|mu\| over the eight points: **8.92 +/- 0.32 D** (range 7.29-10.04).

The long-standing objection to the STEOM dipole is that the single-geometry
def2-SVPD value, 9.807 D, sits *above* the defensible experimental band of
7.20-9.15 D (`dipole_gauge/LOCAL_FIELD_CONVENTION.md`). The thermally averaged
value does not: 8.92 +/- 0.32 D is inside that band, and within 0.24 D of the
TDDFT value 8.68 D that three other observables already preferred.

Treat this as suggestive, not as a resolution: these are def2-SVP and the
9.807 D reference is def2-SVPD, and the calibration series shows \|mu\| moves
23% between double and triple zeta. It does say the gap between STEOM and the
experimental band is smaller for a thermal ensemble than for one relaxed
geometry.

### The transition dipole is steadier than the energy

\|mu\| runs 8.76, 9.51, 10.04 D across the three points while the excitation
energy moves 34 nm. Frame 900 site A and frame 2400 site A agree to 0.01 D
(9.517 vs 9.508) despite a 28 nm difference in their transitions.

## The shift is not the same on the two sites of one frame

This is the finding that matters, and it does not need four frames to state:
the two sites of the *same* geometry want shifts differing by 0.0474 eV =
**382 cm^-1**. Because a constant shift cancels exactly in a difference, the
pipeline necessarily reports the TDDFT detuning unchanged:

| | detuning A-B |
|---|---:|
| TDDFT, and therefore the constant-shift pipeline | 827 cm^-1 |
| STEOM, computed | **1209 cm^-1** |

STEOM makes the dimer **46% more detuned** than the corrected TDDFT does, on
the one frame where both have been computed. That direction reinforces rather
than threatens the paper's conclusion -- more detuning means more localisation --
but it means the quoted `<|Delta|>` is a lower bound if STEOM is the better
method, and the constant-shift construction cannot see the difference at all.

**Weigh it against the noise before believing it.** 382 cm^-1 is only ~1.2x
STEOM's own 312 cm^-1 active-space irreproducibility, and this is one frame.
Three more frames are queued; if the sign and rough size hold across them it is
real, and if the spread stays comparable to 312 cm^-1 it is not distinguishable
from the dressing landing differently each time.

## What the first point says

**The constant is off by +896 cm^-1 here.** Measured shift -0.9567 eV against
the pipeline's -1.06778 eV. For scale, 896 cm^-1 is *larger than the entire*
`<|Delta|>` of ~621 cm^-1 that this paper is about.

The corrected number is also the physically right one:

| | eV | nm |
|---|---:|---:|
| pipeline constant applied to TDDFT | 2.2905 | 541.3 |
| STEOM, actually computed | **2.4016** | **516.3** |
| experiment | | ~515 |

Gauge ratio 1.2533, against 1.5118 for the def2-SVPD reference. Better without
diffuse functions, as the basis series in `gauge_calibration/RESULT.md`
predicted, and consistent with that file's conclusion that the ratio is a
basis artefact rather than a defect.

### Where the 896 cm^-1 comes from

The constant's origin is recorded at `wb97xd3_100fs100/README.md:71`:
`2.366633 - 3.43441587 = -1.06778287` eV. Comparing term by term:

| | STEOM | TDDFT | shift |
|---|---:|---:|---:|
| original calibration | 2.3666 | 3.4344 | -1.06778 |
| v5 frame 900 site A | 2.4016 | 3.3583 | -0.95666 |
| moved by | **+282 cm^-1** | **-614 cm^-1** | **+896 cm^-1** |

The two methods did not move together, which is the whole reason a single
additive constant is suspect. But the two halves have very different standing:

* The **TDDFT** half, -614 cm^-1, is solid. TeraChem was shown tonight to be
  bit-reproducible on identical input (`SELFINTERACTION_VERIFICATION.md`), so
  this is a real response to the changed geometry and embedding.
* The **STEOM** half, +282 cm^-1, is *within* STEOM's own documented
  active-space irreproducibility. The same in-protein site energy has come out
  at 523.9 and 532.6 nm on repeat runs -- a 311 cm^-1 spread -- because the
  active-space dressing can converge onto a different high-doubles EA root.

So the discrepancy is not an artefact: the larger half of it is reproducible.
But no single point can separate "the shift genuinely varies with geometry"
from "STEOM landed on a different dressing this time". Only the spread across
the remaining references can, and that is the number to wait for.

**Do not over-read one point.** A single site of a single frame shows the
constant is wrong at that geometry. It cannot show whether the correct shift is
itself constant -- that is what the remaining points are for, and it is the
question that decides whether a per-frame calibration is needed at all.

## Cost, measured

4.65 h per job, 81 GB peak PNO/PAO scratch at def2-SVP. The 250-300 GB recorded
at `weekend_dipole/run_ddterm_chain.sh:59` is for def2-SVPD and def2-TZVPD, not
this basis. Two jobs fit concurrently on this disk; four do not, and running
four is what drove the scratch burn to 699 GB/h and cost frame 1400 (both sites
terminated at 3.82 h to reclaim space; queued for re-run).

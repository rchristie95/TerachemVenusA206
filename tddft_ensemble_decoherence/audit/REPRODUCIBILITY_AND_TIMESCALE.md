# Reproducibility test, and the fast/slow decomposition it made possible

2026-08-27. The `short-fixedmm` campaign was stopped at 149/2000 frames to run
this test before spending four more days of GPU time on it.

## 1. The numbers are not numerical noise

Two halves, both bit-exact:

| test | method | result |
|---|---|---|
| frame preparation | re-run `prepare_production_frame.py` on frames 5/40/90/140 from the DCD, sha256 the outputs | `geometry.xyz` and `mm_charges.dat` **identical**, 16/16 |
| TDDFT solve | re-run TeraChem on byte-identical inputs already on disk | **0.000e+00 cm^-1** on all 7 roots, all f, all characters |

So the frame-to-frame scatter is what the QM/MM model actually predicts for
these geometries. There is nothing to fix in the plumbing.

## 2. Site A is clean; site B alone is broken

| | site A | site B | 1 ns production A / B |
|---|---|---|---|
| sd of site energy | **494.0** | 953.5 | 519.5 / 515.5 |
| mean abs step per 4 fs | 275.0 | 716.3 | -- |
| ACF lag 1 (4 fs) | 0.75 | 0.42 | -- |
| brightest state character | **(93,95) in 149/149** | (93,96) / (93,95) / (93,97) | -- |
| frames with a 2nd root at f > 0.4 | 4/149 | **51/149** | -- |
| selected f, mean / min | 0.979 / 0.567 | 0.828 / 0.254 | -- |

The production ensemble has sd_B/sd_A = 0.99 -- the two sites are statistically
identical, as they must be for a C2 dimer. The 4 fs site A (494) agrees with
both production sites. Only 4 fs site B is anomalous, at 1.93x.

**Production `<|Delta|>` = 576 is therefore NOT contaminated by this defect.**

Mechanism on site B: root 1 is a dark 93->94 charge-transfer state whose energy
swings over ~0.5 eV frame to frame (2.38 / 1.93 / 2.42 eV at frames 5/40/90)
while the bright state barely moves. When that CT manifold crowds the bright
pi-pi*, the oscillator strength splits across two roots (frame 40: 0.700/0.453;
frame 90: 0.650/0.541) and no single-root pick is well defined. Site A does it
too, but in 4 frames instead of 51.

The f-weighted centroid is **not** a general fix: it rescues site B (954 -> 563)
but degrades site A (494 -> 570), and site B's pi-pi* manifold carries 1.205
total oscillator strength against site A's 1.068 -- it is absorbing intensity
from a state that is not the same excitation.

Also noted, not yet chased: both sites have MM point charges 1.13-1.32 A from a
QM atom despite a nominal 1.8 A boundary exclusion. Fixed across frames, so it
is a systematic offset rather than a noise source.

## 3. What site A alone already delivers: the homogeneous linewidth

The campaign was commissioned for the gap autocorrelation, but the linewidth
needs the SITE autocorrelation, and site A's is clean over all 149 frames
(592 fs).

    sd            494.0 cm^-1
    tau_c           8.5 fs   (first zero crossing at 28 fs)
    sigma*tau_c    0.79      (intermediate, mildly motionally narrowed)
    cumulant HWHM  340 cm^-1     static Gaussian limit would be 581

Against the three previous estimates -- 300.5 (Delta_q.V proxy), 611 (discarded,
root-hopping), 41.5 (Nguyen T2) -- and against Nguyen's fitted couplet widths of
293.51 and 439.25, the computed **340 cm^-1 lands between the two observed
widths**. The proxy was closer than expected; the 611 is dead.

## 4. The consequence: most of the disorder is fast, not static

Block-averaging site A over increasing windows:

| window | 4 fs | 8 | 12 | 20 | 32 | 48 | 80 | 120 |
|---|---|---|---|---|---|---|---|---|
| sd of block means | 494 | 463 | 420 | 379 | 306 | 277 | **169** | 195 |
| fraction of variance | 1.00 | 0.88 | 0.72 | 0.59 | 0.38 | 0.32 | 0.12 | 0.16 |

Independent cross-check from the power spectrum, fraction of variance slower
than a cutoff: 0.34 below 1/(81 fs) -> sd 287; 0.11 below 1/(592 fs) -> sd 162.

81 fs is hbar/(2J) at 2J = 65.6 cm^-1 -- the timescale on which detuning is
"static" from the exciton's point of view.

**The slow, genuinely static site disorder is ~170-290 cm^-1, not 494-519.**
The two methods bracket rather than agree (block 169, bootstrap 107-192, over
only 7 blocks; spectral 287), so the range is the honest answer.

Propagated to the gap: sd_static ~ 240-410, `<|Delta|>_static` ~ **190-330**
cm^-1, against the 576 currently used and the <~115 the couplet width requires.
The discrepancy that [[cd-forward-model-verdict]] put at 5x is 1.7-2.9x.

This is the first internally coherent picture:

| | old | now |
|---|---|---|
| homogeneous HWHM | 300 (proxy) | 340 (computed, from the fast motion) |
| static `<|Delta|>` | 576 | 190-330 |
| predicted couplet width | 576 / 586 | homogeneous alone ~340 |
| observed | 293.5 / 439.3 | 293.5 / 439.3 |

The old picture double-counted: it charged the same motion once as homogeneous
dephasing and again as static detuning. That is exactly the hazard
`cumulant_lineshape.py` was written to avoid.

## Caveats

- 592 fs cannot resolve anything slower than itself, so the slow component here
  is a lower bound in principle. But the 1 ns production sd (519) barely exceeds
  the 592 fs sd (494), so little extra variance appears between 592 fs and 1 ns.
  Different trajectories, so that comparison is suggestive, not exact.
- The ACF is not stationary over the window: first-half lag-1 0.82 against
  second-half 0.69. 149 frames is not enough for a converged ACF.
- Long-window blocks have 4-7 samples. The 120 fs point (195) is above the
  80 fs point (169), which is noise, not structure.

## What more frames would buy

Not the gap -- site B has to be fixed before its data is worth anything. What is
undersampled is precisely the slow tail of site A's ACF, which is where the
fast/slow split (and therefore the whole argument above) is weakest.

## Reproduce

    scratchpad/repro/rerun.sh      # byte-identical TDDFT re-run
    scratchpad/repro/reprep.sh     # re-preparation from the DCD

# The CD couplet width closes, but not for the reason anyone expected

2026-08-27. `scripts/cd_width_with_computed_disorder.py`, output
`results/cd_width_with_computed_disorder.json`.

## Result

With every parameter computed and nothing fitted:

| case | sd(Delta) | `<|Delta|>` | HWHM | w1 | w2 | **mean width** | lobe asym |
|---|---|---|---|---|---|---|---|
| OLD -- static 576, HWHM 366 | 686.5 | 576.4 | 366 | 474.7 | 436.7 | **455.7** | 1.80 |
| old disorder, computed HWHM | 686.5 | 576.4 | 340 | 452.7 | 365.9 | **409.3** | 1.48 |
| **PRIMARY -- both computed** | 224.9 | 188.8 | 340 | 343.2 | 342.3 | **342.7** | 1.30 |
| zero static disorder | ~0 | 0.8 | 340 | 337.0 | 336.8 | **336.9** | 1.28 |
| **observed (Table S3)** | | | | 293.5 | 439.2 | **366.4** | **2.34** |

**342.7 against 366.4 -- 6.5%, with no free parameters.** The width
discrepancy that stood at 5x is gone.

Inputs, all computed elsewhere in the repo: HWHM 340 (cumulant on site A's
592 fs QM/MM gap ACF), sigma_site,slow 159 (what the 1 ns ensemble adds beyond
the 592 fs window), S = 0.253 and omega = 1369.7 (TeraChem CAM-B3LYP/6-31G*
frequencies), J = 32.82. Only the band CENTRE is positioned to experiment, and
that sets location, not width.

## What actually fixed it

Not the disorder correction alone. Both errors mattered, and they compounded:

- the old run used HWHM = 366.4, i.e. the mean of the OBSERVED widths, as the
  homogeneous input -- circular, and too broad;
- and it treated `<|Delta|>` = 576 as static when the site-energy correlation
  time is 8.5 fs and 81 fs of averaging (hbar/2J) leaves ~170-290 standing.

Fixing only the linewidth gives 409 (12% too broad); fixing only the disorder
leaves the circularity. Together they give 343.

## The width no longer constrains the detuning -- in either direction

mean couplet width, cm^-1:

| sigma_site | `<|Delta|>` | HWHM 250 | HWHM 300 | **HWHM 340** | HWHM 400 |
|---|---|---|---|---|---|
| 159 | 179 | 258 | 305 | **343** | 400 |
| 200 | 226 | 266 | 311 | **348** | 405 |
| 250 | 282 | 280 | 323 | **359** | 414 |
| 287 | 324 | 292 | 334 | **369** | 422 |

At the computed HWHM the width moves only 343 -> 369 across the entire
plausible detuning range, and the observed 366 sits inside it throughout.
Compare zero disorder (337) with the primary (343): the static detuning
contributes almost nothing.

So the couplet width is a measurement of the HOMOGENEOUS LINEWIDTH, not of the
disorder. Two consequences:

1. **The claim that the couplet width requires `<|Delta|>` <~ 115 is dead.** It
   was an artifact of feeding the observed width in as the homogeneous input.
2. The width cannot now be cited as evidence FOR the detuning-dominated
   mechanism either. It is silent on the question.

## The mechanism survives, at a smaller margin

The exciton sees the 81 fs-averaged detuning, not the instantaneous one
(tau_c = 8.5 fs is fast against hbar/2J = 81 fs, so the fast part motionally
narrows). That gives `<|Delta|>`_static = 189-324 against 2|J| = 65.6, i.e.
**detuning-dominated by 3-5x rather than 9x**. Still localised, not
delocalised.

## What is now the open discrepancy

The **lobe asymmetry: 1.30 computed against 2.34 observed.** It was always
there and it is now the largest unexplained feature of the CD, no longer
masked by the width problem. Note the old parameters gave 1.80, closer to
experiment -- but from a spectrum that was 24% too broad, so that agreement
was not meaningful.

## A bonus: the experiment's own two linewidths disagree

Using HWHM = 1/(2 pi c T2):

| source | value | implies |
|---|---|---|
| Table S3 couplet widths | 293.5 / 439.3 cm^-1 | T2 = 18.1 / 12.1 fs |
| **our computed HWHM** | **340 cm^-1** | **T2 = 15.6 fs** |
| Table S4 coherence time | 128.0 +- 8.1 fs | HWHM 41.5 cm^-1 |
| instrument autocorrelation | 76.7 fs | HWHM 69.2 cm^-1 |

Nguyen's own two linewidth observables differ by **7-11x**, and our computed
value sits squarely on the width one. Their 128 fs is above the IRF, so it is
not simply instrument-limited; the likely reading is that the two probe
different things, with the couplet width carrying inhomogeneous broadening that
the echo removes. Either way the 128 fs cannot be used as a constraint until
that is resolved.

## Caveats

- sigma_site = 159 comes from comparing the 592 fs window against the 1 ns
  ensemble, which are DIFFERENT trajectories. The sensitivity scan covers
  159-287 and the conclusion holds throughout.
- The cumulant lineshape at sigma*tau_c = 0.79 is intermediate Gaussian /
  Lorentzian, while both the forward model and Nguyen's fit use Lorentzians.
- delta comes out at 105 against their 130.79, but delta was already shown
  unidentifiable (R^2 range 4.1e-4), so that comparison remains meaningless.

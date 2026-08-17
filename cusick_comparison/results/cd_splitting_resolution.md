# Resolving the CD couplet: is the separation the exciton gap?

Lineshape from the published Table-S3 reconstruction: HWHM 293.51/439.25 cm^-1 about 19322 cm^-1, with the fitted peak-amplitude ratio 1.581.

## Model validation

With the fitted amplitude ratio the model reproduces this repo's own reconstruction of the published fit: a latent gap of 261.58 cm^-1 gives extrema separated by 500.2 cm^-1 against the recorded 500.0 (error 0.2). Dropping the amplitude ratio and using equal-amplitude bands instead gives 461.9 -- a 40 cm^-1 error, and it also hides the effect described next.

## The estimator is not monotonic

Peak-to-peak separation as a function of the latent gap does **not** increase
monotonically. It starts high at zero gap (the pure-derivative limit of two
unequal bands), falls to a minimum, and only then rises:

| Latent gap (cm^-1) | Separation, Lorentzian (cm^-1) |
|---|---|
| 0 | 689.8 |
| 100 | 518.3 |
| 200 | 490.7 |
| 300 | 511.7 |
| 400 | 555.8 |
| 500 | 614.9 |
| 600 | 685.1 |
| 700 | 763.5 |

Zero-gap value 689.8; minimum 490.6 at a gap of 190 cm^-1.

Two consequences. First, the inversion is **two-valued** over the relevant
range -- a measured separation does not correspond to a unique gap. Second,
a bisection search silently returns one arbitrary branch; the earlier draft
of this analysis did exactly that and reported a single spurious answer.

| Target separation | Gaps that reproduce it (Lorentzian) |
|---|---|
| Kim's 548 | 69, 385 |
| Repo reconstruction 505 | 122, 279 |

## Verdict

The JPCB round-2 reading identifies Omega with the couplet separation. Feed Omega = 553.9 cm^-1 through the actual lineshape and the couplet comes out separated by **652 cm^-1**, not 548. The claimed parameter-free match is off by 104 cm^-1 (19%), and the agreement it reports exists only because the two quantities were equated rather than computed.

The same conflation is already visible in this repo's own numbers: a latent gap of 261.58 cm^-1 produces extrema separated by 505 cm^-1 -- an inflation of 1.93x. That is the step the 'centrepiece' takes.

So the separation-based reading does not survive. But neither does the strong
form of the derivative reading, which says the separation is *pure* linewidth
and carries nothing: the curve above plainly does depend on the gap, just
non-monotonically and two-valued. The accurate statement is stronger than
either draft:

> The couplet peak-to-peak separation is a **non-monotonic, two-valued and
> lineshape-dominated** function of the exciton gap. It is not a usable
> estimator of the gap in either direction, and no coupling should be read
> off it.

What does carry the coupling is the couplet AMPLITUDE, suppressed by 2|J|/Omega: 0.118 at the JPCB Omega, 0.251 at the constrained gap. That is where the analysis should go.

## What this does NOT settle

- The answer depends on the lineshape. Under a purely Gaussian profile the
  zero-gap separation is 533 cm^-1 and the
  admissible gaps shift substantially. The real profile is Voigt-like and
  was not fitted here.
- Delta remains inconsistent across the repo (253 / 550 / 576 cm^-1). This
  analysis favours the low branch but does not by itself fix Delta.
- The Table-S3 reconstruction is itself a reconstruction: raw observations
  and the fitted baseline were unavailable (noted in summary.json).
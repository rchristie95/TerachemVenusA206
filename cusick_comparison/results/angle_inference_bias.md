# Does multipolar neglect explain the angle discrepancy?

Measured multipolar enhancement over the same 1000 frames: **f = J_TDC/J_PDA = 1.1872 ± 0.0114**. Prefactor self-check against our own unscreened point-dipole coupling: max relative error 7.1e-09.

## The mechanism

Cusick's eq 13 is a point-dipole expression, `nu_bar = mu^2 (1+cos^2 d)/(h c R^3)`.
An analyst who observes the true coupling but models it this way recovers
`kappa_apparent = kappa_true * f`. Because `kappa = 1 + cos^2 delta` *decreases*
with delta, an inflated kappa returns a **smaller** angle. The bias has the
right sign before any number is computed.

## The tandem angle — inferred, so the correction applies

Route: spectroscopic: intersection of eq 15 (point-dipole), eq 19 (1PA shift) and eq 34 (2PA intensity ratios), their Figure 8

| | delta (deg) | kappa |
|---|---|---|
| Cusick, as published | 14–20 | 1.8830–1.9415 |
| after removing f = 1.187 | **37.15–40.04** | 1.5861–1.6354 |
| **this work, computed** | **39.61 ± 1.36** | 1.5192 |

Their angle, corrected only for the multipolar term they omit, lands at 37.1–40.0°. Ours is 39.61 ± 1.36°. The ranges overlap, with no free parameter — f is measured from our own densities.

### …but this is not the explanation

Checked against the published paper, the overlap does not survive as an
explanation. AlphaFold3 predicts **delta = 15°** for dVenus-TD with no spectroscopy and no point-dipole model, and the paper notes that the spectroscopic range "overlaps with the range predicted by AlphaFold3". A purely structural prediction already agrees with their spectroscopic one. Applying our correction would move their value to 37–40° and **break** that agreement.

Their Figure 8 also fixes delta by intersecting three constraints — eq 15
(point-dipole, uses R and mu), eq 19 (the 1PA shift) and eq 34 (2PA intensity
ratios y and C_H). Only eq 15 carries the point-dipole assumption, so
correcting it moves one family of curves, not the intersection.

**The delta gap is therefore a genuine structural disagreement**: our tandem
MD geometry gives 39.6° where AlphaFold3 gives 15°, and their spectroscopy sides with
AlphaFold3. The numbers above are best read as an upper bound on how much of
the gap point-dipole model error could account for — not as a resolution.

## The vdW angle — structural, so the correction does NOT apply

Route: structural, from the dVenus-vdW crystal structure (1myw)

Applying the same correction to their structural 31° would give 47.2°, which agrees with nothing. That is the expected result and a useful control: a structurally measured angle carries no point-dipole model error to remove. The remaining 31° vs 38.16° gap between their crystal geometry and ours is a genuine geometric disagreement — most likely the assumed direction of the transition dipole *within* the chromophore — and multipolar effects cannot explain it.

The two routes are confirmed distinct in the paper, but note this cuts against
the multipolar story rather than for it: their *structural* vdW angle (31°) and
their *spectroscopic* tandem angle (14–20°) differ from each other by more than
either differs from ours after correction. The paper attributes that to a real
difference between the two constructs — linker versus crystal packing, solution
at room temperature versus crystal at low temperature.

## Independent corroboration: their implied transition dipole

Inverting eq 13 for |mu| using *their* J, R and delta:

| Source | implied \|mu\| (D) |
|---|---|
| Cusick tandem | **7.64** |
| Cusick vdW | **7.87** |
| this work, STEOM | 9.81 |
| weekend TDDFT, basis-converged (def2-TZVPD) | 7.84 |
| experiment (extinction + Strickler-Berg) | 7.5–7.9 |

Their numbers imply a transition dipole in the experimental range, not ours.
This is an independent line of evidence for the weekend basis-convergence
result: the STEOM |mu| = 9.8 D is the outlier, and since J scales as |mu|^2
it inflates our coupling by ~1.5x on its own.

## Why their quoted precision cannot absorb this

Their 14–20° range corresponds to a kappa precision of ±1.53%. The multipolar term they neglect is 18.7% — larger by a factor of **12**. The omission is an order of magnitude bigger than the uncertainty they quote, so it cannot be treated as within error.

## Consistency with the 1PA shift

The corrected picture also fixes what looked like a second disagreement. Our Kasha prediction gives an excitonic red shift of −6.16 ± 1.68 cm^-1, against Cusick's ~−23 cm^-1 excitonic residual. But our own independent partner-charge decomposition (TD − TDX) attributes only **6.1 cm^-1** of the 35.3 cm^-1 shift to excitonic coupling, with **15.6 cm^-1** electrostatic. Two independent routes on our side agree to ~1%. The disagreement with Cusick is therefore not about the shift itself but about how much of it is Stark rather than excitonic — and their delta, which sets their split, is the inferred one corrected above.

## Verified against the published paper

Read from JPC A 10.1021/acs.jpca.6c02663 and its SI:

| Claim | Status |
|---|---|
| eq 13 is the point-dipole form assumed here | **confirmed** |
| They apply no dielectric screening | **confirmed** — stated explicitly after eq 13 |
| Their tandem delta is spectroscopic, not structural | **confirmed** |
| Their vdW delta = 31° is structural (1myw crystal) | **confirmed** |
| gamma_0 = 22° is a direct measurement (Omega = 0.709 at 0–0, via eq 7) | **confirmed** |
| Multipolar neglect explains the delta gap | **refuted** — AlphaFold3 independently gives 15° |

One further point in their favour, worth recording because it bears on our own
CD argument: they observe their Omega features separated by 750 cm^-1 against
an H–J splitting of only 60–90 cm^-1, identify the shape as the first
derivative of a Gaussian with extrema separated by FWHM/sqrt(2 ln 2) (i.e. 2 sigma),
and deliberately fit amplitudes rather than reading a splitting off the
separation. That is the same trap our JPCB round-2 draft falls into with the CD
couplet — and they stepped over it.

## Caveats

- f is our own TDC/PDA ratio at our geometry. Using it to correct their
  inference assumes the multipolar enhancement is similar at theirs, which is
  reasonable at comparable R but is an assumption.
- The correction is applied to kappa, i.e. it assumes their mu and R are
  right. Their implied mu differs from ours, so mu and delta are partially
  degenerate in their fit.
- The implied-|mu| result below is unaffected by any of the above and remains
  the most robust finding here.
# Site-energy detuning and the Venus CD splitting profile

## Outcome

The corrected coupling does **not** rule out the experimental CD profile, but
the profile also does not establish that site-energy detuning is the missing
physics.  For a nondegenerate two-site Hamiltonian,

\[
H_i=\begin{pmatrix}E_{A,i}&J_i\\J_i&E_{B,i}\end{pmatrix},\qquad
\Omega_i=\sqrt{\Delta_i^2+4J_i^2},\qquad
\Delta_i=E_{A,i}-E_{B,i},
\]

the spectral separation constrains the total exciton gap `Omega`, not `2J`
alone.  With the corrected 1000-frame coupling, `J = 32.8165 +/- 1.5540
cm-1`, the constrained experimental separation of `261.58 cm-1` is obtained
with

```text
|Delta| = 253.1918 +/- 0.8073 cm-1 = 0.031392 eV.
```

The unconstrained `511.8/521.7 nm` component pair instead requires mean
`|Delta| = 364.9092 cm-1 = 0.045243 eV`.

This is numerically plausible as an environmental site shift, but it changes
the physical interpretation.  At the constrained solution the mean excitonic
mixing factor `2|J|/Omega` is only `0.2509`; the minor-site population in each
eigenstate is `1.60%` and the participation ratio is `1.033`.  The two states
are therefore about `98.4%` localized, and the interaction-induced CD area is
only about `25%` of the resonant-dimer value before any arbitrary normalization.
Reproducing the spacing is not sufficient: the absolute CD strength must also
be credible.

No TeraChem or ORCA calculation was launched for this investigation, and
neither manuscript was modified.

## What the published fit actually constrains

The constrained dVenus fit in Nguyen et al. Table S3 has the following
parameters:

| Component | Center (cm-1) | Peak amplitude | HWHM (cm-1) |
|---|---:|---:|---:|
| Low-energy couplet member | 19191.21 | +3.59 | 293.51 |
| High-energy couplet member | 19452.79 | -2.27 | 439.25 |
| Third band | 20774.00 | -0.34 | 565.79 |

The two couplet centers are `261.58 cm-1` apart, but their widths are larger
than their separation.  In the reconstructed sum, the positive and negative
extrema are about `500 cm-1` apart.  A visual peak-to-peak measurement is
therefore not the Hamiltonian gap.  Conversely, even the resonant corrected-J
model can display extrema separated by about `518 cm-1` when convolved with
these unequal fitted widths, despite its latent gap being only `65.63 cm-1`.

The apparent peak-height imbalance is also less problematic than it first
looks.  The fitted peak-height ratio is `1.58`, but the unequal widths give an
absolute integrated-area ratio of only `1.057`.  Approximately equal and
opposite integrated areas are consistent with an exciton-chirality couplet.
The third band remains outside an S1-only two-state model.

These numbers reconstruct the published fitted function, not the raw CD
observations.  The published baseline and slope coefficients and raw pointwise
uncertainties are unavailable here, so the numerical residuals below are
shape diagnostics, not confidence intervals.

## Extended profile and identifiability sweep

The calculation used all 1000 corrected `J_i`, transition-dipole vectors, and
centroid geometries from the audited coupling archive.  The geometric triple
product has one handedness sign in all 1000 frames, so ensemble sign
cancellation is not responsible for a weak couplet.

A profile fit was performed against the Table-S3 reconstruction with the
experimental low/high widths fixed.  Overall CD scale, the third-band
amplitude, constant baseline, and linear slope were profiled as nuisance
parameters.  With corrected `J`, the best profile shape occurs at
`|Delta| = 229.94 cm-1`, with mean `Omega = 239.14 cm-1`, and has normalized
RMS residual `0.00993` of the reconstructed profile range.  The difference
from the exact center-derived `253.19 cm-1` reflects broad overlap, the unequal
component widths, and nuisance background; it is not evidence for a different
experimental gap.

More importantly, changing assumed `J` while compensating with `Delta` leaves
the normalized residual essentially unchanged:

| Mean assumed J (cm-1) | Best Delta (cm-1) | Mean gap (cm-1) | Normalized RMS |
|---:|---:|---:|---:|
| 8.20 | 238.56 | 239.13 | 0.0099333 |
| 16.41 | 236.86 | 239.13 | 0.0099333 |
| **32.82** | **229.94** | **239.14** | **0.0099331** |
| 65.63 | 199.90 | 239.20 | 0.0099306 |
| 98.45 | 136.11 | 239.42 | 0.0099201 |

The fitted overall amplitude changes to absorb the changing mixing.  This is
the numerical `J/Delta/amplitude` identifiability ridge.  Component centers
alone give the exact ridge

```text
Delta^2 + 4 J^2 = (261.58 cm-1)^2.
```

An absolute interaction-CD calculation, with validated units and an explicit
treatment of intrinsic monomer magnetic-dipole CD, would help break the ridge.
The present relative reconstruction cannot.

## Differential disorder

If the required mean detuning is present but fluctuates normally, the gap and
branches broaden.  Deterministic Gauss-Hermite averaging gives:

| SD(Delta) (cm-1) | Mean gap (cm-1) | SD(gap) (cm-1) | Single-branch SD (cm-1) |
|---:|---:|---:|---:|
| 0 | 261.58 | 0.78 | 0.39 |
| 25 | 261.65 | 24.19 | 12.10 |
| 50 | 261.92 | 48.21 | 24.10 |
| 100 | 263.84 | 93.86 | 46.93 |
| 200 | 286.45 | 162.39 | 81.19 |
| 400 | 392.85 | 272.21 | 136.10 |

The experimental HWHM values (`294/439 cm-1`) are too broad to impose a tight
upper bound on detuning disorder by themselves.  They include homogeneous,
inhomogeneous, intrinsic-CD, and multiband effects.  Interpreting them as pure
dephasing would give artificial effective times of `18.1/12.1 fs`; those are
not comparable directly with the reported `128 fs` thin-film coherence time
or the manuscript's illustrative `60 fs` value.

A second, physically distinct solution exists for a nominal homodimer.  If the
signed detuning has zero ensemble mean but is normally distributed, an SD of
approximately `308 cm-1` is needed to raise the mean gap to `261.58 cm-1` with
the corrected coupling.  This produces gap SD `175 cm-1`, single-branch SD
`87.7 cm-1`, and mean mixing about `0.39`.  Thus a persistent site bias and
zero-mean dynamic/static disorder can both generate the target mean separation,
but they predict very different gap distributions.  A long, matched NVT
site-energy ensemble can distinguish them; one central-frame calculation
cannot.

## What a STEOM correction can and cannot do

The available one-point calibration is

```text
E_STEOM = E_TDDFT - 1.06778287 eV.
```

Applying the same additive shift to both sites and every frame changes the
absolute spectral origin but cancels exactly from `E_A - E_B`.  It cannot
create the required `0.0314 eV` detuning or repair an excessive one.  It must
be applied in energy units, not as a constant wavelength shift.

The STEOM origin `19088.2 cm-1` needs only a common `+233.8 cm-1` (`0.0290 eV`)
translation to coincide with the constrained-fit center `19322 cm-1`; this
alignment has no effect on splitting or mixing.

A site-dependent shift can alter the gap, but a separate arbitrary offset for
site A and site B from one geometry would confound method error with a real
environmental asymmetry.  A defensible differential calibration requires
several matched TDDFT/STEOM site calculations.  With enough points, fit and
cross-validate either

```text
E_STEOM = a + b E_TDDFT
```

or a delta-learning correction based on local electrostatic/environmental
descriptors.  An affine slope `b` rescales TDDFT gaps; a common intercept does
not.

## Existing TDDFT data are diagnostic only

The uncommitted workspace contains a completed 100-frame, 100-fs TDDFT
campaign.  It reports mean gap `0.18305 eV = 1476 cm-1` with sample SD
`0.11092 eV = 895 cm-1`; its common STEOM shift changes the gap by less than
`4.4e-16 eV`.  This cannot be used as the requested result because it:

- samples a short nonequilibrium window rather than the 1-ns NVT ensemble;
- used a different basis;
- did not represent the opposite CR2 with the corrected complete `-1e` RESP
  charge convention; and
- cannot be joined framewise to the corrected retained `J_i` archive.

The corrected-charge one-frame smoke diagnostic gives an apparent gap of
`1916 cm-1`, but state identity remains ambiguous.  If real, it would produce
only `3.4%` mixing and a gap far larger than the observed couplet.  It is a
state-tracking warning, not evidence against the detuning hypothesis.

## Recommended staged calculation

Do **not** launch all 1000 NVT frames yet.

1. **Resolve frame identity.** Restore the exact DCD recorded by the corrected
   `J_i` archive, or recompute corrected `J_i` and dipole geometry on the DCD
   that will supply the TDDFT frames.  Marginal distributions can be explored
   independently, but a scientific `E_A,i/E_B,i/J_i` spectrum requires one
   coordinate identity per row.

2. **Pass a corrected-charge state-tracking pilot.** Use both sites from at
   least five stratified NVT frames (beginning, middle, end, and geometric/J
   extremes).  Retain seven or more roots.  Require physical transition-
   density or NTO overlap, consistent chromophore configuration, dipole
   direction/phase, no root loss, complete partner-CR2 RESP charges, conserved
   boundary charges, and reproducible retry behavior.

3. **Run a 25--50-frame feasibility ensemble.** This costs approximately
   `0.60--1.19 GPU-hours` for energy-only jobs at the measured smoke rate.
   Estimate the mean and block uncertainty of `Delta`, its distribution and
   sign persistence, correlations with `J` and geometry, and the generalized
   CD profile.  The key target is not merely a mean gap near `262 cm-1`: the
   calculation must also retain enough mixing and produce a credible CD area.

4. **Calibrate sparsely with STEOM.** Select approximately 6--12 site/frame
   points spanning TDDFT energy and gap quantiles.  Compare common-offset,
   affine, and descriptor-based corrections with leave-one-frame-out
   validation.  Do not choose the calibration form by its agreement with the
   experimental splitting.

5. **Authorize 1000 frames only after the pilot supports the hypothesis.** At
   the measured two-site smoke runtime, the energy-only 1000-frame campaign is
   about `23.8 GPU-hours`; transition-density/NTO generation and retries will
   add cost.  Production should remain resumable and fail closed on an
   ambiguous state.  The final spectral calculation must use each matched
   tuple `(E_A,i, E_B,i, J_i, mu_A,i, mu_B,i, r_A,i, r_B,i)` and report both
   latent component centers and extrema of the convolved spectrum.

## Reproduction

Run:

```bash
python terachem_site_energy_cd/investigate_splitting_profile.py
```

The machine-readable results and figure are in
`terachem_site_energy_cd/results/splitting_profile_investigation/`.

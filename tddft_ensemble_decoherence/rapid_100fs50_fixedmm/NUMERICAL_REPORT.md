# Numerical report: corrected 100 fs TDDFT/TDC prototype

## Outcome

All 50 retained NVT frames (2--100 fs at 2 fs cadence) have classified A/B
TDDFT states and signed mutual TDC couplings. Six seven-root boundary cases
(frames 31, 32, 33, 42, 43, 44) were preserved and recovered with nine roots.
There are no final failed, replaced, or interpolated frames.

The canonical Gaussian second-cumulant result is **T2* = 11.628 fs**. The
same-geometry STEOM correction is -1.09007669
eV, and the raw and constant-shifted gap coherences agree to
4.441e-16 eV; a constant
shift changes the absolute spectrum but cannot change differential-gap
fluctuations or T2*. No affine calibration is justified from one matched
STEOM geometry.

## Embedding correction

The first prototype used a fresh 12 A residue selection at every frame. Its MM
net charge changed discontinuously and produced an invalid 8.02 fs result.
That output is preserved as an artifact audit. The corrected campaign holds
the frame-0 atom identities fixed: A has 2,510 charges with net -1.3320 e and
B has 2,439 charges with net -0.8158 e at every frame. Coordinates and
minimum-image placement remain instantaneous. The classical MD system itself
is intact and contains 69,456 atoms; link caps occur only in each 44-atom QM
single point.

## Spectral and coupling statistics

- Site A raw energy: 3.2182 +/- 0.0362 eV
- Site B raw energy: 3.3511 +/- 0.1050 eV
- Mean corrected photon wavelength: A 582.78 nm; B 549.57 nm
- Differential gap SD: 0.0924 eV
- Screened signed J: 78.23 +/- 22.02 cm^-1
- Screened Davydov splitting 2|J|: mean 156.46 cm^-1

Link atoms 42--44 are excluded from compact transition densities and TDC.
TDC symmetry passes at 9.28e-15 maximum
relative error; CPU/GPU subset validation passes at
8.52e-16.

## Decoherence and uncertainty

- Canonical cumulant T2*: 11.628 fs
- Sampling-stride T2* (2/4/8 fs): 11.628 / 11.750 / 10.691 fs
- Positive-sequence bath correlation time: 7.583 fs
- First/last 50 fs: first is right-censored at 48 fs; last is 10.468 fs
- Moving-block bootstrap: finite 95% interval
  9.56--51.82
  fs; 235/1000 resamples
  are right-censored
- Exact-density state-path variants: finite median
  11.09 fs, range
  4.20--11.64
  fs, with one censored path. The shortest variants use aggressive
  density-only weighting and introduce large energy discontinuities.
- Direct signed-J propagation first crosses e^-1 at
  11.752 fs, close to the cumulant
  crossing. The gap-only direct trace narrowly avoids the early threshold and
  recurs before crossing at 64.827 fs;
  this finite-window first-crossing sensitivity is not evidence for a robust
  53 fs coupling-induced change.

## Comparison and interpretation

The canonical geometry-dependent result is shorter than the prior 21.8 fs
minimum-image, 24.0 fs PME, and 24.5 +/- 3.7 fs one-ps-block fixed-charge
estimates, and much shorter than the manuscript's phenomenological 60 fs
sensitivity value. The observed difference is driven mainly by the larger
TDDFT differential-gap fluctuation amplitude (SD 0.0924 eV) and B-site state
mixing, with a short positive correlation time of 7.58 fs. Coupling
fluctuations do not materially move the canonical cumulant/full-propagation
crossing, although the direct trace is strongly recurrent in this short
sample.

This is a rapid CAM-B3LYP/6-31G* prototype, not a replacement for the validated
wB97X-D3/6-311G** benchmark. The 100 fs record is too short for independent
1 ps blocks, slow solvation, or inhomogeneous disorder, and its two 50 fs
halves are not stationary. Therefore 11.63 fs is evidence that
geometry-dependent electronic fluctuations may shorten the fixed-charge
estimate, but it is not yet a precise converged T2*. Population relaxation is
not determined here.

Classical moving water point charges are sufficient for this first
electrostatic-dephasing estimate. Quantum waters are not required to compute
the present classical correlation; an explicit-QM first-shell pilot would test
electronic polarization and charge-transfer corrections. Nuclear quantum
effects would require a quantum-nuclear treatment or a separately justified
quantum correction, not merely adding water molecules to snapshot QM regions.

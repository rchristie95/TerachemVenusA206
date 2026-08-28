# 100 fs ωB97X-D3/6-311G** TDDFT/TDC decoherence report

> **Duplicate, 2026-08-27.** 98% identical to `decoherence_calculation_venus/results/REFERENCE_NUMERICAL_REPORT.md` (they describe the same
> campaign and differ only in one caveat paragraph). Neither is git-tracked.

## Scope and completion

The completed dataset contains 100 consecutive, non-interpolated full-dimer
NVT snapshots at 1 fs cadence (1--100 fs). Every snapshot has two independent
embedded site calculations: chromophore A plus its stacked Tyr and chromophore
B plus its stacked Tyr. All 200 TeraChem TDDFTs completed; there were no failed
frames and no retries.

Each site calculation contains 44 QM atoms: 29 atoms of an anionic CR2
chromophore, 12 atoms of the stacked Tyr phenol, and three correctly directed
link hydrogens. The three link atoms are present in the electronic-structure
calculation but excluded from transition-density partitions, compact density
maps, mutual TDC, and the movie. The full MM dimer is not cut.

## Method, state tracking, and calibration

The electronic-structure method is ωB97X-D3/6-311G**, charge −1, singlet, with
the validated fixed local MM point-charge identities in each instantaneous
full-dimer geometry. Seven roots were retained. This was useful for auditing
competing states, but the final target was not selected by a fixed root number.

The retained bright anionic chromophore ππ* state has dominant 93→95
configuration. There is exactly one 93→95 candidate at every frame and site.
Site A changes numerical root between 2 and 3 (50 frames each); site B is root
2 in all 100 frames. An initial tracker incorrectly allowed the distinct
93→96 state to borrow intensity on site B. That provisional path and its
couplings were retained for audit but not used in the final analysis.
Requiring the correct configuration removed 0.15--0.28 eV one-femtosecond
artifacts. All 12 tracking-weight sensitivity paths then agreed at every
frame. Minimum adjacent full-density overlaps are 0.917 (A) and 0.979 (B).
Complete 7×7 full-density overlap matrices were evaluated at 27 ambiguous
edges/shoulders.

The retained ORCA 6.1.1 STEOM reference is 2.366633 eV (523.90 nm,
f = 0.885045). The same-geometry TeraChem reference is 3.43441587 eV, so

    Delta = -1.06778287 eV.

This constant shift changes absolute photon energies but cannot change
E_A−E_B fluctuations or T2*. The measured invariance error is
4.44×10⁻¹⁶ eV. No energy shift or empirical STEOM scale was applied to the
transition densities or couplings.

## Excitation energies and coupling

Raw mean site energies are 3.2158 eV for A and 3.0327 eV for B. After the
constant STEOM shift, they are 2.1480 eV (mean wavelength 577.3 nm) and
1.9649 eV (mean wavelength 632.6 nm). Their mean differential gap is
0.1831 eV with sample SD 0.1109 eV. Wavelength means are reported only as a
spectroscopic aid; all fluctuation statistics use linear energy units.

The signed vacuum TDC has mean 125.22 cm⁻¹ and SD 28.69 cm⁻¹. With the paper's
optical screening convention, J_screened = J_vacuum/1.77, its mean is
70.74 cm⁻¹ and SD 16.21 cm⁻¹. The screened point-dipole diagnostic is only
18.46 cm⁻¹ on average and is not substituted for TDC. The mean screened
Davydov splitting 2|J| is 141.49 cm⁻¹.

Every frame passed A↔B numerical-symmetry and 1000-point CPU/GPU subset checks.
The exact maximum errors are stored in `provenance.json`. Link atoms are absent
from all coupled densities.

## Decoherence

For X(t) = E_A(t)−E_B(t), the mean-subtracted classical Gaussian
second-cumulant gives

    T2* = 41.477 fs.

The positive-sequence integral of the normalized gap correlation is 3.605 fs;
this is a bath-correlation measure, not T2*. Sampling-stride results are
40.852 fs at 2 fs and 38.718 fs at 4 fs, supporting adequate temporal
resolution at 1 fs.

Direct two-site stochastic-Hamiltonian propagation gives an e⁻¹ crossing of
63.514 fs with the full signed J(t), versus 64.075 fs with J fluctuations
removed. Thus coupling fluctuations shift this particular crossing by only
0.56 fs; the differential gap dominates. The direct-propagation coherence and
the Gaussian cumulant are different observables/approximations and their
crossings should not be conflated.

The moving-block bootstrap has 707 finite estimates and 293 estimates
right-censored by the 99 fs available lag. Its finite-only 95% interval is
10.06--62.57 fs and finite median is 17.20 fs. The first 50-frame half gives
9.26 fs, whereas the last half is still above e⁻¹ at 49 fs. This strong segment
dependence means the 41.48 fs whole-window result is a rapid-prototyping
estimate, not a converged stationary value.

Relative to the earlier fixed-charge results (21.8 fs direct minimum image,
24.0 fs PME, and 24.5 ± 3.7 fs across 1 ps blocks), the present whole-window
TDDFT cumulant is longer by about 17--20 fs. It lies below the manuscript's
phenomenological 60 fs sensitivity value. In this 100-fs sample the difference
comes mainly from the geometry-dependent differential-gap amplitude and its
oscillatory correlation, not from J(t). State-character mistakes would instead
have shortened/deformed the result, which is why the 93→95 audit is essential.

## Interpretation and limitations

This calculation estimates ultrafast pure dephasing from a short classical
trajectory. It does not determine population relaxation. A 100-fs window is
too short for independent 1-ps block estimates, slow conformational disorder,
or 1-ns inhomogeneous broadening. No quantum water treatment is included:
waters enter through the validated classical electrostatic embedding.
Nuclear-quantum effects and explicit QM water polarization therefore remain
limitations.

All frame-level values, failed-frame classification (zero failures), selected
roots, assignment diagnostics, signed couplings, checksums, runtimes, GPU
assignment, validation errors, figures, and restart paths are retained in the
manifests and JSON files in this directory. No trajectory or density value was
interpolated.

The primary 50-second movie uses genuine positive and negative volumetric
isosurfaces, not density-point spheres. It follows the existing paper-video
pipeline: 0.30 Å visualization grid, 0.40 Å Gaussian smoothing, normalized
±0.06 marching-cubes surfaces, and 30% CGO transparency. This rasterization is
for display only and does not enter the TDC or decoherence numerics.

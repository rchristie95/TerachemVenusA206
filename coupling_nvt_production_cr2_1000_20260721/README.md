# Audited 1000-frame VenusA206 tandem production ensemble

This package is the production rerun used for manuscript Figure 5. It
contains 1000 snapshots spanning 1 ns and the mutual transition-density
coupling evaluated for every snapshot.

## Molecular dynamics

- OpenMM NVT at 300 K with a 2 fs timestep, Langevin friction 1 ps^-1
- 500,000 integration steps (1 ns total)
- PME, 1.0 nm nonbonded cutoff, 10 A solvent padding, 0.15 M ionic strength
- OpenCL GPU execution with random seed `20260618`
- 1000 protein/chromophore frames, saved every 500 steps = 1 ps
- no artificial CR2--CR2 interface restraint

The raw and PBC-whole trajectories, rolling checkpoint, solvated system,
topology, and final coordinates are external production inputs and are not
version controlled.  Pass their locations through the named parameters of
`run_production_cr2_1000.ps1`, `extract_cr2_transforms.py`, and
`run_coupling_production_cr2_1000.ps1`.  The PBC correction uses only
whole-residue lattice translations; it does not fit or deform any coordinates.

## CR2 parameter audit

Both CR2 sites use the 29-atom anionic production residue from
`anionic_build/monomer_solv.prmtop`, with charge `-0.999999999 e` at each site.
Per chromophore, the audited transplant restored:

- 22 unconstrained bonds and 10 hydrogen-bond constraints
- 56 angles
- 119 periodic proper/improper torsion terms
- 166 nonbonded exclusions and 1--4 exceptions
- production RESP charge, Lennard--Jones, and peptide-edge parameters

The generic XML only allowed OpenMM to construct the initial system. Every
physical term touching CR2 was zeroed and replaced from the AMBER topology
before minimization or dynamics. `cr2_amber_parameter_transplant.json` records
the counts and charges.

## Geometry and rigid-density placement

The cap-masked STEOM density was placed independently at each site by fitting
the 19 shared CR2 heavy atoms. Across all 1000 frames:

- site A fit RMSD: `0.5471 +/- 0.0545 A` (maximum `0.7050 A`)
- site B fit RMSD: `0.5132 +/- 0.0513 A` (maximum `0.6690 A`)
- chromophore separation: `24.6867 +/- 0.3198 A`
  (range `23.7235..26.0772 A`)
- inter-dipole angle: `100.7780 +/- 2.7223 degrees`
  (range `92.2474..108.7192 degrees`)

All spreads above and below are sample standard deviations.

## Full-grid GPU TDC calculation

- NVIDIA GeForce RTX 4080
- 259,277 cap-masked STEOM density points per monomer
- complete unaggregated `259277 x 259277` mutual Coulomb sum for every frame
- optical dielectric `epsilon = 1.77`
- FP32 OpenCL production integrals; no spatial binning or charge aggregation
- link/cap contributions excluded by the precomputed cap Voronoi mask before
  both the coupling calculation and the movie density map

The resulting ensemble, after correcting the reciprocal-distance unit
conversion described below, is:

- TDC coupling: `32.8165 +/- 1.5540 cm^-1` (`n = 1000`)
- median: `32.8062 cm^-1`
- range: `27.5753..37.8988 cm^-1`
- Davydov splitting: `65.6330 +/- 3.1081 cm^-1`
- point-dipole coupling: `27.6424 +/- 1.2832 cm^-1`
- ratio of ensemble means, TDC/PDA: `1.1872`

The original archived GPU output multiplied a Coulomb sum evaluated with
distances in Angstrom by `ANGSTROM_TO_BOHR = 1.889726`. Reciprocal distances
require `BOHR_TO_ANGSTROM = 0.529177`, so the original TDC values were high by
exactly `3.571065`; PDA values were already evaluated in bohr and are unchanged.
The CSV and JSON values in this directory carry the corrected factor
`0.280028520539078`, and the JSON includes an explicit `unit_correction` audit
record. A method-matched ORCA calculation on a separate crystal-derived
geometry gives corrected TDC/PDA values of `6.09/5.36 meV` and
provides a separate implementation check.

Representative full-grid FP64 checks used exactly the same density and
geometry:

| Zero-based frame | FP32 | FP64 | FP32 - FP64 (cm^-1) |
|---:|---:|---:|---:|
| 0 | 29.3944 | 29.3913 | +0.0031 |
| 499 | 33.8863 | 33.8825 | +0.0038 |
| 999 | 34.4638 | 34.4599 | +0.0038 |

## Nguyen-style numerical CD spectra

`manuscript/make_nguyen_style_spectra.py` converts the 1000 frame-specific
couplings and transition-dipole geometries into the same wavelength axis and
`TDX - TD` subtraction order used by Nguyen et al. The site origin is read
from the converged ORCA STEOM electric-dipole spectrum (`19088.2 cm^-1`,
`523.90 nm`, `f = 0.885045444`) rather than aligned to an experimental peak.
At `T2* = 60 fs`, the unit-corrected numerical ensemble gives mean exciton
centres at `522.9847` and `524.7861 nm`.

The Nguyen-form three-Lorentzian fit places the two active components at
`522.8560` and `524.9160 nm` (HWHM `85.40/85.36 cm^-1`). A test component retained
at Nguyen's `481.3709 nm` position has fitted relative amplitude only
`-9.26e-4`; the S1-only two-state model therefore does not predict that third
experimental band. The calculated ordinate is normalized interaction-induced
CD, not absolute molar ellipticity, because intrinsic monomer/TDX CD is not
available from the present electric-dipole data.

- the full 8000-point wavelength grid is regenerated on demand and deliberately
  not committed because it duplicates the compact coupling/geometry inputs
- `spectra_nguyen_style_audit.json`: input hashes, assumptions, fit parameters,
  residuals, and output paths
- `../manuscript/Fig_Spectra_NguyenStyle_Traces.*`: uncoupled and coupled
  interaction-induced traces
- `../manuscript/Fig_Spectra_NguyenStyle_Difference.*`: numerical difference and
  Nguyen-form fit
- `../manuscript/Fig_Spectra_NguyenStyle_Components.*`: component and residual
  diagnostic

## Files

- `coupling_samples.csv`: all frame-level TDC, PDA, angle, separation, and fit values
- `coupling_distribution.json`: summary statistics and all 1000 TDC samples
- `coupling_geometry.npz`: numerical transition-dipole and centroid geometry
- `precision_check_f64/`: representative full-grid FP64 calculations
- `Fig_Coupling_Histogram_1000.png/pdf`: updated coupling histogram
- `cr2_amber_parameter_transplant.json`: force-field transplant audit
- `../videos/tandem_nvt_production_cr2_1000_steom_density.mp4`: 1000-frame
  cap-masked STEOM-density movie (1920x1440, 12 fps; external generated output)

## Checksums

- raw DCD: `5DA1F8B2CE814DD04467B4D121C3BA70CCC8619852D045120D2E73086DF53E5E`
- PBC-whole DCD: `C8285024076091EB14BCEB2C25A93CD73BAB7FB5675CBA41BF42E7C064425C93`
- corrected coupling CSV: `31A28CEFDA1A1C7932D3C8C2E0CDA1EE6103640843728B1580506BC6B45E826A`
- corrected histogram PNG: `3F86C14F37F990E265704E5795F5B118BC61CF3D4C5C2F824D74EC0A5C6F29B9`
- final 1000-frame MP4: `83D79AC1D49F55A5E42C16AA071CDD7DE8DDC91CD90BA8485F156AF946D09553`

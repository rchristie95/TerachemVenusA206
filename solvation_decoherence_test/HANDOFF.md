# Handoff: tandem-Venus solvation and decoherence investigation

> **Coupling correction (8 August 2026):** The shared-scale spectral example
> below has been regenerated with the corrected production J = 32.82 +/- 1.55
> cm-1 in `reference/orca_validation.json`. The electrostatic-gap traces and
> 22--24 fs dephasing estimates do not depend on the corrected TDC conversion.
> Files retained under `spectra_t2_24fs/` are historical pre-correction output
> and are not used for the values below or for the current manuscript.

## Objective

Test whether explicit QM/MM environmental fluctuations can explain the
solvation and electronic-decoherence timescales assumed in
`manuscript/JPCB_tandem_round_2.tex`, and determine the effect on the reconstructed CD
spectrum.

## Work completed

- Continued the retained full-solvent OpenMM trajectory for two contiguous
  4 ps segments, saving all 69,471 atoms every 4 fs.
- Reconstructed the STEOM S1 excited-minus-ground difference density from
  seven NTO pairs and partitioned it onto the 41 physical QM atoms.
- Calculated site-A, site-B and differential electrostatic energy-gap traces
  against the explicit protein, water and ions.
- Converted their correlation functions into a classical second-cumulant
  coherence envelope.
- Checked 1 ps blocks, sampling cadence, ion removal, NTO truncation, dipole
  matching and full OpenMM PME electrostatics.
- Regenerated the shared-scale CD sensitivity panel with the corrected
  production coupling and the calculated PME dephasing time.

## Main result

The environmental bath and the coherence do not decay on the same timescale:

- site-solvation correlation times: approximately 70--80 fs;
- differential bath correlation: approximately 60--70 fs, with oscillatory
  multicomponent behaviour;
- classical electrostatic pure-dephasing time: **21.8 fs**;
- independent full-PME estimate: **24.0 fs**;
- eight 1 ps blocks: median 24.5 fs, SD 3.7 fs.

Protein electrostatics dominate the differential fluctuations. Water alone is
slower and more common-mode, giving approximately 74 fs dephasing in the full
8 ps analysis. Solvation is therefore the source of the stochastic phase
noise, but the solvation memory time is not itself T2*: decoherence also
depends strongly on the fluctuation amplitude.

## Consequence for the corrected-coupling CD sensitivity

Using 24 fs instead of the assumed 60 fs changes the Lorentzian HWHM from 88.5
to 221.1 cm-1 (FWHM 442.1 cm-1). This exceeds the corrected mean Davydov
splitting of 65.6 cm-1. The exciton centres remain at 522.98 and 524.79 nm, and
the CD sign order survives, but the spectrum broadens substantially, the
extrema move outwards, and the unnormalised peak magnitude falls by about
82.7%. Normalizing
the plotted CD hides this intensity loss.

## Important limitations

- This is a classical fixed-charge cumulant estimate, not a rigorous quantum
  bath calculation or direct experimental T2 prediction.
- It omits geometry-dependent excitation energies, polarizable embedding,
  intrachromophore vibrations, non-Condon effects and population relaxation.
- The exact CR2 charge transplant leaves the existing periodic box at net
  -2e. Removing ion contributions changes T2 by only 0.2 fs, but the production
  box should be rebuilt and neutralized before publication-quality work.
- Eight picoseconds resolves the ultrafast decay but not slow protein or bulk
  solvent modes.

## Key files

- `solvation_decoherence_test/NUMERICAL_TEST_REPORT.md` -- full methodology and results.
- `solvation_decoherence_test/validation_8ps/validation_summary.json` -- block and stride tests.
- `solvation_decoherence_test/validation_8ps/solvation_decoherence_validation.png` -- summary figure.
- `solvation_decoherence_test/pme_validation_8ps/summary.json` -- PME cross-check.
- `solvation_decoherence_test/spectra_t2_24fs/` -- archived reconstruction at 24 fs.
- `notes/Fig_T2_CommonScale.pdf` -- corrected-coupling shared-scale sensitivity panel.
- `build_steom_difference_probe.py` -- constructs the STEOM difference-charge probe.
- `analyze_solvation_decoherence.py` -- computes electrostatic gap traces and cumulant coherence.
- `analyze_pme_decoherence.py` -- full-PME validation.
- `validate_solvation_decoherence.py` -- block, cadence and component analysis.

## Recommended next step

Treat 22--24 fs as a sensitivity result rather than immediately replacing the
manuscript value. For a publishable estimate, compute geometry-dependent
electrostatic-embedding QM/MM excitation gaps on a high-cadence subset,
calibrate a tractable TDDFT/sTDA method against retained STEOM points, build a
quantum-corrected spectral density, and regenerate the CD lineshape without
assuming a single Lorentzian T2.

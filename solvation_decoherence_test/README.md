# Solvation/decoherence reproducibility package

This directory retains the compact inputs and audits for the standalone
solvation/decoherence note.  The results are intentionally not incorporated
into `manuscript/JPCB_tandem.tex`.

## Environment and hardware

Create the declared Conda environment with `conda env create -f environment.yml`.
`mdtraj` streams the DCD files; OpenMM evaluates the fixed-charge/PME cross-check;
NumPy/SciPy perform the correlations and cumulant integration.  The mutual TDC
workflow uses a GPU through OpenCL (`pyopencl`) and was validated with FP64
spot checks.  A suitable OpenCL driver/device must be installed separately.

The following large files are external and must be supplied by the user:

- two contiguous full-system trajectories (69,471 atoms), 0--4 ps and 4--8 ps,
  sampled every 4 fs;
- the matching solvated topology/system files for the PME calculation;
- the 1 ns protein/chromophore DCD and topology used to extract 1000 rigid
  chromophore transforms;
- the cap-masked STEOM density grid used by the GPU TDC calculation.

Raw trajectories, rolling checkpoints, OpenMM systems, ORCA scratch, volumetric
density grids, rendered frames, and movies are deliberately not committed.
Link/cap atoms are excluded when the STEOM difference-density probe and TDC
density are constructed.

## Retained compact data

- `run_4ps/gap_timeseries.npz` and `run_4to8ps/gap_timeseries.npz`: primary
  protein/water/ion electrostatic gap traces;
- `validation_8ps/`: block, cadence, component, and visualization audit;
- `pme_validation_8ps/`: independent full-PME trace and summary;
- `steom_difference_probe*.npz/json`: seven-NTO and sensitivity probes;
- `../coupling_nvt_production_cr2_1000_20260721/`: frame-level coupling and
  geometry data plus FP32/FP64 checks.

## Regeneration

From the repository root, the principal compact-input checks are:

```powershell
conda run --no-capture-output -n venus_qmmm python validate_solvation_decoherence.py `
  --inputs solvation_decoherence_test/run_4ps/gap_timeseries.npz `
           solvation_decoherence_test/run_4to8ps/gap_timeseries.npz `
  --out tmp/decoherence_validation

conda run --no-capture-output -n venus_qmmm python manuscript/make_nguyen_style_spectra.py `
  --data-dir coupling_nvt_production_cr2_1000_20260721 `
  --validation-json reference/orca_validation.json --t2-fs 60 `
  --out-dir tmp/spectra_60fs

python notes/make_decoherence_note_figure.py --out-dir notes
```

The first-principles gap extraction and PME checks additionally require the
external high-cadence trajectories described above; see each CLI's `--help`.

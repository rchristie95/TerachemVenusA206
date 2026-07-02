# TerachemVenusA206

QM/MM computational pipeline for studying Davydov (excitonic) coupling in Venus fluorescent protein (YFP A206) dimers using TeraChem and OpenMM.

## Paper

> **Preprint:** [arXiv:2605.00027](https://arxiv.org/abs/2605.00027)

## Overview

This repository contains scripts to:

1. **Build a crystal-contact dimer** of Venus YFP from the monomeric crystal structure (PDB: [1MYW](https://www.rcsb.org/structure/1MYW)) using crystallographic symmetry operations.
2. **Relax the dimer** classically with OpenMM (energy minimisation and NVT MD).
3. **Run a QM/MM pipeline** with TeraChem — ground-state geometry optimisation, TDDFT excited-state calculations (with optional PCM solvation), and Davydov coupling analysis between the two chromophore sites.
4. **Visualise** transition and difference electron densities on the dimer using PyMOL.

## Repository structure

| File | Description |
|------|-------------|
| `1MYW.pdb` | Input crystal structure of Venus YFP (monomer) |
| `build_dimer.py` | Constructs the closest crystallographic dimer from SMTRY symmetry records |
| `venus_dimer.pdb` | Pre-built dimer output (chains A and B) |
| `minimise_dimer.py` | OpenMM classical minimisation of the dimer with trajectory rendering |
| `run_nvt.py` | OpenMM NVT MD simulation of the dimer with trajectory rendering |
| `qmmm_tddft_pipeline.py` | Full standalone QM/MM pipeline: protonation → solvation → QM boundary → TDDFT → Davydov coupling |
| `reproduce_paper.py` | **One-shot orchestrator** reproducing all numerical paper data: TDDFT (TeraChem, GPU) + STEOM-CCSD (ORCA, CPU) + EOM-CCSD(fT)/ADC(2) doubles/triples (Q-Chem) + static and thermal-NVT Davydov coupling. Config block at the top (seed, ε, frames, per-stage cache reuse); writes `paper_data_summary.json` |
| `align_steom_density.py` | Places the STEOM transition density into the dimer-chain coordinate frame (rigid CR2 Kabsch fit) so the STEOM and TDDFT couplings are evaluated at identical geometry |
| `coupling_core.py` | Reusable, OpenMM-free core: transition-density I/O (`read_dx`), the GPU Coulomb coupling routine (`calculate_coupling`), PyMOL site transforms, and excited-state selection. Imported by the pipeline and all analysis scripts below |
| `coupling_ensemble.py` | Conformational sampling of the coupling: `J` over an MD ensemble → mean ± std + histogram |
| `absorption_cd_spectra.py` | Excitonic absorption + circular-dichroism (CD) lineshapes from `J` + dephasing, overlaid on the experimental Davydov splitting |
| `multipole_analysis.py` | Decomposes the TDC-over-PDA enhancement into dipole–dipole / dipole–quadrupole / quadrupole–quadrupole contributions |
| `open_quantum_dynamics.py` | Open-quantum-system dynamics (Lindblad ME + stochastic Schrödinger equation, Debye-screened `J(t)`); regenerates the dynamics figures and runs T₂*/dielectric sensitivity sweeps. Python port of the MATLAB in `LindbladCodes/` |
| `visualise_dimer.pml` | PyMOL script for visualising TDDFT transition/difference densities on the dimer |

## Dependencies

- Python ≥ 3.9
- [OpenMM](https://openmm.org/) ≥ 8.0
- [PDBFixer](https://github.com/openmm/pdbfixer)
- [NumPy](https://numpy.org/), [SciPy](https://scipy.org/), [Matplotlib](https://matplotlib.org/)
- [Numba](https://numba.pydata.org/) and/or [PyOpenCL](https://documen.tician.de/pyopencl/) (GPU Coulomb kernel for the coupling)
- [PyMOL](https://pymol.org/) (for visualisation)
- FFmpeg (optional; for trajectory video rendering)

External, separately licensed QM back-ends (each needed only for its stage):

- [TeraChem](https://www.petachem.com/) — TDDFT reference (GPU)
- [ORCA](https://www.faccts.de/orca/) 6.1 — DLPNO-STEOM-CCSD site energy and transition density (CPU)
- [Q-Chem](https://www.q-chem.com/) — EOM-CCSD(fT)/ADC(2) doubles/triples validation (CPU)

A conda environment covering the open-source dependencies can be created with:

```bash
conda env create -f environment.yml   # or: pip install -r requirements.txt
conda activate venus_qmmm
```

## One-shot reproduction

To regenerate all of the paper's **numerical** data (site energies, doubles/triples
character, and the static + thermal Davydov couplings) in a single command:

```bash
python reproduce_paper.py            # all stages reuse cached heavy outputs by default
```

Controls live in a config block at the top of the script (and mirror to CLI flags):

| Setting | Default | Purpose |
|---------|---------|---------|
| `SEED` | `20260618` | reproducible NVT integrator + random frame selection |
| `EPS` | `1.77` | optical dielectric screening for the coupling |
| `N_FRAMES` | `200` | NVT frames in the thermal-coupling ensemble (`--n-frames`) |
| `COUPLING_BACKEND` | `opencl` | GPU backend for the TDC kernel |
| `REUSE[...]` | `True` | per-stage cache reuse; flip a stage off with `--run-nvt`, `--run-tddft`, `--run-steom`, `--run-eomft`, `--run-density` |

By default every heavy stage (TeraChem TDDFT, the multi-hour ORCA STEOM-CCSD, the Q-Chem
EOM-CCSD(fT)/ADC(2), and the NVT MD) **reuses cached output**; pass the matching `--run-*`
flag to recompute one from scratch. Results are aggregated to `paper_data_summary.json`
and printed as a table. Requires the `TeraChem` conda environment (OpenMM + TeraChem +
PyMOL + the OpenCL coupling backend); STEOM recompute additionally needs ORCA + the
`openmpi416` MPI environment, and EOM/ADC recompute needs Q-Chem.

## Quick start

### 1. Build the dimer

```bash
python build_dimer.py
# Outputs: venus_dimer.pdb
```

### 2. Relax the dimer

```bash
# Energy minimisation
python minimise_dimer.py

# NVT MD (production run)
python run_nvt.py
```

### 3. Run the full QM/MM pipeline

```bash
python qmmm_tddft_pipeline.py --pdb venus_dimer.pdb
```

Key toggles at the top of `qmmm_tddft_pipeline.py`:

| Variable | Default | Description |
|----------|---------|-------------|
| `USE_CHEAP_METHOD` | `False` | Use HF/3-21G for rapid testing instead of ωB97X-D3/6-311G** |
| `USE_FIXED_RUN_SEED` | `True` | Reproduce stochastic steps with the stored seed |

### 4. Visualise densities

```bash
pymol visualise_dimer.pml
```

The script auto-detects the most recent TDDFT output directory and loads transition/difference density isosurfaces onto the dimer.

## Ensemble & spectroscopic analysis

These scripts build on the single-point pipeline to characterise the coupling distribution, its spectroscopic signature, its multipole origin, and the dimer's open-system dynamics. They share `coupling_core.py`, and all but the MD sampler run without OpenMM/TeraChem (they need only NumPy, SciPy, and Matplotlib; `coupling_ensemble.py` additionally needs PyMOL, and its `--mode full` needs TeraChem).

```bash
# 1. Coupling distribution over an MD trajectory (mean ± std + histogram)
python coupling_ensemble.py --traj tc_dimer_nvt/dimer_nvt_trajectory.pdb \
    --workdir tc_tddft_old_current --monomer tc_simple_old/classical_relaxed.pdb \
    --n-frames 100 --mode rigid
#    -> coupling_sampling_out/{coupling_samples.csv, coupling_distribution.json, Fig_Coupling_Histogram.pdf}

# 2. Absorption + CD lineshape (inhomogeneous width taken from step 1's distribution)
python absorption_cd_spectra.py --distribution coupling_sampling_out/coupling_distribution.json
#    -> lineshape_out/{Fig_Absorption_Spectrum.pdf, Fig_CD_Spectrum.pdf, lineshape_data.csv}

# 3. Multipole decomposition of the TDC-over-PDA enhancement (validate first)
python multipole_analysis.py --self-test
python multipole_analysis.py --workdir tc_tddft_old_current \
    --monomer tc_simple_old/classical_relaxed.pdb --dimer venus_dimer.pdb
#    -> multipole_out/{multipole_analysis.csv, Fig_Multipole_Decomposition.pdf}

# 4. Open-system dynamics: regenerate figures + sensitivity sweeps
python open_quantum_dynamics.py --all
#    -> oqs_out/{Fig_Coupling, Fig_SSE_*, Fig_ME_*, Fig_Bloch_Grid, Fig_T2_Sweep, Fig_Dielectric_Sweep}.pdf
```

Extra dependencies for these scripts: [SciPy](https://scipy.org/) and [Matplotlib](https://matplotlib.org/) (`conda install -c conda-forge scipy matplotlib`).

## Citation

If you use this code, please cite:

```
@misc{christie2026nonequilibrium,
  title         = {{Non-Equilibrium Dynamics of the Time-Dependent Excitonic Coupling in Fluorescent Protein Dimers}},
  author        = {Christie, Robson and Murray, Cerys and Kim, Youngchan and Joo, Jaewoo},
  year          = {2026},
  eprint        = {2605.00027},
  archivePrefix = {arXiv},
  primaryClass  = {physics.chem-ph},
  doi           = {10.48550/arXiv.2605.00027},
  url           = {https://arxiv.org/abs/2605.00027}
}
```

## Data availability

This repository ships the **code** plus the two input structures needed to start a
run (`1MYW.pdb` and the pre-built `venus_dimer.pdb`). The large generated artefacts —
TeraChem/ORCA working directories, MD trajectories, transition-density volumetric
files (`.dx`/`.cube`/`.npz`), figures and videos — are **not** tracked (they exceed
150 GB) and are excluded via `.gitignore`; they regenerate deterministically from the
pipeline commands above. TeraChem and ORCA are third-party packages and are not
redistributed here. The working directory names used in the example commands
(`tc_dimer_nvt/`, `tc_tddft_old_current/`, `tc_simple_old/`, …) are created locally
when you run the pipeline.

## License

No license file is currently included. Add a `LICENSE` (e.g. MIT) before publishing.

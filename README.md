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
| `BuildDimer.py` | Constructs the closest crystallographic dimer from SMTRY symmetry records |
| `venus_dimer.pdb` | Pre-built dimer output (chains A and B) |
| `run_dimer_minimise.py` | OpenMM classical minimisation of the dimer with trajectory rendering |
| `run_dimer_nvt.py` | OpenMM NVT MD simulation of the dimer with trajectory rendering |
| `terachem_full_pipeline.py` | Full standalone QM/MM pipeline: protonation → solvation → QM boundary → TDDFT → Davydov coupling |
| `coupling_core.py` | Reusable, OpenMM-free core: transition-density I/O (`read_dx`), the GPU Coulomb coupling routine (`calculate_coupling`), PyMOL site transforms, and excited-state selection. Imported by the pipeline and all analysis scripts below |
| `sample_coupling_md.py` | Conformational sampling of the coupling: `J` over an MD ensemble → mean ± std + histogram |
| `lineshape_cd.py` | Excitonic absorption + circular-dichroism (CD) lineshapes from `J` + dephasing, overlaid on the experimental Davydov splitting |
| `multipole_decomposition.py` | Decomposes the TDC-over-PDA enhancement into dipole–dipole / dipole–quadrupole / quadrupole–quadrupole contributions |
| `oqs_dynamics.py` | Open-quantum-system dynamics (Lindblad ME + stochastic Schrödinger equation, Debye-screened `J(t)`); regenerates the dynamics figures and runs T₂*/dielectric sensitivity sweeps. Python port of the MATLAB in `LindbladCodes/` |
| `visualise_dimer.pml` | PyMOL script for visualising TDDFT transition/difference densities on the dimer |

## Dependencies

- Python ≥ 3.9
- [OpenMM](https://openmm.org/) ≥ 8.0
- [PDBFixer](https://github.com/openmm/pdbfixer)
- [NumPy](https://numpy.org/)
- [TeraChem](https://www.petachem.com/) (licence required; used for QM/MM)
- [PyMOL](https://pymol.org/) (for visualisation)
- FFmpeg (optional; for trajectory video rendering)

A conda environment covering the open-source dependencies can be created with:

```bash
conda create -n venus_qmmm -c conda-forge python=3.11 openmm pdbfixer numpy
conda activate venus_qmmm
```

## Quick start

### 1. Build the dimer

```bash
python BuildDimer.py
# Outputs: venus_dimer.pdb
```

### 2. Relax the dimer

```bash
# Energy minimisation
python run_dimer_minimise.py

# NVT MD (production run)
python run_dimer_nvt.py
```

### 3. Run the full QM/MM pipeline

```bash
python terachem_full_pipeline.py --pdb venus_dimer.pdb
```

Key toggles at the top of `terachem_full_pipeline.py`:

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

These scripts build on the single-point pipeline to characterise the coupling distribution, its spectroscopic signature, its multipole origin, and the dimer's open-system dynamics. They share `coupling_core.py`, and all but the MD sampler run without OpenMM/TeraChem (they need only NumPy, SciPy, and Matplotlib; `sample_coupling_md.py` additionally needs PyMOL, and its `--mode full` needs TeraChem).

```bash
# 1. Coupling distribution over an MD trajectory (mean ± std + histogram)
python sample_coupling_md.py --traj tc_dimer_nvt/dimer_nvt_trajectory.pdb \
    --workdir tc_tddft_old_current --monomer tc_simple_old/classical_relaxed.pdb \
    --n-frames 100 --mode rigid
#    -> coupling_sampling_out/{coupling_samples.csv, coupling_distribution.json, Fig_Coupling_Histogram.pdf}

# 2. Absorption + CD lineshape (inhomogeneous width taken from step 1's distribution)
python lineshape_cd.py --distribution coupling_sampling_out/coupling_distribution.json
#    -> lineshape_out/{Fig_Absorption_Spectrum.pdf, Fig_CD_Spectrum.pdf, lineshape_data.csv}

# 3. Multipole decomposition of the TDC-over-PDA enhancement (validate first)
python multipole_decomposition.py --self-test
python multipole_decomposition.py --workdir tc_tddft_old_current \
    --monomer tc_simple_old/classical_relaxed.pdb --dimer venus_dimer.pdb
#    -> multipole_out/{multipole_decomposition.csv, Fig_Multipole_Decomposition.pdf}

# 4. Open-system dynamics: regenerate figures + sensitivity sweeps
python oqs_dynamics.py --all
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

## License

See [LICENSE](LICENSE) for details.

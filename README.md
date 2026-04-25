# TerachemVenusA206

QM/MM computational pipeline for studying Davydov (excitonic) coupling in Venus fluorescent protein (YFP A206) dimers using TeraChem and OpenMM.

## Paper

> **Preprint:** [arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX)

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

## Citation

If you use this code, please cite:

```
@misc{christie2025venus,
  title  = {<title>},
  author = {Christie, R. and others},
  year   = {2025},
  eprint = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {physics.chem-ph},
  url    = {https://arxiv.org/abs/XXXX.XXXXX}
}
```

## License

See [LICENSE](LICENSE) for details.

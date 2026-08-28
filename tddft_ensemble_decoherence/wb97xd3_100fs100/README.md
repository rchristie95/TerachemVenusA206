# ωB97X-D3/6-311G** 100 fs paired-site campaign

This directory contains the restartable higher-quality TDDFT/TDC campaign
requested for 100 consecutive full-dimer NVT frames at 1 fs cadence. Each
snapshot has two independent embedded QM calculations: one for chromophore A
plus its stacked Tyr phenol and one for chromophore B plus its stacked Tyr
phenol. The classical dimer is never cut; each 44-atom QM single point contains
41 physical atoms and three correctly directed link hydrogens. Link atoms are
excluded from compact transition densities and mutual TDC.

## NVT trajectory

The trajectory is a new continuation from the retained equilibrated checkpoint,
not an interpolation of the earlier 50-frame/2-fs trajectory:

```bash
conda run --no-capture-output -n TeraChem python -u run_nvt.py \
  --pdb tc_tandem_nvt_rerun_20260722_retry2/solvated_protonated.pdb \
  --skip-solvation \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir tc_tandem_nvt_100fs_100frames_20260724 \
  --no-overwrite-workdir \
  --resume-checkpoint tc_tandem_nvt_rerun_20260722_retry2/nvt.chk \
  --timestep-fs 1.0 --nvt-steps 500100 \
  --openmm-trajectory-interval 1 \
  --openmm-trajectory-file tandem_nvt_100fs_100frames.dcd \
  --checkpoint-file nvt.chk --checkpoint-interval 100 \
  --seed 20260618 --skip-final-minimization --no-video
```

The active post-checkpoint integrator timestep was explicitly reapplied and
verified as 1.000000 fs. The resulting DCD has 100 coordinate records covering
1--100 fs, 69,456 atoms, and periodic box data in every frame.

## TDDFT/TDC workflow

```bash
conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/initialize_manifests.py

conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/register_nvt_trajectories.py \
  --short100fs100-wb97xd3-dcd \
  tc_tandem_nvt_100fs_100frames_20260724/tandem_nvt_100fs_100frames.dcd

conda run --no-capture-output -n TeraChem python -u \
  tddft_ensemble_decoherence/scripts/run_rapid_adaptive_100ns50.py \
  --window short100fs100_wb97xd3 \
  --output-dir tddft_ensemble_decoherence/wb97xd3_100fs100 \
  --maximum-workers 3

conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/global_state_tracking_100fs.py \
  --window short100fs100_wb97xd3 \
  --output-dir tddft_ensemble_decoherence/wb97xd3_100fs100

conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/summarize_100fs50.py \
  --window short100fs100_wb97xd3 \
  --output-dir tddft_ensemble_decoherence/wb97xd3_100fs100

conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/audit_paired_site_qmmm.py \
  --window short100fs100_wb97xd3 --expected-frames 100 \
  --result-dir tddft_ensemble_decoherence/wb97xd3_100fs100
```

The electronic-structure method is the retained validated
ωB97X-D3/6-311G** embedded benchmark, with charge −1 and singlet multiplicity.
The same-geometry final STEOM calibration is
\(\Delta = 2.366633 - 3.43441587 = -1.06778287\) eV. This additive shift changes
absolute photon energies but cannot change the A−B gap fluctuations or
decoherence time. Couplings are not energy-shifted or empirically rescaled.

## Movie

The final movie is rendered from all 100 globally tracked measured density
states at 2 frames/s, giving exactly 50 seconds. It uses the fixed
reference-state camera, 90% transparent beta-barrel cartoon, and genuine
30%-transparent red/blue volumetric isosurfaces. Each compact physical-atom
density is rasterised onto a 0.30 Å grid, Gaussian-smoothed by 0.40 Å, and
triangulated at normalized ±0.06 using marching cubes, following the
repository's established paper-style density rendering. No coordinates or
densities are interpolated. The earlier pseudoatom-sphere movie is retained
under its original `fixed_camera_visible` name as a legacy diagnostic but is
not the primary visualisation.

```bash
/home/robson/anaconda3/envs/TeraChem/bin/python \
  tddft_ensemble_decoherence/scripts/render_tddft_density_movie.py \
  --trajectory \
  tc_tandem_nvt_100fs_100frames_20260724/tandem_nvt_100fs_100frames.dcd \
  --topology \
  tc_tandem_nvt_100fs_100frames_20260724/solvated_protonated.pdb \
  --tracking \
  tddft_ensemble_decoherence/wb97xd3_100fs100/global_state_tracking.json \
  --out \
  tddft_ensemble_decoherence/wb97xd3_100fs100/tddft_density_wb97xd3_100fs_50s_fixed_camera_isosurfaces.mp4 \
  --fps 2 --rotate-deg 0 \
  --density-representation surface \
  --surface-spacing 0.30 --surface-smear 0.40 \
  --iso-level 0.06 --density-transparency 0.30
```

## Final result

All 100 paired frames completed with zero failures and zero retries. The final
configuration-tracked cumulant T2* is 41.477 fs; the constant-shifted STEOM
value is identical. The direct signed-J propagation crossing is 63.514 fs,
versus 64.075 fs for gap-only propagation. See `NUMERICAL_REPORT.md` for the
uncertainty and interpretation and `provenance.json` for exact hashes,
software versions, GPU details, and validation tolerances.

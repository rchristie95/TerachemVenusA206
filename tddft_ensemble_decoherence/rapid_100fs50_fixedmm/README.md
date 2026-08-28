# Corrected rapid 100 fs TDDFT/TDC ensemble

This directory contains the restartable 50-frame, 2 fs-cadence prototype with
a fixed frame-0 MM atom selection. The full classical NVT system remains intact;
only the 44-atom QM single-point region is bond-capped. No frame is interpolated
or replaced.

## Exact commands

```bash
python tddft_ensemble_decoherence/scripts/initialize_manifests.py
python tddft_ensemble_decoherence/scripts/register_nvt_trajectories.py \
  --short100fs50-fixedmm-dcd \
  tc_tandem_nvt_100fs_50frames_20260723/tandem_nvt_100fs_50frames.dcd

# Frame 0 creates fixed_embedding_selection.json; frames 1--49 reuse it.
conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/prepare_production_frame.py \
  --topology tc_tandem_nvt_100fs_50frames_20260723/solvated_protonated.pdb \
  --trajectory tc_tandem_nvt_100fs_50frames_20260723/tandem_nvt_100fs_50frames.dcd \
  --frame-index FRAME --output-dir FRAME_DIR \
  --embedding-cache tddft_ensemble_decoherence/cache/embedding_charges.npz \
  --fixed-embedding-selection \
  tddft_ensemble_decoherence/rapid_100fs50_fixedmm/fixed_embedding_selection.json

conda run -n TeraChem python -u \
  tddft_ensemble_decoherence/scripts/run_rapid_adaptive_100ns50.py \
  --window short100fs50_fixedmm \
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm \
  --maximum-workers 3

conda run -n TeraChem python -u \
  tddft_ensemble_decoherence/scripts/global_state_tracking_100fs.py \
  --window short100fs50_fixedmm \
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm

conda run -n TeraChem python \
  tddft_ensemble_decoherence/scripts/summarize_100fs50.py \
  --window short100fs50_fixedmm \
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm

python tddft_ensemble_decoherence/scripts/audit_paired_site_qmmm.py
```

## Result

- Complete frames: 50/50; final failures: 0; interpolation: 0
- Paired site TDDFTs: 100/100 complete; 50 site A + 50 site B
- Per-site QM region: 29-atom CR2 + 12-atom stacked Tyr phenol + 3 link H
- Each A/B pair uses the same instantaneous full-dimer/MM snapshot
- Canonical cumulant T2*: 11.628 fs
- Positive-sequence bath correlation time: 7.583 fs
- Raw mean wavelengths: A 385.31 nm; B 370.35 nm
- Constant-STEOM-shifted mean wavelengths: A 582.78 nm; B 549.57 nm
- Screened signed TDC: mean 78.23 cm^-1,
  SD 22.02 cm^-1

See `NUMERICAL_REPORT.md`, `summary.json`, `provenance.json`, and
`convergence_failure_audit.json`. Large trajectory, raw cube, and scratch data
remain local and are not intended for version control.

`paired_site_qmmm_manifest.csv` is the explicit 100-row individual-site
restart manifest. `paired_site_qmmm_audit.json` verifies the 50 A/50 B split,
QM composition, atom mappings, link-cap directions, link-atom density
exclusions, TeraChem termination, selected compact densities, runtimes, and
output hashes. Archived pre-retry directories are preserved but are not counted
as additional retained calculations.

## PyMOL movie

`tddft_density_100fs.mp4` is a 50-frame, 1600x1200 H.264 movie of the actual
globally tracked TDDFT transition densities. It follows the established
paper-video style: black background, dim rainbow protein cartoon, highlighted
CR2 chromophores, red/blue density phase, barrel-A stabilization, and a slow
90-degree orbit. Link atoms are not displayed in the density. Exact render
settings, hashes, and sub-microangstrom placement-fit diagnostics are recorded
in `tddft_density_100fs.provenance.json`.

`tddft_density_100fs_50s.mp4` is the 50-second presentation version. Each
measured TDDFT frame is held for one second; no optical-flow or density
interpolation is used. Its encoding and source hashes are recorded in
`tddft_density_100fs_50s.provenance.json`.

`tddft_density_100fs_50s_fixed_camera.mp4` is the preferred fixed-camera
50-second presentation version. All trajectory states are fitted to frame 1
using the established barrel-A C-alpha selection, the PyMOL view matrix is then
held fixed, and each of the 50 measured TDDFT frames is shown for one second.
There is no camera orbit, optical flow, or density interpolation. Exact camera,
alignment, input-hash, and encoding details are recorded in
`tddft_density_100fs_50s_fixed_camera.provenance.json`.

`tddft_density_100fs_50s_fixed_camera_visible.mp4` retains that fixed camera
and timing while making the transition density readable through the beta
barrels. The protein cartoon is 90% transparent, the red/blue density spheres
are 2% transparent, and PyMOL transparency mode 2 is used for compositing.
This is a rendering-only change: density coordinates, values, phases, selected
states, and frame timing are unchanged. These visibility settings are now the
defaults in `scripts/render_tddft_density_movie.py` and are recorded in the
movie provenance.

# Rapid 50-frame TDDFT/TDC campaign

This campaign samples 50 intact-MM NVT snapshots at 2 ns spacing over
2--100 ns. OpenMM is completed before electronic-structure work starts, so
TeraChem and OpenMM do not contend for the single RTX 4080.

The rapid-prototype electronic-structure profile is B3LYP/6-31G* with four
TDDFT roots. This profile is intentionally lighter than the validated
wB97X-D3/6-311G** benchmark. The physical bright state is selected by spatial
transition-density overlap, not by a fixed root number. A missing or
root-boundary state is retried with seven roots; SCF/eigensolver failures are
retried with tighter thresholds. Original failed directories are retained.

The same rapid method is first run on the identical retained STEOM reference
geometry and point-charge embedding. Its energy correction is therefore
recomputed as `2.366633 eV - E_B3LYP/6-31G*(reference)`; the older correction
for wB97X-D3/6-311G** is not reused. The constant correction changes absolute
wavelengths but cancels from the A--B gap and cannot change gap-fluctuation
decoherence.

After the MM trajectory finishes, the queued driver runs:

```bash
conda run --no-capture-output -n TeraChem \
  python -u tddft_ensemble_decoherence/scripts/validate_100ns50_trajectory.py
python -u tddft_ensemble_decoherence/scripts/register_nvt_trajectories.py \
  --long100ns50-dcd \
  tc_tandem_nvt_100ns_50frames_20260723/tandem_nvt_100ns_50frames.dcd
conda run --no-capture-output -n TeraChem \
  python -u tddft_ensemble_decoherence/scripts/run_rapid_reference_calibration.py
conda run --no-capture-output -n TeraChem \
  python -u tddft_ensemble_decoherence/scripts/run_rapid_adaptive_100ns50.py \
  --maximum-workers 4
python -u tddft_ensemble_decoherence/scripts/summarize_100ns50.py
```

Concurrency is selected from measured one-, two-, three-, and four-worker
production pilots. A level is accepted only if it completes without failures,
increases throughput by at least 5%, and remains below 90% GPU-memory
occupancy. One TeraChem process per frame is active; A and B site jobs remain
serial within a frame. A campaign lock prevents duplicate schedulers, and only
the parent process writes the manifest.

The MM simulation contains no bond cuts. QM link atoms are added only during
frame extraction along the instantaneous severed-bond vectors. Link atoms
42--44 are excluded from physical transition-density maps and TDC integrals.

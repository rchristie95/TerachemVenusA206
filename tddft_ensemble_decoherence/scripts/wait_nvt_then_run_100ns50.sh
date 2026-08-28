#!/usr/bin/env bash
set -euo pipefail

nvt_pid="${1:?usage: wait_nvt_then_run_100ns50.sh NVT_PYTHON_PID}"
while kill -0 "${nvt_pid}" 2>/dev/null; do
    sleep 30
done

conda run --no-capture-output -n TeraChem \
    python -u tddft_ensemble_decoherence/scripts/validate_100ns50_trajectory.py

python -u tddft_ensemble_decoherence/scripts/register_nvt_trajectories.py \
    --long100ns50-dcd \
    tc_tandem_nvt_100ns_50frames_20260723/tandem_nvt_100ns_50frames.dcd

conda run --no-capture-output -n TeraChem \
    python -u tddft_ensemble_decoherence/scripts/run_rapid_reference_calibration.py \
    --gpu 0

conda run --no-capture-output -n TeraChem \
    python -u tddft_ensemble_decoherence/scripts/run_rapid_adaptive_100ns50.py \
    --gpu 0 --maximum-workers 4

python -u tddft_ensemble_decoherence/scripts/summarize_100ns50.py

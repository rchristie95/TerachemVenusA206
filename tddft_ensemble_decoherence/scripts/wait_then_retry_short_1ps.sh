#!/usr/bin/env bash
set -euo pipefail

campaign_pid="${1:?usage: wait_then_retry_short_1ps.sh CAMPAIGN_PID}"
while kill -0 "${campaign_pid}" 2>/dev/null; do
    sleep 20
done

conda run --no-capture-output -n TeraChem \
    python -u tddft_ensemble_decoherence/scripts/retry_failed_frames.py \
    --window short --start 0 --stop 251 --gpu 0

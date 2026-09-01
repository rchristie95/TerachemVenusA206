#!/usr/bin/env bash
# 5-frame calibration at the final campaign settings: def2-SVP, 7 roots,
# full-system embedding, pointcharges_self_interaction off.
# Measures the real per-frame cost so the frame count per session is set from a
# measurement rather than the ~215 s/site estimate.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
S=calib_v5_status.txt
echo "[$(date '+%F %H:%M:%S')] calibration start: 5 frames, def2-SVP, full embedding" > $S
t0=$(date +%s)
timeout 3h "$PY" scripts/run_ensemble.py --window v5camp --profile svp --roots 7 \
    --start 0 --stop 5 --gpu 0 >> calib_v5.log 2>&1
rc=$?
echo "[$(date '+%F %H:%M:%S')] calibration exit rc=$rc after $(( ($(date +%s)-t0)/60 )) min" >> $S

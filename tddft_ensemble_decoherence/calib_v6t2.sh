#!/usr/bin/env bash
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export OMP_NUM_THREADS=8
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
echo "[$(date '+%F %H:%M:%S')] v6t2 calibration: 2 frames of v6t2a"
time $PY -u scripts/run_ensemble.py --window v6t2a --profile svp --roots 7 --gpu 0 --start 0 --stop 2
echo "[$(date '+%F %H:%M:%S')] calibration exit rc=$?"

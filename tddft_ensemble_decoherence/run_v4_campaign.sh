#!/usr/bin/env bash
# Site-energy campaign on the CORRECTED 507-aa trajectory (window "v4short").
# GPU/TeraChem only -- deliberately does not touch the 8 ORCA ranks or the 37 GB
# the gas-phase STEOM job is holding.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
set +u
source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh
set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
WINDOW=${1:-v4full}; START=${2:-0}; STOP=${3:-25}
echo "[$(date '+%F %H:%M:%S')] $WINDOW $START..$STOP start" >> v4_campaign_status.txt
$PY scripts/run_ensemble.py --window "$WINDOW" --start "$START" --stop "$STOP" >> v4_campaign.log 2>&1
echo "[$(date '+%F %H:%M:%S')] $WINDOW $START..$STOP exit rc=$?" >> v4_campaign_status.txt

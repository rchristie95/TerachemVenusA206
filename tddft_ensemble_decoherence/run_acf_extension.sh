#!/usr/bin/env bash
# Item 6: extend the site-A QM/MM gap ACF beyond the present 592 fs.
# 149 of the 2000 registered frames of tandem_nvt_8ps_4fs.dcd are done (596 fs).
# This continues CONTIGUOUSLY from 149 -- an ACF needs consecutive 4 fs spacing,
# so the usual "shuffle frames not jobs" rule does NOT apply here.
# ~173 s/frame (both sites + TDC). Restartable: completed rows are skipped, so
# it can be stopped at any point and whatever has completed is usable.
# ONE worker only: run_ensemble.py's CLI is deliberately sequential and the
# adaptive driver is documented as the sole manifest writer under concurrency.
set -u
cd "$(dirname "$0")" || exit 1
set +u
source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh
set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
echo "[$(date '+%F %H:%M:%S')] ACF extension start, frames 149->1000" > acf_status.txt
$PY scripts/run_ensemble.py --window short --start 149 --stop 1000 >> acf_extension.log 2>&1
echo "[$(date '+%F %H:%M:%S')] ACF extension exit rc=$?" >> acf_status.txt

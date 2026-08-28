#!/usr/bin/env bash
# Rerun the 2000-frame / 4 fs window with a FIXED MM embedding selection.
#
# WHY: the original `short` campaign rebuilt the MM charge set every frame. That
# is the noise source -- not TDDFT convergence, which was the obvious suspect.
# On the matched 100 fs pair (same frames, same trajectory, only the MM policy
# differing): site-A sd 1028.6 -> 291.9 cm^-1 and root switches 13/49 -> 0/49.
# The 4 fs window was never given a fixed selection (run_ensemble.py only listed
# short100fs50_fixedmm and short100fs100_wb97xd3); it is now.
#
# Frames run IN ORDER so an early stop still leaves a contiguous prefix, which is
# what an autocorrelation needs. Resumable: completed rows are skipped.
# Old varying-MM output is archived in frames/short_varyingmm_ARCHIVE and
# manifests/short_frames_varyingmm_ARCHIVE.csv.
set -u
cd "$(dirname "$0")" || exit 1
# TeraChem must be on PATH or every frame dies with exit 127 (misreported by the
# driver as ground_state_scf_convergence). SetTCVars.sh dereferences
# LD_LIBRARY_PATH, so it needs the guard under set -u.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export PATH="/home/robson/Desktop/TeraChemPython/TeraChem/bin:$PATH"
command -v terachem >/dev/null || { echo "FATAL: terachem not on PATH" | tee -a short_fixedmm_status.txt; exit 1; }
LOG=short_fixedmm_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
note "===== short (4 fs, fixed MM) campaign start (pid $$) ====="
note "GPU $(nvidia-smi --query-gpu=name,memory.used --format=csv,noheader)"
$PY scripts/run_ensemble.py --window short --gpu 0 --roots 7 >> short_fixedmm_run.log 2>&1
rc=$?
done_n=$(awk -F, 'NR>1 && $30=="complete"' manifests/short_frames.csv | wc -l)
note "exit rc=$rc; complete frames: ${done_n}/2000"
note "===== campaign end ====="

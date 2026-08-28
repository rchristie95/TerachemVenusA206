#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export PATH="/home/robson/Desktop/TeraChemPython/TeraChem/bin:$PATH"
command -v terachem >/dev/null || { echo FATAL; exit 1; }
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a beta2_status.txt; }
note "===== constrained S1 opt, targeting the bright state ====="
terachem s1opt_c.in > s1opt_c.out 2>&1
note "opt rc=$? finished=$(grep -c 'Job finished' s1opt_c.out)"
note "===== done ====="

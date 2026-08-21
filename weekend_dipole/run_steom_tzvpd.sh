#!/usr/bin/env bash
# Single-job rerun of STEOM-DLPNO-CCSD/def2-TZVPD, same pattern as run_weekend.sh.
# The plain def2-TZVP job crashed in orca_mdci_mpi (MPI rank abort, signal 6);
# this checks whether the diffuse basis avoids it.
set -u
cd "$(dirname "$0")" || exit 1

ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH

STATUS=steom_tzvpd_status.txt
JOB=steom_def2_tzvpd
HOURS=96

note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }

clean_scratch() {
  # NOTE the glob: ORCA writes ${JOB}.<STAGE>.tmp.<rank> (e.g.
  # steom_def2_tzvpd.PAO_V12.tmp.8), NOT ${JOB}.tmp*. An earlier version of this
  # function matched only "${JOB}.tmp*" and therefore deleted nothing, leaving
  # 299 GB of PNO/PAO scratch behind. Density files are excluded by name and
  # must never be swept.
  find . -maxdepth 1 -name "${JOB}*.tmp*" ! -iname "*densit*" -delete 2>/dev/null
  find . -maxdepth 1 \( -name "${JOB}.bas[0-9]*" -o -name "${JOB}.hostnames" \) -delete 2>/dev/null
}

# systemctl stop kills this script before its cleanup runs, so also sweep on
# SIGTERM/SIGINT. Without this an interrupted run strands its whole scratch.
trap 'clean_scratch; exit 143' TERM INT

note "===== steom_def2_tzvpd start (pid $$) ====="
note "disk $(df --output=avail -BG / | tail -1 | tr -dc '0-9') GB, RAM $(free -g | awk 'NR==2{print $7}') GB"

timeout "${HOURS}h" "$ORCA/orca" "$JOB.inp" > "$JOB.out" 2>&1
rc=$?
if grep -q "ORCA TERMINATED NORMALLY" "$JOB.out"; then
  note "OK      $JOB"
elif [ $rc -eq 124 ]; then
  note "TIMEOUT $JOB after ${HOURS} h"
else
  note "FAIL    $JOB (rc=$rc) -- see $JOB.out"
fi
clean_scratch
note "===== steom_def2_tzvpd complete ====="

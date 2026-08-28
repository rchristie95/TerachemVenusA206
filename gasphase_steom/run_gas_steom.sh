#!/usr/bin/env bash
# Item 1: gas-phase DLPNO-STEOM-CCSD/def2-SVPD control for Table 1.
set -u
cd "$(dirname "$0")" || exit 1

ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH

STATUS=gas_steom_status.txt
JOB=steom_gas_svpd
HOURS=36

note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }

clean_scratch() {
  find . -maxdepth 1 -name "${JOB}*.tmp*" ! -iname "*densit*" -delete 2>/dev/null
  find . -maxdepth 1 \( -name "${JOB}.bas[0-9]*" -o -name "${JOB}.hostnames" \) -delete 2>/dev/null
}
trap 'clean_scratch; exit 143' TERM INT

note "===== $JOB start (pid $$) ====="
note "disk $(df --output=avail -BG / | tail -1 | tr -dc '0-9') GB, RAM $(free -g | awk 'NR==2{print $7}') GB free"

timeout "${HOURS}h" "$ORCA/orca" "$JOB.inp" > "$JOB.out" 2>&1
rc=$?
if grep -q "ORCA TERMINATED NORMALLY" "$JOB.out"; then
  note "OK      $JOB"
  note "$(grep -E '^IROOT= +1:' "$JOB.out" | tail -1)"
elif [ $rc -eq 124 ]; then
  note "TIMEOUT $JOB after ${HOURS} h"
else
  note "FAIL    $JOB (rc=$rc)"
  note "$(grep -iE 'error|not enough memory|aborting' "$JOB.out" | tail -3)"
fi
clean_scratch
note "===== $JOB complete ====="

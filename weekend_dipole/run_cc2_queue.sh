#!/usr/bin/env bash
# CC2 transition-dipole queue (see cc2_def2_svp.inp header for why not adcc).
set -u
cd "$(dirname "$0")" || exit 1
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
LOG=cc2_queue_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }
clean(){ find . -maxdepth 1 -name "$1*.tmp*" ! -iname '*densit*' -delete 2>/dev/null
         find . -maxdepth 1 \( -name "$1.bas[0-9]*" -o -name "$1.hostnames" \) -delete 2>/dev/null; }
trap 'clean cc2_def2_svp; clean cc2_def2_svpd; exit 143' TERM INT
note "===== CC2 dipole queue start (pid $$) ====="
for J in cc2_def2_svp cc2_def2_svpd; do
  if grep -q 'ORCA TERMINATED NORMALLY' "$J.out" 2>/dev/null; then note "SKIP    $J"; continue; fi
  note "START   $J"
  timeout 20h "$ORCA/orca" "$J.inp" > "$J.out" 2>&1; rc=$?
  if grep -q 'ORCA TERMINATED NORMALLY' "$J.out"; then note "OK      $J"
  elif [ $rc -eq 124 ]; then note "TIMEOUT $J"
  else note "FAIL    $J (rc=$rc)"; fi
  clean "$J"
done
note "===== CC2 queue complete ====="

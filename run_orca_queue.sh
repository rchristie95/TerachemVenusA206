#!/usr/bin/env bash
# SEQUENTIAL ORCA queue. Two ORCA jobs were OOM-killed on 2026-08-27 by being run
# concurrently (48 GB + 24 GB + MD on a 62 GB box). Nothing here runs in parallel.
# maxcore is set to the 5000 x 8 = 40 GB that the embedded def2-SVPD STEOM run
# actually completed with, leaving ~22 GB of headroom.
set -u
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
Q=/home/robson/PetaChem/orca_queue_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$Q"; }

run(){ # dir  inputbase
  local d=$1 j=$2
  cd "$d" || return 1
  note "START $j (free $(free -g | awk 'NR==2{print $7}') GB)"
  timeout 36h "$ORCA/orca" "$j.inp" > "$j.out" 2>&1
  if grep -q "ORCA TERMINATED NORMALLY" "$j.out"; then note "OK    $j"
  else note "FAIL  $j :: $(grep -iE 'error|not enough' "$j.out" | tail -1)"; fi
  find . -maxdepth 1 -name "${j}*.tmp*" ! -iname "*densit*" -delete 2>/dev/null
}

note "===== queue start ====="
run /home/robson/PetaChem/gasphase_steom steom_gas_svpd
run /home/robson/PetaChem/region_control  tddft_44
run /home/robson/PetaChem/region_control  tddft_47
note "===== queue complete ====="

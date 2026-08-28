#!/usr/bin/env bash
set -u
cd "$(dirname "$0")" || exit 1
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
echo "[$(date '+%F %H:%M:%S')] region control start" > region_status.txt
for N in 44 47; do
  $ORCA/orca tddft_$N.inp > tddft_$N.out 2>&1
  ok=$(grep -c "ORCA TERMINATED NORMALLY" tddft_$N.out)
  echo "[$(date '+%F %H:%M:%S')] tddft_$N done, normal=$ok" >> region_status.txt
done
echo "[$(date '+%F %H:%M:%S')] region control complete" >> region_status.txt

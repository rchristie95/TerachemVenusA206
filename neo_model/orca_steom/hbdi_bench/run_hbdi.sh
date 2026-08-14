#!/usr/bin/env bash
# HBDI comparison ladder. Resumable: skips any job already terminated normally.
set -u; cd "$(dirname "$0")" || exit 1
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
for j in tddft_camb3lyp tddft_wb97xd3 steom_hbdi; do
  if [ -f $j.out ] && grep -q "ORCA TERMINATED NORMALLY" $j.out; then
    echo "[$(date +%H:%M:%S)] SKIP  $j"; continue
  fi
  echo "[$(date +%H:%M:%S)] START $j"
  "$ORCA/orca" $j.inp > $j.out 2>&1
  grep -q "ORCA TERMINATED NORMALLY" $j.out && echo "[$(date +%H:%M:%S)] OK    $j" || echo "[$(date +%H:%M:%S)] FAIL  $j"
done
echo "[$(date +%H:%M:%S)] hbdi ladder complete"

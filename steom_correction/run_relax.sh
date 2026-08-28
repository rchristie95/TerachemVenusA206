#!/usr/bin/env bash
set -u
cd /home/robson/PetaChem/steom_correction || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
echo "[$(date '+%F %H:%M:%S')] relax start" > relax_status.txt
for m in wb97xd3 camb3lyp; do
  /home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem relax_$m.in > relax_$m.out 2>&1
  echo "[$(date '+%F %H:%M:%S')] $m rc=$?" >> relax_status.txt
done
echo "[$(date '+%F %H:%M:%S')] relax done" >> relax_status.txt

#!/usr/bin/env bash
set -u
cd /home/robson/PetaChem/steom_correction || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
echo "[$(date '+%F %H:%M:%S')] start" > status.txt
/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem tddft_def2svpd.in > tddft_def2svpd.out 2>&1
echo "[$(date '+%F %H:%M:%S')] exit rc=$?" >> status.txt

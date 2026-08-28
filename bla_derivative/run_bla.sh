#!/usr/bin/env bash
set -u
cd /home/robson/PetaChem/bla_derivative || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
echo "[$(date '+%F %H:%M:%S')] BLA derivative start" > bla_status.txt
for g in ref_q0 mode69_q-0.5 mode69_q+0.5 mode77_q-0.5 mode77_q+0.5 mode86_q-0.5 mode86_q+0.5 mode101_q-0.5 mode101_q+0.5; do
  /home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem grad_$g.in > grad_$g.out 2>&1
  echo "[$(date '+%F %H:%M:%S')] $g rc=$?" >> bla_status.txt
done
echo "[$(date '+%F %H:%M:%S')] BLA derivative done" >> bla_status.txt

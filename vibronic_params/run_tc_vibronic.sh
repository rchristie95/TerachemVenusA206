#!/usr/bin/env bash
# GPU vibronic parameters: opt -> numerical Hessian -> excited-state gradient.
# Huang-Rhys factors are computed afterwards from the Hessian + ES gradient
# (TeraChem has no ESD-equivalent module), as S_i = g_i^2 / (2 w_i^3).
set -u
cd "$(dirname "$0")/tc" || exit 1
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"   # SetTCVars.sh dereferences it under set -u
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
TC=/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem
LOG=../tc_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }
# 'FINAL ENERGY' appears in EVERY SCF, so it cannot mark job completion -- a
# killed run looked finished and the pipeline skipped straight to a missing file.
ok(){ grep -q 'Job finished' "$1" 2>/dev/null; }

note "===== TeraChem vibronic run start (pid $$) ====="
note "GPU $(nvidia-smi --query-gpu=name,memory.used --format=csv,noheader)"

if ! ok opt.out; then
  note "START   stage 1: GPU optimisation (from ORCA cycle-69 geometry)"
  timeout 6h "$TC" opt.in > opt.out 2>&1
  ok opt.out && note "OK      stage 1" || { note "FAIL    stage 1"; exit 1; }
fi
# TeraChem writes the optimised geometry under scr_opt/
# optim_geom.xyz is only written on SUCCESS; optim_trajectory.xyz accumulates
# every cycle and is the reliable source. Never accept an empty file -- an empty
# optim_geom.xyz silently clobbered geom.xyz once and killed the run instantly.
for c in scr_opt/optim_geom.xyz scr_opt/optim_trajectory.xyz; do
  [ -s "$c" ] || continue
  python3 - "$c" <<'EX'
import sys
L=open(sys.argv[1]).read().rstrip().split('\n'); n=int(L[0].split()[0]); nb=n+2
nf=len(L)//nb
open('optim.xyz','w').write('\n'.join(L[nb*(nf-1):nb*nf])+'\n')
EX
  [ -s optim.xyz ] && { note "optimised geometry <- $c (last frame)"; break; }
done
[ -f optim.xyz ] || { note "FAIL: no optimised geometry found"; exit 1; }

if ! ok freq.out; then
  note "START   stage 2: numerical Hessian (264 displacements)"
  timeout 10h "$TC" freq.in > freq.out 2>&1
  ok freq.out && note "OK      stage 2" || note "FAIL    stage 2"
fi
if ! ok esgrad.out; then
  note "START   stage 3: excited-state gradient at the FC point"
  timeout 4h "$TC" esgrad.in > esgrad.out 2>&1
  ok esgrad.out && note "OK      stage 3" || note "FAIL    stage 3"
fi
note "===== complete -- run huang_rhys.py to extract S_i, omega_i ====="

#!/usr/bin/env bash
# 50-frame FULL-SYSTEM embedding ensemble, spread across the whole 8 ps
# trajectory (stride 40 frames = 160 fs), frames 50..1999.
#
# WHY SPREAD, NOT CONTIGUOUS. The existing v4full ensemble is 46 usable frames
# at indices 0-49, i.e. only the first 196 fs. That window can see FAST motion
# only. It already answered the primary question -- full-system embedding gives
# <|Delta|> = 574.0 +/- 53.9 cm^-1 against the 12 A production value of 576, so
# the detuning is real and not a truncation artefact. What it cannot give is the
# SLOW/static component, which is the quantity the exciton actually sees
# (tau_c = 8.5 fs vs hbar/2J = 81 fs). Sampling 8 ps instead of 196 fs is the
# missing measurement.
#
# NO OVERALL TIMEOUT, by instruction. A per-frame guard of 2 h is retained: it
# is 8.5x the measured 14.2 min/frame, so it cannot plausibly kill good work,
# but it stops one hung frame from consuming the entire 12 h budget. Remove
# PER_FRAME_GUARD below if that is not wanted.
set -u
cd "$(dirname "$0")" || exit 1
# TeraChem needs its env vars sourced and its own conda python (the system
# python3 has no scipy, which is how the first launch died in 90 s).
set +u
source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh
set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
PER_FRAME_GUARD=2h
STATUS=spread50_status.txt
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }

note "===== spread-50 full-system ensemble start (pid $$) ====="
note "frames: $(wc -l < frames_spread50.txt), stride 40 (160 fs), span 0.2-8.0 ps"
note "no overall timeout; per-frame guard $PER_FRAME_GUARD"

done_n=0; fail_n=0; skip_n=0; i=0; total=$(wc -l < frames_spread50.txt)
t0=$(date +%s)
while read -r f; do
  i=$((i+1))
  summary=$(printf "frames/v4full/frame_%04d/frame_summary.json" "$f")
  if [ -f "$summary" ]; then skip_n=$((skip_n+1)); note "[$i/$total] frame $f already done, skipping"; continue; fi
  fs=$(date +%s)
  timeout "$PER_FRAME_GUARD" "$PY" scripts/run_ensemble.py \
      --window v4full --start "$f" --stop $((f+1)) --gpu 0 >> spread50_run.log 2>&1
  rc=$?
  dt=$(( $(date +%s) - fs ))
  if [ -f "$summary" ]; then
    done_n=$((done_n+1)); note "[$i/$total] frame $f OK (${dt}s)"
  else
    fail_n=$((fail_n+1)); note "[$i/$total] frame $f FAILED rc=$rc (${dt}s)"
  fi
  el=$(( $(date +%s) - t0 ))
  if [ $((i % 5)) -eq 0 ]; then
    note "  progress: $done_n ok, $fail_n failed, $skip_n skipped; elapsed $((el/60)) min; projected total $(( el * total / i / 60 )) min"
  fi
done < frames_spread50.txt

note "===== done: $done_n ok, $fail_n failed, $skip_n skipped, $(( ($(date +%s)-t0)/60 )) min ====="
"$PY" - <<'PY' 2>&1 | tee -a "$STATUS"
import json, glob, statistics as st
rows=[]
for p in sorted(glob.glob('frames/v4full/frame_*/frame_summary.json')):
    d=json.load(open(p)); A,B=d.get('site_A'),d.get('site_B')
    if A and B and A.get('status')=='selected' and B.get('status')=='selected':
        rows.append((d['frame_index'], abs(A['energy_cm-1']-B['energy_cm-1'])))
g=[x for _,x in rows]
early=[x for i,x in rows if i<=49]; late=[x for i,x in rows if i>49]
print(f"\nFULL-SYSTEM ensemble, n={len(g)} frames")
print(f"  all      <|Delta|> = {st.mean(g):8.1f} +/- {st.stdev(g)/len(g)**0.5:5.1f} cm^-1 (sd {st.stdev(g):.1f})")
if early: print(f"  0-196 fs <|Delta|> = {st.mean(early):8.1f} cm^-1  (n={len(early)}, fast motion only)")
if len(late)>1: print(f"  0.2-8 ps <|Delta|> = {st.mean(late):8.1f} +/- {st.stdev(late)/len(late)**0.5:5.1f} cm^-1 (n={len(late)}, spans slow motion)")
print("  reference: 12 A production ensemble reported 576 cm^-1")
PY

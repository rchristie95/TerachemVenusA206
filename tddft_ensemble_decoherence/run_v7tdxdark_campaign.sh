#!/usr/bin/env bash
# dVenus-TDX barrel-2 supervisor: the second QM region, giving Delta and J.
#
# Three workers over interleaved slices, same shape as run_v7tdx_campaign.sh.
# No STEOM stream: the dark region has no chromophore to calibrate.
# Site A is never touched -- this writes only site_Bdark/ and tdc_dark/ inside
# the existing frame directories and reads site A's stored root densities.
set -u
cd "$(dirname "$0")" || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
export OMP_NUM_THREADS=4
LOGS=logs/v7tdxdark; mkdir -p "$LOGS"
STATUS="$LOGS/supervisor.log"; INCIDENTS="$LOGS/incidents.log"
WORKERS=(a b c); MAX_RESTARTS=6; STALL_RESTARTS=3
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }
count() { "$PY" -c 'import csv,sys;print(sum(r["convergence_status"]==sys.argv[2] for r in csv.DictReader(open(sys.argv[1]))))' "manifests/v7tdxdark_$1_frames.csv" "$2" 2>/dev/null || echo 0; }
completed(){ count "$1" complete; }
failed_rows(){ count "$1" failed; }
total_rows(){ echo $(( $(wc -l < "manifests/v7tdxdark_$1_frames.csv") - 1 )); }
start_worker(){ "$PY" -u scripts/run_tdx_darksite.py --manifest-suffix "_$1" --gpu 0 >> "$LOGS/worker_$1.log" 2>&1 & echo $!; }
declare -A PID RESTARTS STALLED LASTDONE RETIRED
for w in "${WORKERS[@]}"; do RESTARTS[$w]=0; STALLED[$w]=0; RETIRED[$w]=0; LASTDONE[$w]=$(completed "$w"); done
cleanup(){ note "SIGTERM -- stopping children"; for w in "${WORKERS[@]}"; do [ -n "${PID[$w]:-}" ] && kill "${PID[$w]}" 2>/dev/null; done; exit 143; }
trap cleanup TERM INT
note "===== TDX barrel-2 supervisor start (pid $$) ====="
for w in "${WORKERS[@]}"; do PID[$w]=$(start_worker "$w"); note "worker $w started (pid ${PID[$w]}, $(completed "$w")/$(total_rows "$w"))"; sleep 20; done
while :; do
  sleep 120; 
  for w in "${WORKERS[@]}"; do
    [ "${RETIRED[$w]}" -eq 1 ] && continue
    kill -0 "${PID[$w]}" 2>/dev/null && continue
    d=$(completed "$w"); t=$(total_rows "$w")
    if [ "$d" -ge "$t" ]; then note "worker $w finished: $d/$t"; RETIRED[$w]=1; continue; fi
    if [ "$d" -gt "${LASTDONE[$w]}" ]; then STALLED[$w]=0; else STALLED[$w]=$(( STALLED[$w]+1 )); fi
    LASTDONE[$w]=$d; RESTARTS[$w]=$(( RESTARTS[$w]+1 ))
    if [ "${STALLED[$w]}" -ge "$STALL_RESTARTS" ] || [ "${RESTARTS[$w]}" -gt "$MAX_RESTARTS" ]; then
      { echo "[$(date '+%F %H:%M:%S')] RETIRING worker $w at $d/$t"; tail -20 "$LOGS/worker_$w.log" | sed 's/^/    /'; } >> "$INCIDENTS"
      RETIRED[$w]=1; continue
    fi
    note "restarting worker $w (restart ${RESTARTS[$w]}, $d/$t)"; PID[$w]=$(start_worker "$w")
  done
  r=0; for w in "${WORKERS[@]}"; do [ "${RETIRED[$w]}" -eq 0 ] && r=1; done; [ "$r" -eq 0 ] && break
done
note "main pass over: a=$(completed a) b=$(completed b) c=$(completed c) (failed: $(failed_rows a)/$(failed_rows b)/$(failed_rows c))"
note "===== TDX barrel-2 complete: $(( $(completed a)+$(completed b)+$(completed c) ))/500 ====="
exit 0

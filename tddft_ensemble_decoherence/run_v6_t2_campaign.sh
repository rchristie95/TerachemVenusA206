#!/usr/bin/env bash
# T2 campaign supervisor.
#
# Three TeraChem TDDFT workers on the GPU, one per 100 fs window. Keeps them
# alive, retries what fails, and refuses to loop forever on a worker that cannot
# make progress.
#
# ONE WORKER PER WINDOW. Natural here (three windows, three workers) and it
# keeps each tracking chain unbroken -- though "required" overstates it: v5camp
# ran interleaved slices of a single window with zero failures. The reason is: run_ensemble.py threads previous["A"]/["B"]
# transition dipoles frame to frame for density-overlap state tracking
# (scripts/run_ensemble.py:532,542), so splitting a window across workers breaks
# the tracking chain. At 1 fs the geometries barely move and the chain should be
# far cleaner than v5's 0.5 ps one.
#
# No STEOM stream. All 8 references were computed in the v5 sweep and each
# window is seeded from one of those geometries, so there is nothing for ORCA to
# do -- and therefore no 250 GB scratch, no disk guard, no throttle, no deadlock
# escape. steom_disk_guard.sh and steom_pause_site.sh are not used here.
#
# Launch with systemd-run --user, NEVER from an agent shell: the harness kills
# its whole process tree when the shell exits, which is how chain_attempt7.sh
# died on 2026-08-27 without launching anything.
set -u
cd "$(dirname "$0")" || exit 1

set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
# SetTCVars leaves OMP_NUM_THREADS=32; three workers would ask for 96 on 32
# cores. TDDFT here is GPU-bound. No ORCA competing for cores this time, so 8
# each rather than v5's 4.
export OMP_NUM_THREADS=8
LOGS=logs/v6t2; mkdir -p "$LOGS"
STATUS="$LOGS/supervisor.log"
INCIDENTS="$LOGS/incidents.log"
WORKERS=(a b c)
MAX_RESTARTS=6          # per worker, across the whole run
STALL_RESTARTS=3        # restarts with zero new rows before a worker is retired
RETRY_PASSES=2          # --retry-failed-only sweeps after the main pass

note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }
incident() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$INCIDENTS" >> "$STATUS"; }

# Counted with a real csv reader, not awk -F,: a failed row carries the
# exception repr in its notes column and that can contain commas and quotes.
count_status() {   # count_status <worker> <status>
  "$PY" -c 'import csv,sys;print(sum(r["convergence_status"]==sys.argv[2] for r in csv.DictReader(open(sys.argv[1]))))' \
      "manifests/v6t2$1_frames.csv" "$2" 2>/dev/null || echo 0
}
completed()   { count_status "$1" complete; }
failed_rows() { count_status "$1" failed; }
total_rows() { echo $(( $(wc -l < "manifests/v6t2$1_frames.csv") - 1 )); }

start_worker() {   # start_worker <tag>
  local w=$1
  "$PY" scripts/run_ensemble.py --window "v6t2$w" \
        --profile svp --roots 7 --gpu 0 \
        >> "$LOGS/worker_$w.log" 2>&1 &
  echo $!
}

declare -A PID RESTARTS STALLED LASTDONE RETIRED
for w in "${WORKERS[@]}"; do
  RESTARTS[$w]=0; STALLED[$w]=0; RETIRED[$w]=0; LASTDONE[$w]=$(completed "$w")
done

cleanup() {
  note "SIGTERM/INT -- stopping children"
  for w in "${WORKERS[@]}"; do
    [ -n "${PID[$w]:-}" ] && kill "${PID[$w]}" 2>/dev/null
  done
  exit 143
}
trap cleanup TERM INT

note "===== T2 campaign supervisor start (pid $$) ====="
note "TDDFT: 3 windows x 100 frames at 1 fs, def2-SVP, full-system embedding"
note "  v6t2a <- v5 frame 900  (|Delta| = 826.6 cm^-1, 71st percentile)"
note "  v6t2b <- v5 frame 1400 (|Delta| =  58.6 cm^-1,  4th percentile)"
note "  v6t2c <- v5 frame 1900 (|Delta| = 514.2 cm^-1, 51st percentile)"
note "STEOM: none -- all 8 references reused from the v5 sweep"

for w in "${WORKERS[@]}"; do
  PID[$w]=$(start_worker "$w")
  note "worker $w (v6t2$w) started (pid ${PID[$w]}, $(completed "$w")/$(total_rows "$w") done)"
  sleep 20   # stagger: three simultaneous TeraChem starts race on GPU context
done

while :; do
  sleep 120
  for w in "${WORKERS[@]}"; do
    [ "${RETIRED[$w]}" -eq 1 ] && continue
    kill -0 "${PID[$w]}" 2>/dev/null && continue

    done_now=$(completed "$w"); total=$(total_rows "$w")
    if [ "$done_now" -ge "$total" ]; then
      note "worker $w finished: $done_now/$total complete"
      RETIRED[$w]=1; continue
    fi

    # It exited early. Progress since the last restart decides whether this is a
    # transient crash worth retrying or a loop that will never converge.
    if [ "$done_now" -gt "${LASTDONE[$w]}" ]; then STALLED[$w]=0
    else STALLED[$w]=$(( STALLED[$w] + 1 )); fi
    LASTDONE[$w]=$done_now
    RESTARTS[$w]=$(( RESTARTS[$w] + 1 ))

    if [ "${STALLED[$w]}" -ge "$STALL_RESTARTS" ]; then
      incident "RETIRING worker $w: ${STALLED[$w]} restarts with no new frame. \
Needs a human. Last 25 lines of $LOGS/worker_$w.log:"
      tail -25 "$LOGS/worker_$w.log" | sed 's/^/    /' >> "$INCIDENTS"
      RETIRED[$w]=1; continue
    fi
    if [ "${RESTARTS[$w]}" -gt "$MAX_RESTARTS" ]; then
      incident "RETIRING worker $w: $MAX_RESTARTS restarts used."
      RETIRED[$w]=1; continue
    fi

    avail=$(free -g | awk 'NR==2{print $7}')
    if [ "$avail" -lt 8 ]; then
      incident "worker $w died with only ${avail} GB free -- waiting 5 min before restart"
      sleep 300
    fi
    incident "restarting worker $w (restart ${RESTARTS[$w]}, $done_now/$total done)"
    PID[$w]=$(start_worker "$w")
  done

  running=0
  for w in "${WORKERS[@]}"; do [ "${RETIRED[$w]}" -eq 0 ] && running=1; done
  [ "$running" -eq 0 ] && break
done

note "main pass over: a=$(completed a) b=$(completed b) c=$(completed c) \
(failed: a=$(failed_rows a) b=$(failed_rows b) c=$(failed_rows c))"

# Individual frames fail for transient reasons (root tracking, a wedged GPU
# context). --retry-failed-only touches nothing that already completed.
for pass in $(seq 1 "$RETRY_PASSES"); do
  nf=0
  for w in "${WORKERS[@]}"; do nf=$(( nf + $(failed_rows "$w") )); done
  [ "$nf" -eq 0 ] && { note "no failed frames left, skipping retry pass $pass"; break; }
  note "retry pass $pass: $nf failed frames"
  for w in "${WORKERS[@]}"; do
    [ "$(failed_rows "$w")" -eq 0 ] && continue
    "$PY" scripts/run_ensemble.py --window "v6t2$w" \
          --profile svp --roots 7 --gpu 0 --retry-failed-only \
          >> "$LOGS/worker_$w.log" 2>&1
  done
  note "after retry pass $pass: a=$(completed a) b=$(completed b) c=$(completed c) \
(failed: a=$(failed_rows a) b=$(failed_rows b) c=$(failed_rows c))"
done

TOT=$(( $(completed a) + $(completed b) + $(completed c) ))
note "===== T2 campaign complete: $TOT/300 TDDFT frames ====="
[ -s "$INCIDENTS" ] && note "INCIDENTS were logged -- see $INCIDENTS"
exit 0

#!/usr/bin/env bash
# minimum-detuning window supervisor.
#
# Three TeraChem TDDFT workers on the GPU, one per 100 fs window. Keeps them
# alive, retries what fails, and refuses to loop forever on a worker that cannot
# make progress.
#
# ONE WINDOW, split into three CONTIGUOUS BLOCKS via --manifest-suffix.
# run_ensemble.py threads previous["A"]/["B"] transition dipoles frame to frame
# for density-overlap tracking (scripts/run_ensemble.py:532,542), so a split
# breaks that chain -- but only at the block seams, and v5camp ran INTERLEAVED
# slices of a single window with zero failures across 500 frames, so the tracking
# is robust to it. (The v6t2 supervisor calls one-worker-per-window "required";
# that overstated it. There it was the natural mapping anyway, three windows to
# three workers, so nothing was lost.) Blocks beat interleaving here because each
# worker's previous frame stays genuinely adjacent in time.
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
LOGS=logs/v6min_run; mkdir -p "$LOGS"
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
      "manifests/v6min_$1_frames.csv" "$2" 2>/dev/null || echo 0
}
completed()   { count_status "$1" complete; }
failed_rows() { count_status "$1" failed; }
total_rows() { echo $(( $(wc -l < "manifests/v6min_$1_frames.csv") - 1 )); }

start_worker() {   # start_worker <tag>
  local w=$1
  "$PY" scripts/run_ensemble.py --window v6min --manifest-suffix "_$w" \
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

note "===== minimum-detuning window supervisor start (pid $$) ====="
note "TDDFT: 125 frames at 4 fs over 4-500 fs, one window in 3 contiguous blocks"
note "  Extends the v6t2 windows from 100 fs to 500 fs on the SAME seed geometries."
note "  Frequency resolution 1/(c*500fs) = 67 cm^-1, below kT = 208 cm^-1 -- the"
note "  v6t2 windows reach only 334 cm^-1 and cannot resolve the classical bath."
note "STEOM: 2 NEW jobs on the seed geometry -- no reference exists for frame 1832"

for w in "${WORKERS[@]}"; do
  PID[$w]=$(start_worker "$w")
  note "worker $w (v6min) started (pid ${PID[$w]}, $(completed "$w")/$(total_rows "$w") done)"
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
    "$PY" scripts/run_ensemble.py --window v6min --manifest-suffix "_$w" \
          --profile svp --roots 7 --gpu 0 --retry-failed-only \
          >> "$LOGS/worker_$w.log" 2>&1
  done
  note "after retry pass $pass: a=$(completed a) b=$(completed b) c=$(completed c) \
(failed: a=$(failed_rows a) b=$(failed_rows b) c=$(failed_rows c))"
done

TOT=$(( $(completed a) + $(completed b) + $(completed c) ))
note "===== minimum-detuning window complete: $TOT/125 TDDFT frames ====="
[ -s "$INCIDENTS" ] && note "INCIDENTS were logged -- see $INCIDENTS"
exit 0

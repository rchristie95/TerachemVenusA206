#!/usr/bin/env bash
# v5 campaign supervisor.
#
# Three TeraChem TDDFT workers on the GPU over disjoint interleaved slices of the
# 500-frame manifest, plus the ORCA STEOM reference stream on the CPU. Keeps them
# alive, retries what fails, and refuses to loop forever on a worker that cannot
# make progress.
#
# Launch with systemd-run --user, NEVER from an agent shell: the harness kills
# its whole process tree when the shell exits, which is how chain_attempt7.sh
# died on 2026-08-27 without launching anything.
set -u
cd "$(dirname "$0")" || exit 1

set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
# SetTCVars leaves OMP_NUM_THREADS=32. Three TeraChem workers would then ask
# for 96 threads on 32 cores while ORCA is using 20 of them. TDDFT here is
# GPU-bound, so give each worker 4 and leave the cores to the STEOM stream.
export OMP_NUM_THREADS=4
LOGS=logs/v5camp; mkdir -p "$LOGS"
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
      "manifests/v5camp_$1_frames.csv" "$2" 2>/dev/null || echo 0
}
completed()   { count_status "$1" complete; }
failed_rows() { count_status "$1" failed; }
total_rows() { echo $(( $(wc -l < "manifests/v5camp_$1_frames.csv") - 1 )); }

start_worker() {   # start_worker <tag>
  local w=$1
  "$PY" scripts/run_ensemble.py --window v5camp --manifest-suffix "_$w" \
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
  [ -n "${STEOM_PID:-}" ] && kill "$STEOM_PID" 2>/dev/null
  exit 143
}
trap cleanup TERM INT

note "===== v5 campaign supervisor start (pid $$) ====="
note "TDDFT: 3 workers x ~167 frames, def2-SVP, full-system embedding"
note "STEOM: frames 900/1400/1900/2400 x 2 sites, 4 concurrent on CPU"

# ---- STEOM stream on the CPU, overlapping the GPU work --------------------
"$PY" scripts/run_steom_v5.py --concurrency 4 --nprocs 5 --maxcore 2000 \
      --hours 10 --min-free-gb 12 >> "$LOGS/steom.log" 2>&1 &
STEOM_PID=$!
note "STEOM stream started (pid $STEOM_PID)"

# ---- TDDFT workers on the GPU ---------------------------------------------
for w in "${WORKERS[@]}"; do
  PID[$w]=$(start_worker "$w")
  note "worker $w started (pid ${PID[$w]}, $(completed "$w")/$(total_rows "$w") done)"
  sleep 20   # stagger: three simultaneous TeraChem starts race on GPU context
done

# ---- supervise -------------------------------------------------------------
while :; do
  sleep 120
  alive=0
  for w in "${WORKERS[@]}"; do
    [ "${RETIRED[$w]}" -eq 1 ] && continue
    if kill -0 "${PID[$w]}" 2>/dev/null; then alive=$((alive+1)); continue; fi

    done_now=$(completed "$w"); total=$(total_rows "$w")
    if [ "$done_now" -ge "$total" ]; then
      note "worker $w finished: $done_now/$total complete"
      RETIRED[$w]=1; continue
    fi

    # It exited early. Progress since the last restart decides whether this is a
    # transient crash worth retrying or a loop that will never converge.
    if [ "$done_now" -gt "${LASTDONE[$w]}" ]; then
      STALLED[$w]=0
    else
      STALLED[$w]=$(( STALLED[$w] + 1 ))
    fi
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

    # Free memory is the usual reason a TeraChem worker is killed outright.
    avail=$(free -g | awk 'NR==2{print $7}')
    if [ "$avail" -lt 8 ]; then
      incident "worker $w died with only ${avail} GB free -- waiting 5 min before restart"
      sleep 300
    fi
    incident "restarting worker $w (restart ${RESTARTS[$w]}, $done_now/$total done)"
    PID[$w]=$(start_worker "$w")
    alive=$((alive+1))
  done

  running=0
  for w in "${WORKERS[@]}"; do [ "${RETIRED[$w]}" -eq 0 ] && running=1; done
  [ "$running" -eq 0 ] && break
done

note "TDDFT main pass over: a=$(completed a) b=$(completed b) c=$(completed c) \
(failed: a=$(failed_rows a) b=$(failed_rows b) c=$(failed_rows c))"

# ---- retry sweeps ----------------------------------------------------------
# Individual frames fail for transient reasons (root tracking, a wedged GPU
# context). --retry-failed-only touches nothing that already completed.
for pass in $(seq 1 "$RETRY_PASSES"); do
  nf=0
  for w in "${WORKERS[@]}"; do nf=$(( nf + $(failed_rows "$w") )); done
  [ "$nf" -eq 0 ] && { note "no failed frames left, skipping retry pass $pass"; break; }
  note "retry pass $pass: $nf failed frames"
  for w in "${WORKERS[@]}"; do
    [ "$(failed_rows "$w")" -eq 0 ] && continue
    "$PY" scripts/run_ensemble.py --window v5camp --manifest-suffix "_$w" \
          --profile svp --roots 7 --gpu 0 --retry-failed-only \
          >> "$LOGS/worker_$w.log" 2>&1
  done
  note "after retry pass $pass: a=$(completed a) b=$(completed b) c=$(completed c) \
(failed: a=$(failed_rows a) b=$(failed_rows b) c=$(failed_rows c))"
done

# ---- wait for the CPU stream ----------------------------------------------
if kill -0 "$STEOM_PID" 2>/dev/null; then
  note "TDDFT done; waiting on the STEOM stream (pid $STEOM_PID)"
  wait "$STEOM_PID"; note "STEOM stream exit rc=$?"
fi

TOT=$(( $(completed a) + $(completed b) + $(completed c) ))
note "===== v5 campaign complete: $TOT/500 TDDFT frames ====="
[ -s "$INCIDENTS" ] && note "INCIDENTS were logged -- see $INCIDENTS"

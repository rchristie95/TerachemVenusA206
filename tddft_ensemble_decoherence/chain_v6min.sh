#!/usr/bin/env bash
# Queued behind the v6lf (500 fs) campaign. Two jobs, in order:
#
#   1. HOLE CLOSE. check_lf_overlap.py showed the 500 fs MD does NOT reproduce
#      the 100 fs run -- OpenMM CUDA is not bit-reproducible across runs even
#      with identical seed frame, velocity seed and integrator seed. The v6lf
#      manifests had been built to skip DCD frames 0-24 on the assumption they
#      would duplicate v6t2, so each v6lf trajectory has a genuine 0-100 fs
#      hole. 75 frames (25 x 3 windows), ~45 min.
#
#   2. MINIMUM-DETUNING WINDOW. Same 500 fs / 4 fs protocol seeded on v5 frame
#      1832, |Delta| = 0.170 cm^-1 against 2|J| = 49.1 -- the most degenerate
#      configuration in the 500-frame ensemble, and the one the manuscript
#      already features as its extreme delocalised case. Unlike every other
#      window this seed has NO existing STEOM reference, so two new ORCA jobs
#      run on the CPU alongside the GPU stream.
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
LOG=logs/v6min; mkdir -p "$LOG"
S="$LOG/chain.log"
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$S"; }
cnt() { $PY -c "import csv,sys;print(sum(sum(r['convergence_status']=='complete' for r in csv.DictReader(open(f'manifests/{sys.argv[1]}{w}_frames.csv'))) for w in 'abc'))" "$1" 2>/dev/null || echo -1; }

note "===== queued: waiting on the v6lf campaign ====="
while systemctl --user is-active --quiet v6lf-chain.service; do sleep 120; done
N=$(cnt v6lf); note "v6lf finished at $N/300"
[ "$N" -lt 300 ] && { note "ABORT: v6lf incomplete"; exit 1; }

# ---- 1. hole close ---------------------------------------------------------
note "--- job 1: closing the 0-100 fs hole in the v6lf windows ---"
if ! $PY scripts/initialize_v6lf_manifests.py --first-frame 0 >> "$S" 2>&1; then
  note "ABORT: could not extend the v6lf manifests"; exit 1
fi
note "v6lf manifests extended to 125 rows; rerunning the supervisor (completed rows are skipped)"
./run_v6_lf_campaign.sh >> "$LOG/holeclose.log" 2>&1
M=$(cnt v6lf); note "hole close done: $M/375 v6lf frames complete"

# ---- 2. minimum-detuning window -------------------------------------------
note "--- job 2: minimum-detuning window, v5 frame 1832 ---"
note "building 500 fs MD"
if ! bash ../tandem_v6_t2_campaign/run_md_ext_named.sh 1832 -1208659.943436 20260904 min \
        >> "$LOG/md_w1832.log" 2>&1; then
  note "ABORT: MD failed for frame 1832"; exit 1
fi
if ! $PY scripts/initialize_v6min_manifests.py >> "$S" 2>&1; then
  note "ABORT: v6min manifest initialisation failed"; exit 1
fi

# STEOM for the seed geometry, on the CPU, overlapping the GPU stream. This seed
# has no existing reference; every other window reused one from the v5 sweep.
note "starting STEOM for frame 1832 (2 jobs, CPU) -- no reference exists for this seed"
$PY scripts/run_steom_v5.py --frames 1832 --concurrency 1 --nprocs 10 \
    --maxcore 3000 --hours 10 --min-free-gb 12 --min-disk-gb 150 \
    >> "$LOG/steom.log" 2>&1 &
STEOM_PID=$!

note "starting the v6min TDDFT campaign (125 frames, 3 block workers, ~3.5 h)"
./run_v6_min_campaign.sh >> "$LOG/tddft.log" 2>&1

if kill -0 "$STEOM_PID" 2>/dev/null; then
  note "TDDFT done; waiting on STEOM (pid $STEOM_PID)"
  wait "$STEOM_PID"; note "STEOM exit rc=$?"
fi
note "===== queue complete ====="
exit 0

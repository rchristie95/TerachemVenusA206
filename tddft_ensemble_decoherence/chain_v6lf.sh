#!/usr/bin/env bash
# Waits for the v6t2 (1 fs) campaign to finish, then builds and runs the 500 fs
# low-frequency extension. Launched under systemd so it survives the session.
#
# Refuses to start the extension if v6t2 did not actually complete: an aborted or
# retired run means something needs a human, and chaining more GPU hours onto a
# broken result is worse than stopping.
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
LOG=logs/v6lf; mkdir -p "$LOG"
S="$LOG/chain.log"
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$S"; }

count() {  # count <prefix>
  $PY -c "import csv,sys;print(sum(sum(r['convergence_status']=='complete' for r in csv.DictReader(open(f'manifests/{sys.argv[1]}{w}_frames.csv'))) for w in 'abc'))" "$1" 2>/dev/null || echo -1
}

note "===== chain armed: waiting on v6t2 ====="
while systemctl --user is-active --quiet v6t2.service; do sleep 120; done
N=$(count v6t2)
note "v6t2 unit finished at $N/300 complete"
if [ "$N" -lt 300 ]; then
  note "ABORT: v6t2 ended incomplete ($N/300). Not chaining 8 h of GPU onto a broken run."
  exit 1
fi

# ---- MD: 500 fs from each seed, 4 fs saving -------------------------------
note "building 500 fs MD windows"
set -- "900 -1209516.140068 20260901" "1400 -1209328.412976 20260902" "1900 -1209262.681868 20260903"
for spec in "$@"; do
  set -- $spec
  note "  MD seed frame $1"
  if ! bash ../tandem_v6_t2_campaign/run_md_ext.sh "$1" "$2" "$3" >> "$LOG/md_w$1.log" 2>&1; then
    note "ABORT: MD failed for seed $1 (see $LOG/md_w$1.log)"; exit 1
  fi
done

# ---- overlap check: did the MD reproduce the first 100 fs? ----------------
note "checking the 100 fs overlap against the v6t2 trajectories"
$PY scripts/check_lf_overlap.py >> "$S" 2>&1 || note "  (overlap check reported a mismatch -- see above; windows remain usable separately)"

# ---- manifests -------------------------------------------------------------
if ! $PY scripts/initialize_v6lf_manifests.py >> "$S" 2>&1; then
  note "ABORT: manifest initialisation failed"; exit 1
fi

# ---- run -------------------------------------------------------------------
note "starting the low-frequency TDDFT campaign (300 frames, ~8 h)"
exec ./run_v6_lf_campaign.sh

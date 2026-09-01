#!/usr/bin/env bash
# Second pass over the STEOM stream.
#
# run_steom_v5.py runs its eight jobs once. Anything the disk guard sheds, or
# that hits the per-job timeout, is left incomplete -- the driver skips a job
# only when its output says ORCA TERMINATED NORMALLY, so re-running it picks up
# exactly the gaps and nothing else.
#
# Deliberately a separate unit rather than an edit to run_v5_campaign.sh: bash
# reads a script incrementally, so editing one that is already running corrupts
# it mid-execution.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
LOG=logs/v5camp/steom_completion.log; mkdir -p logs/v5camp
note() { echo "[$(date '+%F %H:%M:%S')] $*" >> "$LOG"; }

done_count() { ls steom_v5/frame_*/site_*/steom_result.json 2>/dev/null | \
               xargs -r grep -l '"status": "complete"' 2>/dev/null | wc -l; }

note "waiting for the first STEOM pass to finish"
while pgrep -f "run_steom_v5.py" >/dev/null; do sleep 300; done
note "first pass over: $(done_count)/8 complete"

for pass in 1 2; do
  n=$(done_count)
  [ "$n" -ge 8 ] && { note "all 8 complete, no pass $pass needed"; break; }
  # ONE job at a time. A DLPNO-STEOM run of this class writes 250-300 GB of
  # PNO/PAO scratch (weekend_dipole/run_ddterm_chain.sh:59), which on a 437 GB
  # disk permits exactly one. Running four concurrently is what drove the burn
  # to 699 GB/h and cost frame 1400.
  note "completion pass $pass: $n/8 done, re-running the remainder one at a time"
  "$PY" scripts/run_steom_v5.py --concurrency 1 --nprocs 10 --maxcore 3000 \
        --hours 8 --min-free-gb 12 --min-disk-gb 150 >> logs/v5camp/steom.log 2>&1
  note "completion pass $pass finished: $(done_count)/8 complete"
done
note "STEOM stream final: $(done_count)/8 complete"

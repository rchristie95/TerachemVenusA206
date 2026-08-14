#!/usr/bin/env bash
# Autonomous weekend queue. Nothing here touches the network.
#
# GOAL. Two things the manuscript needs:
#   (a) Is |mu| = 9.6 D a basis artefact? Extinction and Strickler-Berg give
#       7.5-7.9 D, and J ~ |mu|^2, so this is a 48% lever on the coupling and
#       2.3x on the CD rotational strength -- the largest open systematic.
#   (b) The ladder's EOM-CCSD row is a 6-31G number sitting among def2-SVP
#       rows. Recomputing it at def2-SVP closes a referee-visible hole.
#
# DESIGN FOR UNATTENDED RUNNING
#   - Every job is independent. A failure is logged and the queue MOVES ON;
#     one bad job cannot eat the weekend or block the rest.
#   - Per-job wall-clock timeout, so a hung or thrashing job is killed rather
#     than starving everything behind it.
#   - Scratch is deleted after EVERY job. A previous run filled a 915 GB disk
#     with 467 GB of ORCA integral scratch and took out two unrelated jobs;
#     that must not recur. Density files are never touched.
#   - Disk and memory are checked before each job; the queue stops cleanly
#     rather than thrashing if either runs short.
#   - Resumable: a job whose .out already shows ORCA TERMINATED NORMALLY is
#     skipped, so restarting continues rather than repeating.
#   - Cheap jobs run first, so partial results exist early.
#
# Launch:  ./launch_weekend.sh          (user systemd unit, survives logout)
# Watch:   journalctl --user -u weekend -f
# Collect: python3 collect_weekend.py

set -u
cd "$(dirname "$0")" || exit 1

ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH

STATUS=weekend_status.txt
MIN_DISK_GB=60          # refuse to start a job below this
MIN_RAM_GB=20

note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }

free_disk_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }
free_ram_gb()  { free -g | awk 'NR==2{print $7}'; }

# Remove ORCA integral scratch but NEVER anything with a density in the name.
clean_scratch() {
  find . -maxdepth 1 -name "*.tmp*" ! -iname "*densit*" -delete 2>/dev/null
  find . -maxdepth 1 \( -name "*.bas[0-9]*" -o -name "*.hostnames" \) -delete 2>/dev/null
}

# run_job <name> <timeout_hours>
run_job() {
  local job=$1 hours=$2
  if [ -f "$job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$job.out"; then
    note "SKIP    $job (already complete)"; return 0
  fi
  local d r; d=$(free_disk_gb); r=$(free_ram_gb)
  if [ "$d" -lt "$MIN_DISK_GB" ]; then note "STOP    only ${d} GB disk free, need ${MIN_DISK_GB} -- halting queue"; return 2; fi
  if [ "$r" -lt "$MIN_RAM_GB" ]; then note "SKIP    $job (only ${r} GB RAM free)"; return 0; fi

  note "START   $job (timeout ${hours} h, ${d} GB disk, ${r} GB RAM)"
  timeout "${hours}h" "$ORCA/orca" "$job.inp" > "$job.out" 2>&1
  local rc=$?
  if grep -q "ORCA TERMINATED NORMALLY" "$job.out"; then
    note "OK      $job"
  elif [ $rc -eq 124 ]; then
    note "TIMEOUT $job after ${hours} h -- moving on"
  else
    note "FAIL    $job (rc=$rc) -- moving on; see $job.out"
  fi
  clean_scratch
  return 0
}

note "===== weekend queue start (pid $$) ====="
note "disk $(free_disk_gb) GB, RAM $(free_ram_gb) GB"

# Cheapest first so partial results exist early.
run_job tddft_def2_svp    2  || exit 0
run_job tddft_def2_svpd   3  || exit 0
run_job tddft_def2_tzvp   6  || exit 0
run_job tddft_def2_tzvpd 10  || exit 0

# The referee-facing gap.
run_job eomccsd_neutral_svp 20 || exit 0

# Correlated sweep. def2-SVPD reproduces the production setting.
run_job steom_def2_svp   12 || exit 0
run_job steom_def2_svpd  16 || exit 0
run_job steom_def2_tzvp  24 || exit 0

clean_scratch
note "===== weekend queue complete ====="
/home/robson/anaconda3/envs/adcc_env/bin/python collect_weekend.py >> "$STATUS" 2>&1 || true

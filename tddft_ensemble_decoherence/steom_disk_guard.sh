#!/usr/bin/env bash
# Resource throttle for the STEOM stream.
#
# THE BINDING FACT, already recorded at weekend_dipole/run_ddterm_chain.sh:59:
# a DLPNO-STEOM run of this class writes 250-300 GB of PNO/PAO scratch, and that
# chain guarded a SINGLE job with a 280-350 GB disk precondition. The footprint
# is set by the QM basis and the 134 correlated electrons, not by the 75,035
# point charges, so it cannot be reduced by trimming the embedding.
#
# On a 437 GB disk that means exactly ONE STEOM job at a time. This campaign was
# launched at four, which is why the burn hit 699 GB/h once frame 900 site A
# reached the EOM dressing.
#
# The cap is therefore on JOBS (one site of one frame), not on frames -- a frame
# is two jobs and would need 500-600 GB.
#
# SIGSTOP rather than kill: PNO/PAO integrals have no restart value, so killing a
# job at 60 GB throws away hours AND does not reduce pressure, because
# run_steom_v5.py's thread pool immediately starts the next queued job in the
# freed slot. Stopping holds the work and the slot. Thresholds on the symptom do
# not work either: an earlier version resumed on "swap is calm now", which was
# only true because it had just paused something, and it flapped straight back
# into thrash. Shed on the load, not on the symptom.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
LOG=logs/v5camp/disk_guard.log; mkdir -p logs/v5camp
MAX_JOBS=${MAX_JOBS:-2}       # concurrent STEOM jobs
# Measured on the completed frame 900 site A: peak scratch 81 GB at def2-SVP,
# NOT the 250-300 GB the repo note records -- that figure is for def2-SVPD and
# def2-TZVPD. Two jobs is ~170 GB steady, which leaves room for one job's
# transient during EOM intermediate generation (the spike that reached
# 699 GB/h with four running). Three would not.
PAUSE_GB=${PAUSE_GB:-140}     # disk emergency: shed regardless of the cap
RESUME_GB=${RESUME_GB:-230}   # 2 x 81 GB peak plus transient headroom
DWELL=${DWELL:-600}
note() { echo "[$(date '+%F %H:%M:%S')] $*" >> "$LOG"; }
free_gb() { df -BG --output=avail / | tail -1 | tr -dc 0-9; }

# Read state from the RANK PROCESS TABLE, not a bookkeeping file: a file drifts
# out of step across a restart, the process table cannot. 'T' is SIGSTOPped.
job_states() {   # -> "<frame> <site> running|paused", highest frame first
  for fd in $(ls -d steom_v5/frame_* 2>/dev/null | sort -r); do
    for s in A B; do
      tgt=$(readlink -f "$fd/site_$s" 2>/dev/null) || continue
      live=0; stopped=0
      for p in $(pgrep -f "orca_mdci|orca_leanscf" 2>/dev/null); do
        [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$tgt" ] || continue
        if [ "$(ps -o state= -p "$p" 2>/dev/null | tr -d ' ')" = "T" ]; then
          stopped=$((stopped+1)); else live=$((live+1)); fi
      done
      [ $((live + stopped)) -eq 0 ] && continue
      f=$(basename "$fd" | sed 's/frame_//')
      if [ "$live" -gt 0 ]; then echo "$f $s running"; else echo "$f $s paused"; fi
    done
  done
}
running_jobs() { job_states | awk '$3=="running"{print $1" "$2}'; }
paused_jobs()  { job_states | awk '$3=="paused"{print $1" "$2}'; }

declare -A PAUSED_AT
note "throttle start (max ${MAX_JOBS} JOB(s), shed<${PAUSE_GB} GB, resume>${RESUME_GB} GB, dwell ${DWELL}s, $(free_gb) GB free)"
while systemctl --user is-active --quiet v5camp.service; do
  avail=$(free_gb)
  nrun=$(running_jobs | wc -l)
  cap=$MAX_JOBS
  [ "$avail" -lt "$PAUSE_GB" ] && cap=0

  # Deadlock escape: scratch is reclaimed only when a job COMPLETES, so hitting
  # the floor with everything already paused can never recover on its own. This
  # is the one case where killing is right -- it trades one job for the rest.
  # run_steom_v5.py re-runs anything lacking ORCA TERMINATED NORMALLY, so the
  # casualty is picked up by steom_completion.sh.
  if [ "$avail" -lt "$PAUSE_GB" ] && [ "$nrun" -eq 0 ] && [ -n "$(paused_jobs)" ]; then
    set -- $(paused_jobs | head -1); df=$1; ds=$2
    note "DEADLOCK: ${avail} GB free, all jobs paused -> killing frame_${df}/site_${ds}"
    ./steom_pause_site.sh cont "$df" "$ds" >> "$LOG" 2>&1
    sleep 5
    tgt=$(readlink -f "steom_v5/frame_${df}/site_${ds}")
    for p in $(pgrep -f "orca|mpi"); do
      [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$tgt" ] && kill -TERM "$p" 2>/dev/null
    done
    sleep 25
    note "reclaimed to $(free_gb) GB free"
    sleep 120; continue
  fi

  if [ "$nrun" -gt "$cap" ]; then
    set -- $(running_jobs | head -1); vf=$1; vs=$2
    note "over cap: ${nrun} running, cap ${cap} (${avail} GB free) -> pausing frame_${vf}/site_${vs}"
    ./steom_pause_site.sh stop "$vf" "$vs" >> "$LOG" 2>&1
    PAUSED_AT[${vf}_${vs}]=$(date +%s)
  elif [ "$nrun" -lt "$MAX_JOBS" ] && [ "$avail" -gt "$RESUME_GB" ]; then
    # Lowest frame first: the earliest-paused job is the furthest along, and its
    # `timeout Nh` wrapper keeps counting the whole time it is stopped.
    paused_jobs | sort -n | while read -r rf rs; do
      since=${PAUSED_AT[${rf}_${rs}]:-0}
      if [ $(( $(date +%s) - since )) -ge "$DWELL" ]; then
        note "under cap: ${nrun} running (${avail} GB free) -> resuming frame_${rf}/site_${rs}"
        ./steom_pause_site.sh cont "$rf" "$rs" >> "$LOG" 2>&1
        break
      fi
    done
  fi
  sleep 120
done
# Never exit leaving a SIGSTOPped ORCA behind: it would hang until reboot.
paused_jobs | while read -r rf rs; do
  note "shutdown: resuming frame_${rf}/site_${rs}"; ./steom_pause_site.sh cont "$rf" "$rs" >> "$LOG" 2>&1
done
note "throttle stop ($(free_gb) GB free)"

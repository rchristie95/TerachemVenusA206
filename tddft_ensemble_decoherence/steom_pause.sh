#!/usr/bin/env bash
# Pause or resume one STEOM job pair without losing its scratch.
#
# SIGSTOP rather than kill: the PNO/PAO integral files have no restart value, so
# killing a job at 30 GB throws away every hour of it. Stopping halts the writes
# and holds the work. Every process of the pair is signalled, ranks and their
# mpirun/orca/timeout parents together, so no rank is left waiting on a stopped
# peer. Identification is by /proc/PID/cwd, not a pgrep pattern -- a pattern
# broad enough to catch the ranks also matches this script's own command line.
#
# Caveat: the `timeout Nh` wrapper keeps counting while stopped, so a long pause
# eats the job's remaining budget.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
ACTION=${1:?usage: steom_pause.sh <stop|cont> <frame>}
FRAME=${2:?usage: steom_pause.sh <stop|cont> <frame>}
SIG=$([ "$ACTION" = stop ] && echo STOP || echo CONT)
n=0
for s in A B; do
  tgt=$(readlink -f "steom_v5/frame_${FRAME}/site_$s" 2>/dev/null) || continue
  [ -z "$tgt" ] && continue
  for p in $(pgrep -f "orca|mpi"); do
    [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$tgt" ] || continue
    kill -"$SIG" "$p" 2>/dev/null && n=$((n+1))
  done
done
echo "[$(date '+%F %H:%M:%S')] SIG$SIG -> $n processes of frame_${FRAME}" \
  | tee -a logs/v5camp/steom_pause.log

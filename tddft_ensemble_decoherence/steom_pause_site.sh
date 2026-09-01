#!/usr/bin/env bash
# SIGSTOP/SIGCONT one SITE of one STEOM frame. See steom_pause.sh for why
# stopping beats killing and why /proc/PID/cwd is used instead of a pgrep pattern.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
ACTION=${1:?usage: steom_pause_site.sh <stop|cont> <frame> <A|B>}
FRAME=${2:?}; SITE=${3:?}
SIG=$([ "$ACTION" = stop ] && echo STOP || echo CONT)
tgt=$(readlink -f "steom_v5/frame_${FRAME}/site_${SITE}") || exit 1
n=0
for p in $(pgrep -f "orca|mpi"); do
  [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$tgt" ] || continue
  kill -"$SIG" "$p" 2>/dev/null && n=$((n+1))
done
echo "[$(date '+%F %H:%M:%S')] SIG$SIG -> $n procs of frame_${FRAME}/site_${SITE}" \
  | tee -a logs/v5camp/steom_pause.log

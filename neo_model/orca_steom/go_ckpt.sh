#!/bin/bash
# Checkpoint-aware ORCA STEOM runner.  Usage: go_ckpt.sh <input.inp>
# PRESERVES <base>.gbw + ALL scratch (never deletes) so a re-run can reuse:
#   * SCF: ORCA AutoStart reads a matching <base>.gbw for a single-point calculation.
#   * MDCI/EOM: ORCA 6.1 does not accept a job-level `Restart true` keyword in %mdci.
#     If an EOM root exhausts MaxIter, keep the failed output for provenance and rerun
#     through the normal ORCA driver with the manual-recommended larger MaxIter. The
#     Davidson solver performs its own internal subspace restarts at the NDav limit.
# Snapshots <base>.ckpt.gbw every few min so a mid-run crash/watchdog-kill still leaves a usable
# SCF checkpoint. Preserve the failed output before rerunning through the normal ORCA driver.
cd /home/robson/PetaChem/neo_model/orca_steom
ORCABIN=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
OMPI=/home/robson/anaconda3/envs/openmpi416
export PATH="$OMPI/bin:$ORCABIN:$PATH"
export LD_LIBRARY_PATH="$OMPI/lib:$ORCABIN:$LD_LIBRARY_PATH"

INP="${1:?need input.inp}"; BASE="${INP%.inp}"; OUT="${BASE}.out"
echo "[ckpt start $(date '+%F %T')] $INP | AutoStart .gbw present: $([ -f ${BASE}.gbw ] && echo YES || echo no)" | tee -a "$OUT"

"$ORCABIN/orca" "$INP" >> "$OUT" 2>&1 &
PID=$!; echo "$PID" > ".orca_${BASE}_pid"

# periodic SCF-orbital checkpoint snapshot
( while kill -0 "$PID" 2>/dev/null; do [ -f "${BASE}.gbw" ] && cp -f "${BASE}.gbw" "${BASE}.ckpt.gbw" 2>/dev/null; sleep 180; done ) &
SNAP=$!
# RAM watchdog (preserve the box; scratch is kept so the run is resumable)
( while kill -0 "$PID" 2>/dev/null; do
    a=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
    [ "$a" -lt 3000 ] && { echo "[WATCHDOG $(date +%T)] MemAvailable ${a}MB<3000 -> kill ORCA (scratch kept for resume)" | tee -a "$OUT"; pkill -9 orca; break; }
    sleep 30
  done ) &
WD=$!

wait "$PID"; RC=$?
kill "$SNAP" "$WD" 2>/dev/null
[ -f "${BASE}.gbw" ] && cp -f "${BASE}.gbw" "${BASE}.ckpt.gbw"
echo "[ckpt done $(date '+%F %T')] rc=$RC | SCF checkpoint=${BASE}.ckpt.gbw + scratch preserved" | tee -a "$OUT"
echo "============================================================" | tee -a "$OUT"
grep -aiE "TERMINATED NORMALLY|error termination|equations failed|restart" "$OUT" | tail -4 | tee -a "$OUT"
echo "--- FINAL STEOM-CCSD ABSORPTION SPECTRUM (bright = max fosc) ---" | tee -a "$OUT"
awk '/ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS/{n=NR}{a[NR]=$0}END{for(i=n;i<=n+11&&i<=NR;i++)print a[i]}' "$OUT" | tee -a "$OUT"

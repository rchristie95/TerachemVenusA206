#!/bin/bash
# Checkpoint-aware ORCA STEOM runner.  Usage: go_ckpt.sh <input.inp>
# PRESERVES <base>.gbw + ALL scratch (never deletes) so a re-run can resume:
#   * SCF  : ORCA AutoStart auto-reads <base>.gbw (skips the SCF) when the basename matches.
#   * CCSD : on a RESUME, add  '%mdci Restart true end'  to the input to reuse stored amplitudes
#            from the preserved scratch (best-effort for DLPNO; verify it prints "restart").
# Snapshots <base>.ckpt.gbw every few min so a mid-run crash/watchdog-kill still leaves a usable
# checkpoint. To resume after any failure: just re-run the SAME input with go_ckpt.sh.
cd /home/robson/PetaChem/neo_model/orca_steom
ORCABIN=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
OMPI=/home/robson/anaconda3/envs/openmpi416
export PATH="$OMPI/bin:$ORCABIN:$PATH"
export LD_LIBRARY_PATH="$OMPI/lib:$ORCABIN:$LD_LIBRARY_PATH"

INP="${1:?need input.inp}"; BASE="${INP%.inp}"; OUT="${BASE}.out"
echo "[ckpt start $(date '+%F %T')] $INP | AutoStart .gbw present: $([ -f ${BASE}.gbw ] && echo YES || echo no) | Restart-flag in input: $(grep -qiE 'Restart[[:space:]]+true' $INP && echo YES || echo no)" | tee -a "$OUT"

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
echo "[ckpt done $(date '+%F %T')] rc=$RC | checkpoint=${BASE}.ckpt.gbw + scratch preserved (resume: re-run same input)" | tee -a "$OUT"
echo "============================================================" | tee -a "$OUT"
grep -aiE "TERMINATED NORMALLY|error termination|equations failed|restart" "$OUT" | tail -4 | tee -a "$OUT"
echo "--- FINAL STEOM-CCSD ABSORPTION SPECTRUM (bright = max fosc) ---" | tee -a "$OUT"
awk '/ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS/{n=NR}{a[NR]=$0}END{for(i=n;i<=n+11&&i<=NR;i++)print a[i]}' "$OUT" | tee -a "$OUT"

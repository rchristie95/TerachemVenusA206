#!/bin/bash
# Parallel DLPNO-STEOM-CCSD runner.  Usage: go_par.sh <input.inp>
# 16 ranks FAILS IP-EOM on the 44-atom model (over-parallelized); use <=8.
cd /home/robson/PetaChem/neo_model/orca_steom
ORCABIN=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
OMPI=/home/robson/anaconda3/envs/openmpi416
export PATH="$OMPI/bin:$ORCABIN:$PATH"
export LD_LIBRARY_PATH="$OMPI/lib:$ORCABIN:$LD_LIBRARY_PATH"

INP="${1:-steom_svp_par8.inp}"
BASE="${INP%.inp}"
OUT="${BASE}.out"
NP=$(awk '/nprocs/{for(i=1;i<=NF;i++) if($i=="nprocs") print $(i+1)}' "$INP")
: > "$OUT"
echo "[parallel STEOM-CCSD start $(date '+%F %T')] input=$INP nprocs=$NP" | tee -a "$OUT"

"$ORCABIN/orca" "$INP" >> "$OUT" 2>&1 &
ORCA_PID=$!
echo "$ORCA_PID" > ".orca_${BASE}_pid"
echo "[ORCA driver PID=$ORCA_PID]" | tee -a "$OUT"

# RAM watchdog: kill ORCA if MemAvailable drops below 3 GB (prevent box crash).
( while kill -0 "$ORCA_PID" 2>/dev/null; do
    avail=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
    if [ "$avail" -lt 3000 ]; then
      echo "[WATCHDOG $(date '+%T')] MemAvailable ${avail}MB < 3000MB -> KILLING ORCA to save the box" | tee -a "$OUT"
      pkill -9 orca 2>/dev/null; break
    fi
    sleep 30
  done ) &
WD_PID=$!

wait "$ORCA_PID"; RC=$?
kill "$WD_PID" 2>/dev/null
echo "[parallel STEOM-CCSD done $(date '+%F %T')] orca_exit=$RC (NB: check 'error termination' below; driver RC is unreliable)" | tee -a "$OUT"
echo "============================================================" | tee -a "$OUT"
grep -iE "TERMINATED NORMALLY|error termination|aborting|equations failed" "$OUT" | tail -4 | tee -a "$OUT"
echo "--- FINAL STEOM-CCSD ABSORPTION SPECTRUM (bright = max fosc) ---" | tee -a "$OUT"
awk '/ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS/{n=NR} {a[NR]=$0} END{for(i=n;i<=n+11 && i<=NR;i++) print a[i]}' "$OUT" | tee -a "$OUT"

#!/bin/bash
# Set ORCABIN to your ORCA install dir (e.g. /home/robson/orca_6_1_0), then:  ./run_orca.sh steom_svp.inp
set -e
ORCABIN="${ORCABIN:-/home/robson/orca}"
export PATH="$ORCABIN:$PATH"
export LD_LIBRARY_PATH="$ORCABIN:$LD_LIBRARY_PATH"
# OpenMPI that matches the ORCA build must also be on PATH/LD_LIBRARY_PATH (see README)
INP="${1:-steom_svp.inp}"; OUT="${INP%.inp}.out"
echo "[orca start $(date)] $INP -> $OUT"
"$ORCABIN/orca" "$INP" > "$OUT" 2>&1
echo "[orca done $(date)]  bright state:"
grep -iE "STATE|nm|fosc|oscillator" "$OUT" | tail -20

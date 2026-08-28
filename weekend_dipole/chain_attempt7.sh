#!/usr/bin/env bash
# Item 3 must wait for item 1: the serial TZVPD run needs 43-53 GB and the
# gas-phase SVPD job is holding up to 48 GB. Poll for item 1, then start.
set -u
cd "$(dirname "$0")" || exit 1
GAS=/home/robson/PetaChem/gasphase_steom/gas_steom_status.txt
echo "[$(date '+%F %H:%M:%S')] waiting for item 1 to release memory" > chain_status.txt
while ! grep -q "complete" "$GAS" 2>/dev/null; do sleep 300; done
echo "[$(date '+%F %H:%M:%S')] item 1 done; starting TZVPD attempt 7" >> chain_status.txt
cp steom_def2_tzvpd_attempt7.inp steom_def2_tzvpd.inp
exec ./run_steom_tzvpd.sh

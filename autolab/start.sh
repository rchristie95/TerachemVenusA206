#!/usr/bin/env bash
# Launch the compute loop detached, so it survives a closed terminal, a dead
# session and a rate limit. Everything numerical runs in the TeraChem env
# (OpenMM 8.4 + CUDA); ORCA needs its own LD_LIBRARY_PATH.
#
#   ./autolab/start.sh              # drain the queue and exit
#   ./autolab/start.sh --idle-forever   # stay up, waiting for new work
#   touch autolab/STOP              # halt after the current job
#
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${AUTOLAB_PYTHON:-$HOME/anaconda3/envs/TeraChem/bin/python}"
ORCA="$REPO/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"

cd "$REPO"

if [[ -f autolab/STOP ]]; then
  echo "autolab/STOP exists; refusing to start. Remove it first." >&2
  exit 1
fi

if pgrep -f "autolab/driver.py" >/dev/null; then
  echo "driver already running (pid $(pgrep -f 'autolab/driver.py' | tr '\n' ' ')); not starting a second one." >&2
  exit 1
fi

export LD_LIBRARY_PATH="$ORCA:${LD_LIBRARY_PATH:-}"
export PATH="$ORCA:$PATH"

mkdir -p autolab/runs autolab/results
nohup "$PY" autolab/driver.py "$@" >> autolab/driver.out 2>&1 &
echo "driver started, pid $!"
echo "  log:    tail -f $REPO/autolab/driver.log"
echo "  status: $PY autolab/ledger.py status"
echo "  halt:   touch $REPO/autolab/STOP"

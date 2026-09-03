#!/usr/bin/env bash
# The triple- and quadruple-zeta rungs of the transition-dipole ladder.
#
# Runs AFTER the dense sweep and with a SINGLE worker. Both are needed:
#   * def2-TZVPD is 1292 basis functions with F shells. Sharing the card with
#     three SVP workers left 1736 MW of GPU memory and the job died silently
#     right after basis setup. This rung needs the card to itself.
#   * TeraChem's default 'generate' guess dies at this size with "Too many
#     electrons for guess generation code" (guess.cpp:58); the templates now
#     carry 'guess sad'.
#
# def2-SVPD is already complete (48/48 site-frames, <|mu|> = 9.370 +/- 0.626 D
# against production def2-SVP's 9.40), so diffuse functions are settled and this
# chain answers only the remaining question: does |mu| come down with zeta?
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
LOGS=logs/ladder; mkdir -p "$LOGS"
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOGS/heavy_chain.log"; }

note "waiting for the dense sweep"
while systemctl --user is-active v8dense-chain >/dev/null 2>&1; do sleep 180; done
note "dense sweep finished; GPU free"
sleep 60

set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export OMP_NUM_THREADS=16          # single worker, so give it the cores

note "starting def2-TZVPD, 8 frames x 2 sites, one worker"
"$PY" -u scripts/run_basis_ladder.py 0 a 1 tzvpd >> "$LOGS/heavy_tzvpd.log" 2>&1
note "TZVPD exit rc=$?"

note "starting def2-QZVPD, 1 frame x 2 sites (convergence indicator)"
"$PY" -u scripts/run_basis_ladder.py 0 a 1 qzvpd >> "$LOGS/heavy_qzvpd.log" 2>&1
note "QZVPD exit rc=$?"
note "===== heavy ladder rungs complete ====="
exit 0

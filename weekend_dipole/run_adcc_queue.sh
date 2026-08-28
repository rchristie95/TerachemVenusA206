#!/usr/bin/env bash
# ADC(2) dipole triangulation -- the independent third opinion on the transition
# dipole that decides whether J stays 32.82 or drops toward 21 cm^-1.
# Runs on the weekend_dipole geometry + field, so it is directly comparable to
# the STEOM and TDDFT numbers already in this directory.
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/adcc_env/bin/python
GEOM=geom_cthrp.xyz                       # relaxed C-O 1.263 -- NOT the MM-min geom
MM=/home/robson/PetaChem/neo_model/mm_light_cthrp.dat   # identical to field.pc
LOG=adcc_queue_status.txt
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }

# libxm pagefiles can grow large; sweep stale ones and refuse to start without
# headroom, since the GPU campaign shares this disk.
FREE_GB=$(df --output=avail -BG . | tail -1 | tr -dc '0-9')
if [ "$FREE_GB" -lt 120 ]; then echo "REFUSING: only ${FREE_GB} GB disk free" >&2; exit 1; fi
note "===== adcc ADC(2) dipole queue start (pid $$) ====="
note "disk ${FREE_GB} GB free (libxm pagefiles land in adcc_scratch/)"
note "RAM $(free -g | awk 'NR==2{print $7}') GB available"

# adcc needs FAR more than the naive o^2v^2 estimate: def2-svp with all 387
# virtuals peaked at 33.9 GB and was cgroup-OOM-killed. Freeze the HIGHEST
# virtuals down to a constant 250 active in every basis -- this (a) fits in
# ~14 GB, (b) keeps the comparison clean because the active space is identical
# across bases, and (c) preserves the diffuse manifold, which is LOW-lying
# (SVPD LUMO sits at 0.94 eV) and so is never among the frozen high virtuals.
# Caveat: truncating high virtuals biases the excitation ENERGY blue; the
# transition DIPOLE, which is what this queue exists to measure, is dominated
# by occupied -> low-virtual character and is far less sensitive.
#     basis        nfc nfv  nstates  hours    (nvir: svp 387, svpd 609, tzvpd 1069)
JOBS=("def2-svp    26  187  3        12"
      "def2-svpd   26  409  3        18"
      "def2-tzvpd  26  869  3        24")

for J in "${JOBS[@]}"; do
  set -- $J; B=$1; NFC=$2; NFV=$3; NST=$4; HRS=$5
  OUT="adcc_${B//-/_}.out"
  if grep -qE '^\[[0-9]+s\] ALL DONE' "$OUT" 2>/dev/null; then note "SKIP    $B (already complete)"; continue; fi
  note "START   $B (frozen_core=$NFC frozen_virtual=$NFV nstates=$NST, ${HRS}h cap)"
  timeout "${HRS}h" "$PY" run_adcc_dipole.py "$GEOM" "$MM" "$B" "$NFC" "$NFV" "$NST" > "$OUT" 2>&1
  rc=$?
  if grep -qE '^\[[0-9]+s\] ALL DONE' "$OUT"; then
    note "OK      $B -> $(grep -E '^BRIGHT STATE' "$OUT" | head -1)"
  elif [ $rc -eq 124 ]; then note "TIMEOUT $B after ${HRS} h"
  else note "FAIL    $B (rc=$rc) -- see $OUT"; fi
done
note "===== adcc queue complete ====="

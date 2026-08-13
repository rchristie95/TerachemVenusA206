#!/usr/bin/env bash
# Neutral gas-phase benchmark ladder -- detached and resumable.
#
# Every rung uses the SAME geometry (neutral_chromophore.xyz), the SAME basis
# (def2-SVP) and the same neutral closed-shell chromophore, so the ladder is
# method-matched throughout. The neutral form is used because the isolated
# ANION has no bound electron-attached state in a diffuse basis, which breaks
# the EA-EOM sector that STEOM needs.
#
# Rungs:
#   1  TDDFT CAM-B3LYP        (ORCA)
#   2  TDDFT wB97X-D3         (ORCA)
#   3  DLPNO-STEOM-CCSD       (ORCA)   <- method of record
#   4  canonical EOM-CCSD     (ORCA)
#   5  ADC(2)                 (adcc)   <- independent hierarchy, replaces the
#                                         Q-Chem EOM-CCSD(fT) rung, which cannot
#                                         be run (licence expired 2026-07-25 and
#                                         ORCA has no excited-state triples)
#
# RESILIENCE. Launch with launch.sh, which runs this as a USER SYSTEMD UNIT.
# setsid+nohup is NOT enough -- the agent harness kills its whole process tree
# on exit and a setsid child is still a descendant. systemd reparents the job
# outside that tree. Nothing here touches the network. Each rung is skipped if
# its output already shows normal termination, so an interruption resumes
# rather than restarting.

set -u
cd "$(dirname "$0")" || exit 1

ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
PY=/home/robson/anaconda3/envs/adcc_env/bin/python
STATUS=ladder_status.txt

stamp() { date "+%Y-%m-%d %H:%M:%S"; }
note()  { echo "[$(stamp)] $*" | tee -a "$STATUS"; }

orca_done() { [ -f "$1.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$1.out"; }

run_orca() {
  local job=$1
  if orca_done "$job"; then note "SKIP  $job (already terminated normally)"; return 0; fi
  note "START $job"
  "$ORCA/orca" "$job.inp" > "$job.out" 2>&1
  if orca_done "$job"; then note "OK    $job"; else note "FAIL  $job (see $job.out)"; fi
}

note "===== ladder start (pid $$) ====="
run_orca tddft_camb3lyp
run_orca tddft_wb97xd3
run_orca steom_neutral_check
run_orca eomccsd_neutral

if [ -f adc2_neutral.json ]; then
  note "SKIP  adc2 (adc2_neutral.json exists)"
else
  note "START adc2"
  "$PY" adc2_neutral.py > adc2_neutral.out 2>&1
  if [ -f adc2_neutral.json ]; then note "OK    adc2"; else note "FAIL  adc2 (see adc2_neutral.out)"; fi
fi

note "===== ladder complete ====="
"$PY" collect_ladder.py >> "$STATUS" 2>&1 || true

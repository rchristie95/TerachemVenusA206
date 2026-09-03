#!/usr/bin/env bash
# Dense-segment chain: wait for the basis ladder, then MD, then the QM sweep.
#
# 1 ps at 2 fs seeded from v5camp frame 1756 -- the ensemble MEDIAN detuning, so
# a typical configuration rather than the near-degenerate tail the T2 windows
# used. The v6t2 windows are 100 fs and resolve only to 334 cm^-1, coarser than
# kT = 208; 1 ps resolves to 33 cm^-1, which is what the second-order cumulant
# needs to separate the fast lineshape from the static disorder without the two
# being fitted against one another.
#
# Runs AFTER the ladder rather than beside it: both are GPU bound, so sharing
# the card would only split it.
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
LOGS=logs/v8dense; mkdir -p "$LOGS"
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOGS/chain.log"; }

note "waiting for the basis ladder to finish"
while systemctl --user is-active basis-ladder >/dev/null 2>&1; do sleep 120; done
note "ladder finished; starting MD"

# 0.25 fs timestep as in the T2 windows; 4000 steps = 1 ps; save every 8 = 2 fs.
"$PY" -u ../run_nvt.py \
  --pdb ../tc_tandem_v5_campaign/solvated_protonated.pdb \
  --amber-cr2-prmtop ../anionic_build/monomer_solv.prmtop \
  --workdir ../tc_tandem_v8_dense \
  --skip-solvation \
  --initial-dcd ../tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd \
  --initial-dcd-frame 1756 \
  --expect-initial-potential-kj -1206937.727741 \
  --npt-steps 0 --timestep-fs 0.25 --nvt-steps 4000 \
  --openmm-trajectory-interval 8 \
  --openmm-trajectory-file tandem_v8_1ps_2fs.dcd \
  --checkpoint-file nvt.chk --checkpoint-interval 800 \
  --seed 20260903 --no-video --skip-final-minimization \
  >> "$LOGS/md.log" 2>&1
rc=$?; note "MD exit rc=$rc"
[ $rc -ne 0 ] && { note "MD failed, stopping"; exit 1; }

"$PY" scripts/initialize_v8dense_manifest.py >> "$LOGS/chain.log" 2>&1 || {
  note "manifest init failed"; exit 1; }
"$PY" -c "
import csv
rows=list(csv.DictReader(open('manifests/v8dense_frames.csv'))); cols=list(rows[0])
for i,w in enumerate('abc'):
    with open(f'manifests/v8dense_{w}_frames.csv','w',newline='') as f:
        c=csv.DictWriter(f,fieldnames=cols); c.writeheader(); c.writerows(rows[i::3])
print('  split into a/b/c')" >> "$LOGS/chain.log" 2>&1

note "starting the QM sweep, 3 workers"
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export OMP_NUM_THREADS=4
for w in a b c; do
  "$PY" -u scripts/run_ensemble.py --window v8dense --manifest-suffix "_$w" \
        --profile svp --roots 7 --gpu 0 >> "$LOGS/worker_$w.log" 2>&1 &
  sleep 20
done
wait
note "===== v8dense complete ====="
exit 0

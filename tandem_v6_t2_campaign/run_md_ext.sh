#!/usr/bin/env bash
# 500 fs low-frequency window at 4 fs cadence, from the same seed geometry as
# the matching v6t2 window.
#
#   2000 steps x 0.25 fs = 500 fs, DCD every 16 steps = 4 fs = 125 frames.
#
# Same seed frame, same velocity seed and same integrator seed as the 100 fs run,
# so the first 100 fs should reproduce it exactly. That is CHECKED, not assumed:
# frame 24 here (t = 100 fs) must hash to frame 99 of the v6t2 DCD. If it does,
# the two datasets splice into one 500 fs record; if not, this stands as its own
# window and only the >100 fs frames are used either way.
set -eu
cd /home/robson/PetaChem || exit 1
FRAME=$1; EXPECT=$2; VSEED=$3
PY=/home/robson/anaconda3/envs/TeraChem/bin/python

$PY -u run_nvt.py \
  --pdb tc_tandem_v5_campaign/solvated_protonated.pdb \
  --skip-solvation \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir "tc_tandem_v6_lf_w${FRAME}" \
  --initial-dcd tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd \
  --initial-dcd-frame "$FRAME" \
  --expect-initial-potential-kj "$EXPECT" \
  --npt-steps 0 \
  --nvt-steps 2000 \
  --timestep-fs 0.25 \
  --temperature-k 300.0 \
  --openmm-trajectory-interval 16 \
  --openmm-trajectory-file "tandem_v6_w${FRAME}_500fs_4fs.dcd" \
  --protein-topology-file protein_topology.pdb \
  --checkpoint-interval 0 \
  --skip-final-minimization \
  --seed "$VSEED" \
  --no-video

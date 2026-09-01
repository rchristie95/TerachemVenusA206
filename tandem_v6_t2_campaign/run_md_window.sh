#!/usr/bin/env bash
# One 100 fs / 1 fs high-cadence window, seeded on a v5 frame that already has
# both STEOM sites computed.
#
#   400 steps x 0.25 fs = 100 fs, DCD every 4 steps = 1 fs = 100 frames.
#
# The topology PDB is the v5 run's OWN solvated_protonated.pdb with
# --skip-solvation. It cannot be re-solvated: addSolvent is unseeded (--seed
# reaches only the Langevin integrator), so a rebuild would place different
# waters and the atom ordering would no longer match the DCD being seeded from.
#
# Usage: run_md_window.sh <frame> <expected_PE_kJ> <velocity_seed>
set -eu
cd /home/robson/PetaChem || exit 1
FRAME=$1; EXPECT=$2; VSEED=$3
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
WD=tc_tandem_v6_t2_w${FRAME}

$PY -u run_nvt.py \
  --pdb tc_tandem_v5_campaign/solvated_protonated.pdb \
  --skip-solvation \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir "$WD" \
  --initial-dcd tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd \
  --initial-dcd-frame "$FRAME" \
  --expect-initial-potential-kj "$EXPECT" \
  --npt-steps 0 \
  --nvt-steps 400 \
  --timestep-fs 0.25 \
  --temperature-k 300.0 \
  --openmm-trajectory-interval 4 \
  --openmm-trajectory-file "tandem_v6_w${FRAME}_100fs_1fs.dcd" \
  --protein-topology-file protein_topology.pdb \
  --checkpoint-interval 0 \
  --skip-final-minimization \
  --seed "$VSEED" \
  --no-video

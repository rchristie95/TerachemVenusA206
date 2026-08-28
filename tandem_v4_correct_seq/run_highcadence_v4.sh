#!/usr/bin/env bash
# High-cadence continuation on the CORRECTED 507-aa trajectory, for the site-energy
# ACF / homogeneous linewidth. Mirrors the reference tandem_nvt_8ps_4fs.dcd exactly:
#   2 fs timestep, --openmm-trajectory-interval 2 => a frame every 4 fs
#   --nvt-steps is a TOTAL, not an increment: the checkpoint is at step 500000,
#   so 504000 runs 4000 NEW steps => 8 ps => 2000 frames (this is why the
#   original reference command also used 504000)  (reference header: nset=2000, nsavc=2)
# Resumes the 1 ns checkpoint rather than restarting from the PDB, which both
# preserves the thermal state (run_nvt.py skips minimisation and velocity
# reseeding on resume) and avoids the documented stale-box trap -- loadCheckpoint
# restores the production box, whereas solvated_protonated.pdb carries the pre-NPT one.
set -u
cd /home/robson/PetaChem || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
S=tandem_v4_correct_seq/highcadence_status.txt
echo "[$(date '+%F %H:%M:%S')] high-cadence start" > $S
$PY -u run_nvt.py \
  --pdb tc_tandem_v4_correctseq/solvated_protonated.pdb \
  --skip-solvation \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir tc_tandem_v4_highcadence \
  --resume-checkpoint tc_tandem_v4_correctseq/nvt.chk \
  --nvt-steps 504000 \
  --openmm-trajectory-interval 2 \
  --openmm-trajectory-file tandem_v4_8ps_4fs.dcd \
  --protein-topology-file protein_topology.pdb \
  --checkpoint-file nvt.chk --checkpoint-interval 1000 \
  --skip-final-minimization \
  --seed 20260618 --no-video >> tandem_v4_correct_seq/highcadence.log 2>&1
echo "[$(date '+%F %H:%M:%S')] high-cadence exit rc=$?" >> $S

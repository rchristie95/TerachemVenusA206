#!/usr/bin/env bash
# Does freezing the X-H stretches change the coordinate that carries the gap?
#
# constraints=HBonds removes X-H stretch motion at ~3000 cm^-1, inside the band
# a ~12 fs T2 is sensitive to. Two 500 fs runs from the same seed frame, saving
# every 1 fs, identical but for the constraint set. 500 rather than 100 fs
# because the bridge C-C stretch has a ~20-30 fs period and 100 fs would give
# only ~4 of them.
#
# The comparison is on the METHINE BRIDGE LENGTH and BLA, not on an
# electrostatic proxy. The v5 bridge analysis put 25.9% of the detuning variance
# on the bridge and only 4.1% on electrostatics, and a Stark proxy holds the
# chromophore geometry rigid -- so it is blind by construction to exactly the
# channel this control is about.
set -eu
cd /home/robson/PetaChem
PY=/home/robson/anaconda3/envs/TeraChem/bin/python

run () {  # run <tag> <extra flags...>
  local tag=$1; shift
  $PY -u run_nvt.py \
    --pdb tc_tandem_v5_campaign/solvated_protonated.pdb --skip-solvation \
    --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
    --workdir "tc_tandem_v6_ctrl_${tag}" \
    --initial-dcd tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd --initial-dcd-frame 900 \
    --npt-steps 0 --nvt-steps 2000 --timestep-fs 0.25 --temperature-k 300.0 \
    --openmm-trajectory-interval 4 --openmm-trajectory-file "ctrl_${tag}.dcd" \
    --protein-topology-file protein_topology.pdb --checkpoint-interval 0 \
    --skip-final-minimization --seed 20260901 --no-video "$@"
}

run constrained
run unconstrained --no-hbond-constraints
echo "CONTROL_MD_DONE"

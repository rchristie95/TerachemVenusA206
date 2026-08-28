#!/usr/bin/env bash
# Rebuild the tandem MD with the CORRECT 511-aa Addgene #250938 sequence.
#
# Structure: tandem_dimer_v4_af3tails.pdb
#   = tandem_dimer_v2.pdb (correct linker Gly)
#   + LGMDELYK restored at both barrel C-termini (PDBFixer)
#   + those two tails' conformations taken from AlphaFold3 (fold_dvenus_td511_s1,
#     model_0), superposed on the 15 preceding barrel residues only (RMSD 0.42-0.44 A).
#     AF3's INTER-DOMAIN REGISTER IS DELIBERATELY NOT IMPORTED -- it is 11.7 A from
#     crystal, and the coupling depends on the A206 crystal register. Verified:
#     CR2-CR2 separation is 25.388 A both before and after grafting.
#
# Force field left at the defaults (amber14-all / tip3pfb, 2 fs) to match the
# trajectory this replaces. NPT added because this is a fresh box, unlike the
# original which resumed from an already-equilibrated checkpoint.
#
# The -2 e neutralisation bug cannot recur: run_nvt.py now seeds the CR2 template
# with its true -1 e BEFORE addSolvent counts ions, and aborts if |net| > 1e-4 e.
set -u
cd /home/robson/PetaChem || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
S=tandem_v4_correct_seq/md_v4_status.txt
echo "[$(date '+%F %H:%M:%S')] v4 MD start" > $S
$PY -u run_nvt.py \
  --pdb tandem_v4_correct_seq/tandem_dimer_v4_af3tails.pdb \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir tc_tandem_v4_correctseq \
  --add-solvent --padding-a 10.0 --ionic-strength-m 0.15 \
  --npt-steps 100000 \
  --nvt-steps 500000 \
  --openmm-trajectory-interval 500 \
  --openmm-trajectory-file tandem_v4_1ns.dcd \
  --protein-topology-file protein_topology.pdb \
  --checkpoint-file nvt.chk --checkpoint-interval 1000 \
  --seed 20260618 --no-video >> tandem_v4_correct_seq/md_v4.log 2>&1
echo "[$(date '+%F %H:%M:%S')] v4 MD exit rc=$?" >> $S

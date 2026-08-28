#!/usr/bin/env bash
# Ensemble coupling, corrected-sequence trajectory vs the old truncated one.
# Uses steom_transdens_specnorm_oldframe.npz for BOTH, not the manuscript's
# cap-masked density: only 4 of the 7 NTO pairs still have cubes on disk
# (98.45% vs 99.046% occupation), so the cap-masked density cannot be rebuilt
# exactly. Absolute values therefore will not equal the published 32.8, but the
# NEW/OLD ratio is clean because the density is identical on both sides.
# GPU only -- the ORCA queue owns the RAM.
set -u
cd /home/robson/PetaChem || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
D=neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz
S=tandem_v4_correct_seq/coupling_v4_status.txt
L=tandem_v4_correct_seq/coupling_v4.log
: > "$L"; echo "[$(date '+%F %H:%M:%S')] start" > "$S"

do_one(){ # tag traj topo
  local tag=$1 traj=$2 topo=$3
  echo "[$(date '+%F %H:%M:%S')] $tag transforms" >> "$S"
  $PY extract_cr2_transforms.py --traj "$traj" --topology "$topo" \
      --monomer tc_simple_old/classical_relaxed.pdb \
      --out "tandem_v4_correct_seq/tf_$tag.npz" >> "$L" 2>&1 || { echo "  transforms FAILED" >> "$S"; return 1; }
  echo "[$(date '+%F %H:%M:%S')] $tag coupling" >> "$S"
  $PY coupling_dcd_steom.py --density "$D" \
      --transforms "tandem_v4_correct_seq/tf_$tag.npz" \
      --out "tandem_v4_correct_seq/J_$tag.npz" --epsilon 1.77 >> "$L" 2>&1 \
      && echo "  $tag OK" >> "$S" || echo "  $tag FAILED" >> "$S"
}

do_one new tc_tandem_v4_correctseq/tandem_v4_1ns.dcd tc_tandem_v4_correctseq/solvated_protonated.pdb
do_one old tc_tandem_nvt_v2/tandem_nvt_v2_1000.dcd  tc_tandem_nvt_v2/solvated_protonated.pdb
echo "[$(date '+%F %H:%M:%S')] done" >> "$S"

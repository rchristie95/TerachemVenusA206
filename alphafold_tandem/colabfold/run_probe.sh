#!/usr/bin/env bash
# Probe the barrel-barrel register by perturbing the TRUNK, not the seed.
#
# The 14-seed AlphaFold3 sweep showed the seed does nothing (ANOVA p = 0.064,
# 0.930): the register is fixed upstream, by the sequence and its MSA. So the
# two perturbations that CAN move it are MSA depth and a different architecture.
#
#   A. venusdimer, AF2-MULTIMER: dock two Venus barrels with NO linker at all.
#      Independent architecture AND no tether to bias the arrangement. If AF2
#      finds the crystal interface, the crystal register is a real docking
#      minimum and AF3's tandem answer is something the linker/trunk imposes.
#   B. tandem511, AF2 monomer, MSA-depth ladder. Shallow MSAs expose alternative
#      modes that a full MSA collapses onto one.
set -u
cd "$(dirname "$0")"
export PATH="/home/robson/localcolabfold/colabfold-conda/bin:$PATH"
LOG=probe_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }
note "===== register probe start (pid $$) ====="
note "GPU $(nvidia-smi --query-gpu=name,memory.used --format=csv,noheader)"

note "--- A: AF2-multimer docking of two Venus barrels, no linker ---"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 5 --num-seeds 3 \
  --num-recycle 3 in/venusdimer.fasta out_dimer_multimer >> probe_run.log 2>&1
note "A done rc=$?  $(ls out_dimer_multimer/*.pdb 2>/dev/null | wc -l) pdb"

note "--- B: tandem511 monomer, MSA-depth ladder ---"
for MM in "512:1024" "64:128" "32:64" "16:32"; do
  tag=$(echo $MM | tr ':' '_')
  colabfold_batch --model-type alphafold2_ptm --num-models 5 --num-seeds 2 \
    --num-recycle 3 --max-msa "$MM" in/tandem511.fasta "out_tandem_msa$tag" >> probe_run.log 2>&1
  note "B max-msa $MM done rc=$?  $(ls out_tandem_msa$tag/*.pdb 2>/dev/null | wc -l) pdb"
done
note "===== register probe end ====="

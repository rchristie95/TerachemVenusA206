#!/bin/bash
# Candidate-2 detuning campaign: the kill test for the J=30 hypothesis.
#
# 40 frames (every other frame of the 80-frame full-system continuation),
# processed in a seeded-shuffle order so an early stop gives an unbiased
# time sample. Each frame: prepare (site A + site B QM/MM inputs, embedding
# policy identical to the v2 baseline: full-system, link-only exclusion,
# partner CR2 RESP, boundary charge conserved) then sequential TeraChem
# CAM-B3LYP/6-311+G** on the one GPU. Resumable: completed sites are skipped.

set -u
source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh

REPO=/home/robson/PetaChem
SED=$REPO/terachem_site_energy_cd
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
TC=/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem

# Every other frame of 0..79, shuffled with a fixed seed (python RNG --
# `shuf --random-source=<(yes SEED)` is degenerate, see runtime-gotchas).
FRAMES=$($PY -c "
import random
f = list(range(0, 80, 2))
random.Random(20260813).shuffle(f)
print(' '.join(map(str, f)))")

echo "frame order: $FRAMES"
for IDX in $FRAMES; do
    DIR=$SED/results/cand2_camb3lyp_frame_$(printf %04d "$IDX")
    if [ ! -d "$DIR" ]; then
        $PY "$SED/scripts/prepare_candidate2_frame.py" \
            --topology "$REPO/tc_candidate2_ff19sb_opc/solvated_protonated.pdb" \
            --trajectory "$REPO/tc_candidate2_ff19sb_opc/candidate2_fullsystem.dcd" \
            --frame-index "$IDX" \
            --output-dir "$DIR" \
            --embedding-cache "$SED/results/cand2_embedding_charges.npz" \
            --amber-cr2-prmtop "$REPO/anionic_build/monomer_solv.prmtop" \
            --full-system-embedding \
            --link-only-exclusion \
            --retain-partner-cr2-charges \
            --conserve-boundary-residue-charge || { echo "PREPARE FAILED frame $IDX"; exit 1; }
    fi
    $PY "$SED/scripts/launch_jobs.py" "$DIR" --terachem "$TC" --gpu 0 \
        || { echo "TERACHEM FAILED frame $IDX"; exit 1; }
    echo "=== frame $IDX complete $(date '+%H:%M:%S') ==="
done
echo "=== campaign complete $(date) ==="

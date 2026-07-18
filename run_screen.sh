#!/bin/bash
# TDDFT screen of candidate STEOM QM regions. Runs each on the GPU sequentially.
# Serial STEOM keeps running on CPU throughout (no contention).
source /home/robson/anaconda3/etc/profile.d/conda.sh
conda activate TeraChem
export TeraChem=/home/robson/Desktop/TeraChemPython/TeraChem
export PATH="$TeraChem/bin:$PATH"
export LD_LIBRARY_PATH="$TeraChem/lib:$LD_LIBRARY_PATH"

# wait for any terachem already on the GPU (the 54-atom run) to finish
while pgrep -f 'terachem tddft' >/dev/null; do sleep 10; done

for d in cr2only phenol tyrsc tyrfull hisonly phenolhis phenolser; do
  cd /home/robson/PetaChem/tc_screen_$d || continue
  sed -i 's/^cisnumstates .*/cisnumstates 8/' tddft.in   # enough roots to bracket the bright state
  echo "[screen $d start $(date +%T)]"
  terachem tddft.in > tddft.out 2>&1
  echo "[screen $d done  $(date +%T)  rc=$?]"
done
echo "SCREEN COMPLETE $(date +%T)"

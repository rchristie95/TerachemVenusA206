#!/bin/bash
# Check 2: TDDFT on the 44-atom region, standard (6-311G**) vs diffuse (6-311++G**),
# with transition-density output, to separate method / region / diffuse effects on J.
source /home/robson/anaconda3/etc/profile.d/conda.sh
conda activate TeraChem
export TeraChem=/home/robson/Desktop/TeraChemPython/TeraChem
export PATH="$TeraChem/bin:$PATH"; export LD_LIBRARY_PATH="$TeraChem/lib:$LD_LIBRARY_PATH"
for d in tc_tddft_44 tc_tddft_44_diff; do
  cd /home/robson/PetaChem/$d
  echo "[TDDFT $d start $(date '+%T')]"
  terachem td.in > td.out 2>&1
  echo "[TDDFT $d done $(date '+%T') rc=$?]"
  for f in scr/transdens_*.dx; do [ -f "$f" ] && cp -f "$f" .; done
  echo "  transdens files: $(ls transdens_*.dx 2>/dev/null | tr '\n' ' ')"
done
echo "TDD COMPLETE $(date '+%T')"

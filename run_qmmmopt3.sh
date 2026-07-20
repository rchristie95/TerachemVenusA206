#!/bin/bash
cd /home/robson/PetaChem/tc_qmmm_opt
source /home/robson/anaconda3/etc/profile.d/conda.sh; conda activate TeraChem
echo "[opt start $(date +%H:%M:%S)]"
timeout 100000 terachem opt.in > opt.out 2>&1; echo "[opt exit $? $(date +%H:%M:%S)]"
cp scr_opt/optim_geom.xyz qm_opt.xyz 2>/dev/null && echo "got qm_opt.xyz"

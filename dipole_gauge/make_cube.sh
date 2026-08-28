#!/usr/bin/env bash
# Item 2 (rigorous): transition density cube for the BRIGHT root 3 of the (entered as 5: TDA false => 2*state-1 selects the X+Y vector; entering 3 gives state 2 and yields cos=-0.60)
# gauge-verified CAM-B3LYP/def2-TZVPD run (TDA false, f_vel/f_len = 0.995).
set -u
cd "$(dirname "$0")" || exit 1
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
echo "[$(date '+%F %H:%M:%S')] orca_plot start" > cube_status.txt
printf '4\n100\n7\ny\n5\n11\n12\n' | $ORCA/orca_plot tddft_def2-TZVPD.gbw -i > orca_plot.log 2>&1
echo "[$(date '+%F %H:%M:%S')] orca_plot exit rc=$? " >> cube_status.txt
ls -la *.cube >> cube_status.txt 2>&1

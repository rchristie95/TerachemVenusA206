#!/bin/bash
cd /home/robson/PetaChem/neo_model/orca_steom
ORCABIN=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH="$ORCABIN:$LD_LIBRARY_PATH"
sed 's/%pal nprocs 16 end/%pal nprocs 1 end/' steom_svp.inp > steom_svp_serial.inp
echo "[orca STEOM-CCSD serial start $(date '+%H:%M:%S')]"
exec "$ORCABIN/orca" steom_svp_serial.inp

#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
# ORCA needs the matching OpenMPI on PATH and LD_LIBRARY_PATH, not just its own
# directory -- otherwise orca_startup_mpi dies on libmpi.so.40 immediately.
OMPI=/home/robson/anaconda3/envs/openmpi416
ORCADIR=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export PATH="$OMPI/bin:$ORCADIR:${PATH}"
export LD_LIBRARY_PATH="$OMPI/lib:$ORCADIR:${LD_LIBRARY_PATH:-}"
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a gauge_status.txt; }
note "===== TDDFT gauge diagnostic ====="
for B in def2-SVPD def2-TZVPD; do
  note "$B running"
  $ORCA tddft_${B}.inp > tddft_${B}.out 2>&1
  note "$B rc=$? terminated_normally=$(grep -c 'ORCA TERMINATED NORMALLY' tddft_${B}.out)"
done
note "===== done ====="

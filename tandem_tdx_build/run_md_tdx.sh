#!/usr/bin/env bash
# dVenus-TDX production NVT: 1 ns, 500 frames, 2000 MD steps between frames.
#
#   timestep                    1.0 fs
#   steps between saved frames  2000      -> 2 ps per frame
#   frames                      500
#   NVT steps                   1000000   -> 1.000 ns
#
# Those four numbers only close on a 1 fs timestep: 500 x 2000 x 1 fs = 1 ns.
# The v1-v4 TD trajectories ran at 2 fs, where the same frame count and frame
# stride would give 2 ns, so this run is NOT a drop-in timestep match to them --
# it matches their SAMPLING (500 frames over 1 ns is v4's 1000-frame ns
# thinned by two).  Set DT=2.0 STEPS_PER_FRAME=1000 for a 2 fs run with the
# same 500 frames and the same 1 ns.
#
# Structure: tandem_dimer_tdx_dark2.pdb, 508 residues, from
#   python3 tandem_tdx_build/build_tdx_structure.py \
#       --from-af3 alphafold_tandem/af3_tdx_server_run/fold_dvenus_tdx_dark2_model_0.cif
# which keeps tandem_dimer_v4_af3tails.pdb's barrels, linker and A206 crystal
# register and replaces barrel 2's mature CR2 with the Gly-Gly that Addgene
# #250939's deleted chromophore tyrosine (Y339del) leaves behind, taken from an
# AlphaFold model of the 510-aa TDX sequence.  NOT Gly-Cys-Gly: the deposit
# deletes the tyrosine, it does not mutate it to cysteine -- see
# alphafold_tandem/TDX_CONSTRUCT.md, which records the check against the deposit.
# The surviving chromophore and every atom outside the grafted pair are untouched
# to 0.000 A, so a TD-vs-TDX difference is attributable to the missing
# chromophore and not to a different starting geometry.
#
# The two junction bonds either side of the grafted pair come out at 1.48 and
# 1.58 A against an ideal 1.33 A: the crystal leaves a slot a three-residue
# triad occupied and the shortened strand spans 0.4 A more than it.  That strain
# is local and is what the 200 ps NPT equilibration below is there to take out;
# check the first frames rather than assuming it did.
#
# Seed and force field are deliberately those of the v4 TD run (amber14-all /
# tip3pfb, seed 20260618): TD and TDX are meant to differ in the construct only.
#
# The dark barrel holds no CR2, so this needs the run_nvt.py that requires at
# least one CR2 rather than exactly two, and --restrain-interface is not
# available (there is no second chromophore centroid to tether).
set -u

ROOT=${ROOT:-/home/robson/PetaChem}
PY=${PY:-/home/robson/anaconda3/envs/TeraChem/bin/python}
DARK=${DARK:-2}
DT=${DT:-1.0}
STEPS_PER_FRAME=${STEPS_PER_FRAME:-2000}
FRAMES=${FRAMES:-500}
NVT_STEPS=$((STEPS_PER_FRAME * FRAMES))
NPT_STEPS=${NPT_STEPS:-200000}          # 200 ps at 1 fs, matching v4's 200 ps

cd "$ROOT" || exit 1
TAG=tdx_dark${DARK}
S=tandem_tdx_build/md_${TAG}_status.txt
printf '[%s] %s MD start: %s steps x %s fs = %s ps, %s frames every %s steps\n' \
  "$(date '+%F %H:%M:%S')" "$TAG" "$NVT_STEPS" "$DT" \
  "$(awk -v n=$NVT_STEPS -v d=$DT 'BEGIN{printf "%.1f", n*d/1000}')" \
  "$FRAMES" "$STEPS_PER_FRAME" > "$S"

$PY -u run_nvt.py \
  --pdb tandem_tdx_build/tandem_dimer_${TAG}.pdb \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --workdir tc_tandem_${TAG} \
  --add-solvent --padding-a 10.0 --ionic-strength-m 0.15 \
  --npt-steps "$NPT_STEPS" \
  --timestep-fs "$DT" \
  --nvt-steps "$NVT_STEPS" \
  --openmm-trajectory-interval "$STEPS_PER_FRAME" \
  --openmm-trajectory-file tandem_${TAG}_1ns.dcd \
  --protein-topology-file protein_topology.pdb \
  --checkpoint-file nvt.chk --checkpoint-interval "$STEPS_PER_FRAME" \
  --seed 20260618 --no-video >> "tandem_tdx_build/md_${TAG}.log" 2>&1
rc=$?
echo "[$(date '+%F %H:%M:%S')] $TAG MD exit rc=$rc" >> "$S"
exit $rc

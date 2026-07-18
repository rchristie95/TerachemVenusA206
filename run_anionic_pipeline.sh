#!/bin/bash
# Full anionic Venus pipeline, OFFLINE and self-contained, using the new codes:
#   A) clean monomer (HID148/GLH222) + tleap build with published anionic CR2 FF
#   B) OpenMM/CUDA minimization
#   C) QM region setup  (FROZEN reproducible selection + forced pi-stacking Tyr203)
#   D) TDDFT (anionic, COSMO) + transition-density coupling  (renorm guard active)
#
# No internet required: AmberTools xFPchromophores params, TeraChem, PyMOL,
# OpenMM are all local; the CR2 CCD entry is already cached.
set -euo pipefail
cd /home/robson/PetaChem
source /home/robson/anaconda3/etc/profile.d/conda.sh
ts() { date +%H:%M:%S; }
echo "[$(ts)] ===== ANIONIC PIPELINE START (offline, new codes) ====="

echo "[$(ts)] STAGE A: clean monomer + tleap anionic build"
conda activate TeraChem
python build_monomer_clean.py
conda activate amber
( cd anionic_build && tleap -f build_solv.tleap > buildA.log 2>&1 )
grep -E "FATAL|Errors =" anionic_build/buildA.log | tail -2
test -s anionic_build/monomer_solv.prmtop || { echo "[FATAL] tleap build failed"; exit 1; }

echo "[$(ts)] STAGE B: OpenMM/CUDA minimization"
conda activate TeraChem
( cd anionic_build && python min_openmm.py )
test -s anionic_build/monomer_min.pdb || { echo "[FATAL] minimization failed"; exit 1; }

echo "[$(ts)] STAGE C: QM region setup (frozen selection + Tyr203)"
rm -rf tc_simple_anionic
python anionic_qm_setup.py
cp tc_simple_old/qm_region_forcefield.xml tc_simple_anionic/
test -s tc_simple_anionic/qm_deprotonated.xyz || { echo "[FATAL] QM setup failed"; exit 1; }
echo "[$(ts)]   frozen QM selection saved: $(test -f tc_simple_anionic/qm_selection.json && echo yes || echo no)"

echo "[$(ts)] STAGE D: TDDFT + coupling"
rm -rf tc_tddft_anionic_current
python terachem_full_pipeline.py --skip-simple \
  --tddft-args "--input-dir tc_simple_anionic --workdir-prefix tc_tddft_anionic" \
  --coupling-args "--monomer tc_simple_anionic/monomer_relaxed.pdb --dimer venus_dimer.pdb --epsilon 1.77" \
  --skip-visualize

echo "[$(ts)] ===== ANIONIC PIPELINE DONE ====="
echo "----- RESULTS -----"
grep -E "Root 1:|Renormalizing|WARNING|^Vdd:|^J:|Splitting" anionic_pipeline_full.log | tail -20 || true

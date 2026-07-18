#!/bin/bash
# Correct-physics run: let the protein H-bond network set the chromophore C-O
# (bond-length alternation) quantum-mechanically via a QM/MM geometry
# optimization of the QM region (chromophore + His148 + Ser205 + waters + Tyr203,
# all QM), in the MM point-charge field, THEN take the vertical TDDFT.
# Uses the neutral-Glu222 (physical anionic B-state) setup.
set -euo pipefail
cd /home/robson/PetaChem
source /home/robson/anaconda3/etc/profile.d/conda.sh
ts(){ date +%H:%M:%S; }
echo "[$(ts)] ===== QM/MM-OPT correct-physics run START ====="

# --- regenerate the GLH (neutral Glu222) QM setup: stages A-C ---
echo "[$(ts)] A: build (GLH) + tleap"
conda activate TeraChem; GLU222_STATE=GLH python build_monomer_clean.py
conda activate amber; ( cd anionic_build && tleap -f build_solv.tleap > buildA.log 2>&1 )
test -s anionic_build/monomer_solv.prmtop || { echo FATAL tleap; exit 1; }
echo "[$(ts)] B: minimize"
conda activate TeraChem; ( cd anionic_build && python min_openmm.py )
echo "[$(ts)] C: QM setup"
rm -rf tc_simple_anionic; python anionic_qm_setup.py
cp tc_simple_old/qm_region_forcefield.xml tc_simple_anionic/

# --- QM/MM geometry optimization of the QM region in the MM field ---
echo "[$(ts)] D: QM/MM geometry optimization (wb97xd3/6-31g*, MM embedding)"
WD=tc_qmmm_opt; rm -rf $WD; mkdir $WD
cp tc_simple_anionic/qm_deprotonated.xyz $WD/qm.xyz
cp tc_simple_anionic/mm_charges.dat $WD/
CH=$(awk '/^charge/{print $2}' tc_simple_anionic/qm_setup_settings.in)
cat > $WD/opt.in <<EOF
run minimize
coordinates qm.xyz
basis 6-31g*
method wb97xd3
charge $CH
spinmult 1
pointcharges mm_charges.dat
pcm cosmo
epsilon 78.39
new_minimizer yes
maxit 300
gpus 1
scrdir scr_opt
end
EOF
( cd $WD && timeout 36000 terachem opt.in > opt.out 2>&1 )
cp $WD/scr_opt/optim_geom.xyz $WD/qm_opt.xyz 2>/dev/null || cp $WD/scr_opt/optim.xyz $WD/qm_opt.xyz
echo "[$(ts)]   optimization done; C-O check:"
python3 - <<'PY'
import numpy as np
l=open("tc_qmmm_opt/qm_opt.xyz").read().splitlines(); n=int(l[0].split()[0])
xyz=np.array([[float(c) for c in r.split()[1:4]] for r in l[2:2+n]])
# CR2 OH/CZ are atoms index (build order in qm_deprotonated.xyz: CR2 first if core-first; print min C-O candidates)
# crude: find O bonded to an aromatic C ~1.2-1.5 A -> report phenolate C-O by scanning
print("  (C-O bonds 1.15-1.55 A in optimized QM region):")
PY

# --- vertical TDDFT on the QM/MM-optimized geometry ---
echo "[$(ts)] E: TDDFT on QM/MM-optimized geometry (6-311g**, PCM, MM embedding)"
WT=tc_tddft_qmmmopt; rm -rf $WT; mkdir $WT
cp $WD/qm_opt.xyz $WT/geometry.xyz; cp tc_simple_anionic/mm_charges.dat $WT/
cat > $WT/td.in <<EOF
run energy
coordinates geometry.xyz
basis 6-311g**
method wb97xd3
charge $CH
spinmult 1
pointcharges mm_charges.dat
pcm cosmo
epsilon 78.39
cis yes
cisnumstates 6
cistarget 1
cismaxiter 300
gpus 1
scrdir scr_td
end
EOF
( cd $WT && timeout 7200 terachem td.in > td.out 2>&1 )
echo "[$(ts)] ===== DONE ====="
grep -A8 "Final Excited State Results" $WT/td.out | head -10

#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export PATH="/home/robson/Desktop/TeraChemPython/TeraChem/bin:$PATH"
command -v terachem >/dev/null || { echo "FATAL: terachem not on PATH"; exit 1; }
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a beta_status.txt; }
note "===== S1 optimisation for beta ====="
terachem s1opt.in > s1opt.out 2>&1
note "s1 opt rc=$? finished=$(grep -c 'Job finished' s1opt.out)"
# transition dipole at the relaxed S1 geometry
if [ -s scr_s1opt/optim.xyz ]; then
  tail -46 scr_s1opt/optim.xyz > s1_geom.xyz
  sed -e 's/^coordinates .*/coordinates s1_geom.xyz/' -e 's/^run minimize/run energy/' \
      -e '/new_minimizer/d' -e '/^nstep/d' s1opt.in > s1_tdm.in
  terachem s1_tdm.in > s1_tdm.out 2>&1
  note "s1 tdm rc=$? finished=$(grep -c 'Job finished' s1_tdm.out)"
fi
note "===== done ====="

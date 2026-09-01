#!/usr/bin/env bash
# Is the 1-9 cm^-1 spread between the flag-on and flag-off runs caused by the
# flag, or is it just the run-to-run noise floor at 75,035 point charges?
# Rerun the identical flag-off input and compare it with itself.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence/selfint_check2 || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export OMP_NUM_THREADS=4
sed 's|^scrdir .*|scrdir scr_full_off2|' full_off.in > full_off2.in
timeout 45m /home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem full_off2.in > full_off2.out 2>&1
echo "rc=$?"

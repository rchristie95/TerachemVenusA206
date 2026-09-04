#!/usr/bin/env bash
# Low-frequency spectral density: stride-20 (80 fs) over the existing 8 ps /
# 4 fs trajectory, 100 frames, resolving 4.2 cm^-1 with Nyquist 208 cm^-1.
#
# WHY. v8dense (1 ps at 2 fs) cannot see below 33 cm^-1, and the quantum
# cumulant's two residuals -- the static width the band needs and the
# under-predicted peak-to-peak Stokes shift -- both point there. sigma cannot
# settle it: variance genuinely at 4 cm^-1 is reported at 33 by a 1 ps window,
# giving the same sigma but a very different lineshape, because 33 cm^-1
# motionally narrows where 4 cm^-1 is effectively static. Only a resolved
# spectrum in 4-33 cm^-1 separates them.
#
# The 33-208 cm^-1 overlap with v8dense is deliberate: the two spectral
# densities are cross-checked there before being spliced, rather than joined
# on faith.
#
# Settings MATCH the 99 v4full frames already computed -- validated profile,
# 7-root density-overlap tracking -- or the two halves would not be comparable.
# Frames are assigned round-robin, so each worker's own progress spans the full
# 8 ps and a partial run still yields a usable, sparser spectrum.
set -u
cd "$(dirname "$0")" || exit 1

if ! nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null; then
  echo "nvidia-smi FAILED -- rebuild and MOK-sign the RTX 4080 module first." >&2
  exit 1
fi
for w in a b c; do
  [ -f "manifests/v4full_lf${w}_frames.csv" ] || { echo "missing manifest lf$w" >&2; exit 1; }
done

systemd-run --user --unit=lf-campaign --same-dir /bin/bash -c '
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
cd /home/robson/PetaChem/tddft_ensemble_decoherence
mkdir -p logs/lf
echo "[$(date +%F\ %H:%M:%S)] launching three workers" | tee -a logs/lf/chain.log
export OMP_NUM_THREADS=4
for w in a b c; do
  "$PY" -u scripts/run_ensemble.py --window v4full --manifest-suffix "_lf$w" \
        --profile validated --roots 7 --gpu 0 >> logs/lf/lf_$w.log 2>&1 &
  sleep 20
done
wait
echo "[$(date +%F\ %H:%M:%S)] lf campaign complete" | tee -a logs/lf/chain.log
'
echo "launched as lf-campaign; follow with:"
echo "    tail -f logs/lf/lf_a.log"

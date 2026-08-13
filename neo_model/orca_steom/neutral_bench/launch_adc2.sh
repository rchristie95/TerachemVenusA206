#!/usr/bin/env bash
# ADC(2) alone, detached. TMPDIR is pinned to a visible directory inside the
# run folder rather than /tmp, because the previous attempt died writing pyscf
# integrals to /tmp when the root filesystem filled.
cd "$(dirname "$0")" || exit 1
export TMPDIR="$PWD/adc2_scratch"
export PYSCF_TMPDIR="$TMPDIR"
export OMP_NUM_THREADS=8
/home/robson/anaconda3/envs/adcc_env/bin/python adc2_neutral.py

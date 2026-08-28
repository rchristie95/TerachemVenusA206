#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
pilot_dir="$repo_root/tddft_ensemble_decoherence/reference_pilot"
output="$pilot_dir/reference_current.out"
terachem_bin="${TC_PATH:-/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem}"

if [[ -s "$output" ]] && rg -q 'Job finished:' "$output"; then
    echo "[skip] completed reference pilot: $output"
    exit 0
fi

cd "$pilot_dir"
start_epoch=$(date +%s)
CUDA_VISIBLE_DEVICES=0 "$terachem_bin" reference_current.in > reference_current.out 2>&1
end_epoch=$(date +%s)
printf '%s\n' "$((end_epoch - start_epoch))" > runtime_seconds.txt
rg -q 'Job finished:' reference_current.out
echo "[done] reference pilot runtime $((end_epoch - start_epoch)) s on GPU 0"

#!/usr/bin/env bash
# Resume everything that was running when the machine went down on 2026-09-03.
#
#   ./resume_all.sh          launch under systemd (survives the agent shell)
#   ./resume_all.sh --check  report state and exit, launch nothing
#
# Order matters: v8dense finishes on three workers, then the heavy basis rungs
# take the card ALONE. def2-TZVPD at 1292 basis functions died silently when it
# shared the GPU with three SVP workers.
set -u
cd "$(dirname "$0")" || exit 1
PY=/home/robson/anaconda3/envs/TeraChem/bin/python

count(){ "$PY" -c "
import csv,sys
c=t=f=0
for w in 'abc':
    try: r=list(csv.DictReader(open(f'manifests/$1_{w}_frames.csv')))
    except FileNotFoundError: continue
    c+=sum(x['convergence_status']=='complete' for x in r)
    f+=sum(x['convergence_status']=='failed' for x in r); t+=len(r)
print(f'{c}/{t} complete, {f} failed')"; }

echo "=== GPU"
if ! nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null; then
  echo "  nvidia-smi FAILED. If the kernel was updated, the RTX 4080 module needs"
  echo "  rebuilding and MOK-signing from /usr/src (580.126.09 open module)."
  echo "  Do NOT apt-install the skewed dkms package. Nothing below will run."
  exit 1
fi

echo "=== state"
for m in v8dense v7tdx v7tdxdark; do printf "  %-10s %s\n" "$m" "$(count $m)"; done
LAD=$("$PY" -c "
import json,glob,collections
r=[x for f in glob.glob('results/basis_ladder_*.json') for x in json.load(open(f))]
print(dict(collections.Counter(x['basis'] for x in r if x.get('status')=='complete')))")
echo "  ladder     $LAD"

echo "=== partial frame directories (block retry if present)"
find frames -maxdepth 2 -type d -name 'frame_*' 2>/dev/null | while read -r d; do
  [ -f "$d/preparation.json" ] || echo "  PARTIAL, delete before resuming: $d"
done
[ "${1:-}" = "--check" ] && exit 0

systemd-run --user --unit=resume-chain --same-dir /bin/bash -c '
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
cd /home/robson/PetaChem/tddft_ensemble_decoherence
mkdir -p logs/resume
note(){ echo "[$(date +%F\ %H:%M:%S)] $*" | tee -a logs/resume/chain.log; }

note "finishing v8dense on three workers"
export OMP_NUM_THREADS=4
for w in a b c; do
  "$PY" -u scripts/run_ensemble.py --window v8dense --manifest-suffix "_$w" \
        --profile svp --roots 7 --gpu 0 >> logs/resume/v8dense_$w.log 2>&1 &
  sleep 20
done
wait
note "v8dense done"

note "def2-TZVPD, one worker, card to itself"
export OMP_NUM_THREADS=16
"$PY" -u scripts/run_basis_ladder.py 0 a 1 tzvpd >> logs/resume/tzvpd.log 2>&1
note "TZVPD rc=$?"
"$PY" -u scripts/run_basis_ladder.py 0 a 1 qzvpd >> logs/resume/qzvpd.log 2>&1
note "QZVPD rc=$?"
note "===== resume chain complete ====="
'
echo "=== launched as resume-chain; follow with:"
echo "    journalctl --user -u resume-chain -f    or    tail -f logs/resume/chain.log"

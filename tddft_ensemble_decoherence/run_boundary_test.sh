#!/usr/bin/env bash
# 5-frame QM/MM boundary-treatment test.
#
# QUESTION. The production embedding (window v4full) cuts three bonds per site,
# caps with link H, and then simply DELETES the three severed-bond partners from
# the point-charge field with no charge redistribution. Two of the three cuts go
# through peptide bonds, so deleting LEU:65:C (~+0.6 e) while retaining its
# carbonyl O (~-0.57 e, ~2.4 A away) leaves an orphaned amide dipole right at the
# QM boundary. Does any of that actually move the observable, Delta = E_A - E_B?
#
# Prior expectation: probably not much, because the boundary is built identically
# at both sites (MM net charge agrees to 13 decimals), so symmetric errors cancel
# in a DIFFERENCE. This measures that rather than assuming it.
#
# VARIANTS (QM region identical in all three -- 44 atoms -- so the density-overlap
# state reference stays valid):
#   v4full      baseline: link-only exclusion, charge deleted        [already run]
#   v4bnd_rc    link-only exclusion + charge redistributed over residue
#   v4bnd_excl  geometric 3.0 A exclusion (removes the whole orphaned
#               amide) + charge redistributed
#
# Not tested here: moving the backbone cuts to C-alpha, and putting the first
# solvation shell in QM. Both change the QM atom count, which invalidates the
# transition-density state reference -- a larger job, not a flag.
set -u
cd "$(dirname "$0")" || exit 1
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
PY=/home/robson/anaconda3/envs/TeraChem/bin/python
STATUS=boundary_test_status.txt
note() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$STATUS"; }

note "===== boundary test start (pid $$) ====="
note "waiting for spread50 to release the GPU"
while systemctl --user is-active spread50 >/dev/null 2>&1; do sleep 60; done
while pgrep -x terachem >/dev/null 2>&1; do sleep 60; done
note "GPU free; starting"

for w in v4bnd_rc v4bnd_excl; do
  while read -r f; do
    s=$(printf "frames/%s/frame_%04d/frame_summary.json" "$w" "$f")
    [ -f "$s" ] && { note "$w frame $f already done"; continue; }
    t0=$(date +%s)
    timeout 2h "$PY" scripts/run_ensemble.py --window "$w" --start "$f" --stop $((f+1)) --gpu 0 \
        >> boundary_test.log 2>&1
    rc=$?
    if [ -f "$s" ]; then note "$w frame $f OK ($(( $(date +%s)-t0 ))s)"
    else note "$w frame $f FAILED rc=$rc"; tail -3 boundary_test.log | sed 's/^/      /' | tee -a "$STATUS"; fi
  done < boundary_test_frames.txt
done

note "===== comparison ====="
"$PY" - <<'PYEOF' 2>&1 | tee -a "$STATUS"
import json, glob, statistics as st
def gaps(w):
    out={}
    for p in glob.glob(f'frames/{w}/frame_*/frame_summary.json'):
        d=json.load(open(p)); A,B=d.get('site_A'),d.get('site_B')
        if A and B and A.get('status')=='selected' and B.get('status')=='selected':
            out[d['frame_index']]=abs(A['energy_cm-1']-B['energy_cm-1'])
    return out
base, rc, ex = gaps('v4full'), gaps('v4bnd_rc'), gaps('v4bnd_excl')
frames=[int(x) for x in open('boundary_test_frames.txt').read().split()]
print(f"\n{'frame':>7}{'baseline':>11}{'+redistrib':>12}{'delta':>9}{'3.0A excl':>12}{'delta':>9}")
print("-"*60)
d1,d2=[],[]
for f in frames:
    b=base.get(f)
    if b is None: continue
    r,e=rc.get(f),ex.get(f)
    s=f"{f:>7}{b:>11.1f}"
    if r is not None: s+=f"{r:>12.1f}{r-b:>+9.1f}"; d1.append(r-b)
    else: s+=f"{'-':>12}{'-':>9}"
    if e is not None: s+=f"{e:>12.1f}{e-b:>+9.1f}"; d2.append(e-b)
    else: s+=f"{'-':>12}{'-':>9}"
    print(s)
print()
if d1: print(f"  charge redistribution : mean shift {st.mean(d1):+7.1f} cm^-1, max |shift| {max(map(abs,d1)):6.1f}")
if d2: print(f"  3.0 A exclusion       : mean shift {st.mean(d2):+7.1f} cm^-1, max |shift| {max(map(abs,d2)):6.1f}")
print(f"\n  for scale: <|Delta|> ~ 574-631 cm^-1, ensemble s.e. ~ +/-54-72 cm^-1")
print("  a shift small against that s.e. means the boundary treatment does not")
print("  matter for Delta, and effort belongs on polarisation / QM region size.")
PYEOF
note "===== boundary test complete ====="

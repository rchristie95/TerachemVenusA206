#!/usr/bin/env bash
# Does `pointcharges_self_interaction no` change an EXCITATION energy?
#
# Physically it cannot: the MM-MM Coulomb sum contains no wavefunction, so it is
# a constant offset on every state of a fixed field. This checks TeraChem's
# implementation agrees, because the entire 500-frame campaign runs with it off.
#
# Two passes. Pass 1 uses a 2500-charge truncation of the real frame-900 field:
# both runs see the IDENTICAL field and differ only in the flag, so it isolates
# the flag's semantics and is cheap. Pass 2 repeats it at the full 75,035
# charges; an earlier attempt there dumped core, which is worth pinning down
# but does not affect the campaign either way.
set -u
cd /home/robson/PetaChem/tddft_ensemble_decoherence || exit 1
# SetTCVars.sh does not export TC_EXE outside an interactive shell, so the binary
# is named by absolute path here, the same way run_ensemble.py does it.
set +u; source /home/robson/Desktop/TeraChemPython/TeraChem/SetTCVars.sh; set -u
export OMP_NUM_THREADS=4
SRC=steom_v5/frame_0900/site_A
D=selfint_check2; rm -rf $D; mkdir -p $D
cp $SRC/geometry.xyz $D/

python3 - <<'PY'
import pathlib
src = pathlib.Path("steom_v5/frame_0900/site_A/mm_charges.dat")
d = pathlib.Path("selfint_check2")
lines = src.read_text().splitlines()
n = int(lines[0].split()[0]); body = lines[2:2+n]
(d/"mm_full.dat").write_text(f"{n}\nfull field\n" + "\n".join(body) + "\n")
k = 2500
(d/"mm_small.dat").write_text(f"{k}\nfirst {k} charges of the frame-900 field\n" + "\n".join(body[:k]) + "\n")
print(f"wrote full={n} small={k}")
PY

mkin() {  # mkin <name> <chargefile> <yes|no>
cat > "$D/$1.in" <<IN
coordinates geometry.xyz
run energy
cis yes
basis def2-svp
method wb97xd3
charge -1
spinmult 1
pointcharges $2
pointcharges_self_interaction $3
scrdir scr_$1
cisnumstates 7
cismaxiter 200
cismax 500
cistransdensity no
cisdiffdensity no
cistarget 7
scf diis+a
threall 1.0e-13
gpus 1
end
IN
}
mkin small_on  mm_small.dat yes
mkin small_off mm_small.dat no
mkin full_on   mm_full.dat  yes
mkin full_off  mm_full.dat  no

cd $D
for j in small_on small_off full_on full_off; do
  echo "[$(date '+%F %H:%M:%S')] $j start"
  timeout 45m /home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem $j.in > $j.out 2>&1
  echo "[$(date '+%F %H:%M:%S')] $j exit rc=$?"
done
cd ..

python3 - <<'PY'
import re, pathlib
HARTREE_CM = 219474.6313632
def read(p):
    t = pathlib.Path(p).read_text(errors="replace")
    tot = None
    for m in re.finditer(r"FINAL ENERGY:\s+(-?\d+\.\d+)", t): tot = float(m.group(1))
    ex = [float(x) for x in re.findall(r"^\s*1\s+(\d+\.\d+)\s+", t, re.M)]
    st = re.findall(r"Root\s+1\s+.*?(\d+\.\d+)\s+eV", t)
    return tot, t
def states(p):
    t = pathlib.Path(p).read_text(errors="replace")
    rows = re.findall(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+", t, re.M)
    return t
print("=== self-interaction flag verification, frame 900 site A ===")
for tag in ["small", "full"]:
    try:
        on = pathlib.Path(f"selfint_check2/{tag}_on.out").read_text(errors="replace")
        off = pathlib.Path(f"selfint_check2/{tag}_off.out").read_text(errors="replace")
    except OSError as e:
        print(f"  {tag}: missing output ({e})"); continue
    def final(t):
        m = re.findall(r"FINAL ENERGY:\s+(-?\d+\.\d+)", t)
        return float(m[-1]) if m else None
    def exc(t):
        # TeraChem 'Root   N   E(a.u.)   E(eV)   osc' block
        m = re.findall(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$", t, re.M)
        return [(int(a), float(c), float(d)) for a, b, c, d in m]
    fo, ff = final(on), final(off)
    eo, ef = exc(on), exc(off)
    print(f"  --- {tag} field ---")
    print(f"    total energy  ON {fo}   OFF {ff}   diff {None if None in (fo,ff) else f'{fo-ff:+.8f} au'}")
    if eo and ef and len(eo) == len(ef):
        d = [(a[1]-b[1])*8065.544 for a, b in zip(eo, ef)]
        print(f"    excitation energies, ON-OFF (cm^-1): "
              + " ".join(f"{x:+.4f}" for x in d[:7]))
        print(f"    max |delta| = {max(abs(x) for x in d):.5f} cm^-1  -> "
              + ("SAFE" if max(abs(x) for x in d) < 0.5 else "NOT SAFE"))
    else:
        print(f"    could not pair excitation blocks (ON {len(eo)} rows, OFF {len(ef)} rows)")
PY

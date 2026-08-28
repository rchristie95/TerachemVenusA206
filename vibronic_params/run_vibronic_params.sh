#!/usr/bin/env bash
# Stage 1: gas-phase GS opt+freq.  Stage 2: ESD vertical gradient -> Huang-Rhys.
# Restartable: stage 1 resumes from the last trajectory geometry, but only after
# validating it (the first attempt's geometry was destroyed by MM-charge
# collapse, so a blind restart would have resumed from garbage).
set -u
cd "$(dirname "$0")" || exit 1
ORCA=/home/robson/PetaChem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH=$ORCA:/home/robson/anaconda3/envs/openmpi416/lib:${LD_LIBRARY_PATH:-}
export PATH=/home/robson/anaconda3/envs/openmpi416/bin:$PATH
PY=/home/robson/anaconda3/bin/python3
LOG=vibronic_params_status.txt
note(){ echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$LOG"; }

avail=$(free -g | awk 'NR==2{print $7}')
if [ "$avail" -lt 40 ]; then note "REFUSING: only ${avail} GB RAM"; exit 1; fi

# --- restart guard: accept vib_gs.xyz only if it is physically sane ----------
if [ -f vib_gs.xyz ] && [ ! -f vib_gs.hess ]; then
  if $PY - <<'CHK'
import numpy as np,sys
def rd(f):
    L=open(f).read().split('\n'); n=int(L[0])
    return np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2+n]])
try: a,b=rd('geom_cthrp.xyz'),rd('vib_gs.xyz')
except Exception: sys.exit(1)
if a.shape!=b.shape: sys.exit(1)
d=np.linalg.norm(a-b,axis=1)
D=np.linalg.norm(b[:,None,:]-b[None,:,:],axis=2)+np.eye(len(b))*9
sys.exit(0 if (d.max()<0.8 and D.min()>0.6) else 1)
CHK
  then note "RESTART from vib_gs.xyz (validated)"; cp vib_gs.xyz restart_geom.xyz
       sed -i 's|geom_cthrp.xyz|restart_geom.xyz|' vib_gs.inp
  else note "REJECT vib_gs.xyz (unphysical) -- restarting from geom_cthrp.xyz"
       rm -f vib_gs.xyz; fi
fi

note "===== vibronic parameter run start (pid $$) ====="
note "RAM ${avail} GB free; gas-phase optimisation (field.pc deliberately absent)"
if ! grep -q 'ORCA TERMINATED NORMALLY' vib_gs.out 2>/dev/null; then
  note "START   stage 1: gas-phase GS opt + freq"
  timeout 30h "$ORCA/orca" vib_gs.inp > vib_gs.out 2>&1
  if grep -q 'ORCA TERMINATED NORMALLY' vib_gs.out; then
    nimag=$(grep -c '\*\*\*imaginary mode\*\*\*' vib_gs.out || true)
    note "OK      stage 1 (${nimag} imaginary modes)"
    [ "${nimag:-0}" -gt 0 ] && note "WARNING imaginary modes -> Huang-Rhys factors unreliable"
  else note "FAIL    stage 1"; exit 1; fi
fi
cat > vib_esd.inp <<'INP'
! CAM-B3LYP def2-SVP def2/J ESD(ABS) TightSCF RIJCOSX
%tddft NRoots 5 IRoot 1 end
%esd GSHESSIAN "vib_gs.hess"
     HESSFLAG VG
end
%maxcore 8000
%pal nprocs 4 end
* xyzfile -1 1 vib_gs.xyz
INP
note "START   stage 2: ESD vertical gradient -> Huang-Rhys per mode"
timeout 20h "$ORCA/orca" vib_esd.inp > vib_esd.out 2>&1
grep -q 'ORCA TERMINATED NORMALLY' vib_esd.out && note "OK      stage 2" || note "FAIL    stage 2"
note "===== complete ====="

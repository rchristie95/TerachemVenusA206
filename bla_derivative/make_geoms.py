#!/usr/bin/env python3
"""Displace along dimensionless normal coordinates q to get d = d(delta_mu)/dq.

For dimensionless q the Cartesian displacement is
    dx_i[bohr] = q * L_i * sqrt(1/omega_au) / sqrt(m_i[m_e])
with L the orthonormal MASS-WEIGHTED eigenvector.

Everything is done at the FREQUENCY geometry, gas phase, camb3lyp/6-31gs -- the
same conditions the modes were computed in -- so |d| and |delta_mu_0| come from
one consistent point and their RATIO (which sets the HT fraction) is meaningful.
"""
import numpy as np, re
from pathlib import Path
SRC=Path('/home/robson/PetaChem/vibronic_params/tc')
BOHR=0.529177210903; AMU=1822.888486
txt=(SRC/'scr_freq/Mass.weighted.modes.dat').read_text().split('\n')
modes={}; cur=None; buf=[]
for l in txt:
    m=re.match(r'=== Mode (\d+): ([\d.]+) cm\^-1 ===', l.strip())
    if m:
        if cur: modes[cur]=np.array(buf,float)
        cur=(int(m.group(1)),float(m.group(2))); buf=[]
    elif l.strip():
        try: buf.append(float(l))
        except ValueError: pass
if cur: modes[cur]=np.array(buf,float)
masses=[float(l.split()[1]) for l in (SRC/'scr_freq/Atom.mass.dat').read_text().split('\n')[1:] if l.strip()]
masses=np.array(masses)
g=(SRC/'optim.xyz').read_text().split('\n')
n=int(g[0].split()[0]); els=[l.split()[0] for l in g[2:2+n]]
xyz=np.array([[float(x) for x in l.split()[1:4]] for l in g[2:2+n]])
print(f'{len(modes)} modes, {n} atoms, masses {masses[:3]}...')
TARGETS=[1226.9,1675.2,1389.6,1091.2]
sel=[]
for (i,w),L in modes.items():
    for t in TARGETS:
        if abs(w-t)<1.0: sel.append((i,w,L))
print('selected modes:', [(i,round(w,1)) for i,w,_ in sel])
for i,w,L in sel:
    L=L/np.linalg.norm(L)
    wau=w/219474.6314
    amp=np.sqrt(1.0/wau)                       # mass-weighted amplitude for q=1
    disp=(L.reshape(-1,3)*amp)/np.sqrt(masses*AMU)[:,None]*BOHR   # Angstrom
    for q in (-0.5,0.5):
        new=xyz+q*disp
        with open(f'mode{i}_q{q:+.1f}.xyz','w') as f:
            f.write(f'{n}\nmode {i} {w:.1f} cm-1, dimensionless q={q:+.2f}\n')
            for e,c in zip(els,new): f.write(f'{e} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n')
    print(f'  mode {i} ({w:.1f} cm^-1): max Cartesian displacement at q=0.5 = {np.abs(0.5*disp).max():.4f} A')
import shutil; shutil.copy(SRC/'optim.xyz','ref_q0.xyz')

#!/usr/bin/env python3
"""ADC(2) transition dipoles via adcc, on the SAME geometry + MM field as the
STEOM and TDDFT runs in this directory, so the numbers are directly comparable.

This exists because the transition dipole is the critical path: J scales as
|mu|^2, and STEOM (9.80 D) vs TDDFT (7.84 D) is (9.80/7.84)^2 = 1.56x, which
is essentially the whole 1.61x discrepancy Cusick 2026 reports. STEOM's dipole
is basis-converged (SVP 9.797 -> SVPD 9.744, -0.54%) while TDDFT's is not
(SVP 8.707 -> TZVPD 7.839, -10%), so the gap is a METHOD gap. ADC(2) is the
independent third opinion that breaks the tie.

NB: use weekend_dipole/geom_cthrp.xyz (relaxed C-O 1.263), NOT
neo_model/qm_cthrp.xyz (MM-min geom) -- they differ by 0.080 A RMSD and are
not interchangeable. The MM fields ARE identical between the two locations.

Usage: run_adcc_dipole.py <qm.xyz> <mm.dat|gas> <basis> <n_frozen_core> <n_frozen_virtual> <nstates>
Env:   adcc_env
"""
import os, sys, time
import numpy as np
import adcc
from pyscf import gto, scf, qmmm

t0 = time.time()
geom, mmf, basis = sys.argv[1], sys.argv[2], sys.argv[3]
nfc = int(sys.argv[4]); nfv = int(sys.argv[5]); nst = int(sys.argv[6])
mmf = None if mmf == 'gas' else mmf
AU2D, EV, NM = 2.5417464519, 27.211386245988, 1239.841984

mol = gto.M(atom=geom, basis=basis, charge=-1, spin=0, max_memory=20000, verbose=3)
mf = scf.RHF(mol).density_fit()
if mmf:
    D = open(mmf).readlines(); nmm = int(D[0].split()[0])
    mm = np.array([[float(x) for x in D[2 + i].split()] for i in range(nmm)])
    mf = qmmm.mm_charge(mf, mm[:, 1:4], mm[:, 0], unit='Angstrom')
    print(f"MM: {nmm} charges, net {mm[:,0].sum():+.4f} e", flush=True)
mf.conv_tol = 1e-10
mf.kernel()
if not mf.converged:
    mf = mf.newton(); mf.kernel()
print(f"[{time.time()-t0:.0f}s] RHF E={mf.e_tot:.8f} conv={mf.converged} "
      f"nao={mol.nao} basis={basis}", flush=True)

# adcc's default "standard" allocator is unbounded: def2-svp with all virtuals
# peaked at 33.9 GB, and with 250 active virtuals it still climbed to the 44 GB
# cgroup cap, dying ~11 min into ADC(2) (the SCF takes 24 s -- it is the tensor
# pool, not the SCF). "libxm" maps tensor blocks to disk pagefiles instead, so
# RAM stays bounded and the NVMe absorbs the rest. Scratch must NOT be tmpfs.
# libxm (out-of-core) bounded RAM but was disk-bound: def2-svp burned the full
# 8 h cap without converging. Back to the in-core allocator with a smaller
# problem instead; the cgroup cap is the safety net.

adcc.set_n_threads(28)
kw = dict(n_singlets=nst, conv_tol=1e-5)
if nfc: kw['frozen_core'] = nfc
if nfv: kw['frozen_virtual'] = nfv
print(f"adcc ADC(2) basis={basis} frozen_core={nfc} frozen_virtual={nfv} "
      f"n_singlets={nst} threads={adcc.get_n_threads()}", flush=True)
state = adcc.adc2(mf, **kw)
import resource
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576
print(f"[{time.time()-t0:.0f}s] ADC(2) converged   peak RSS {peak:.1f} GB", flush=True)
print(state.describe())

print(f"\n===== ADC(2)/{basis} TRANSITION DIPOLES =====")
print(f"{'state':6} {'E(eV)':>8} {'nm':>8} {'f':>9} {'|mu|(D)':>9}   {'DX':>8} {'DY':>8} {'DZ':>8}")
best = (0.0, None)
for i in range(len(state.excitation_energy)):
    e = float(state.excitation_energy[i])
    f = float(state.oscillator_strength[i])
    tdm = np.asarray(state.transition_dipole_moment[i], dtype=float)
    mu = np.linalg.norm(tdm) * AU2D
    print(f"S{i+1:<5} {e*EV:8.4f} {NM/(e*EV):8.1f} {f:9.5f} {mu:9.3f}   "
          f"{tdm[0]:8.4f} {tdm[1]:8.4f} {tdm[2]:8.4f}", flush=True)
    if f > best[0]:
        best = (f, (i + 1, e * EV, NM / (e * EV), f, mu))
if best[1]:
    i, ee, wl, f, mu = best[1]
    print(f"\nBRIGHT STATE: S{i}  {ee:.4f} eV  {wl:.1f} nm  f={f:.5f}  |mu|={mu:.3f} D")
    print(f"  reference on this geometry/field:")
    print(f"    STEOM/SVP   501.5 nm  9.797 D      STEOM/SVPD  532.8 nm  9.744 D")
    print(f"    TDDFT/SVP   385.5 nm  8.707 D      TDDFT/SVPD  399.3 nm  8.239 D")
    print(f"    TDDFT/TZVPD 400.9 nm  7.839 D      experiment  ~515 nm")
print(f"[{time.time()-t0:.0f}s] ALL DONE", flush=True)


#!/usr/bin/env python3
"""
Test the user's hypothesis: the NEO-optimized geometry (quantum bridging proton)
has a longer phenolate C-O than the classical-optimized geometry, red-shifting
the absorption.

Usage: run_opt.py cl   |   run_opt.py neo
  cl  : classical geometry optimization (all nuclei classical)
  neo : NEO geometry optimization (His148 HD1 quantum, no epc)
Both then take a CLASSICAL vertical TDDFT on their OWN relaxed geometry, so the
two excitation energies differ only by the geometry (the NQE structural effect).
Logs C-O / N-H / H..O / N..O per optimization cycle to opt_<mode>.dat.
"""
import json, sys, numpy as np
mode = sys.argv[1]
idx = json.load(open("idx_full.json"))
HD1, O, CZ, ND1 = idx["HD1"], idx["O"], idx["CZ"], idx["ND1"]
XC, EBAS, NBAS = "b3lyp", "6-31g*", "pb4d"
BOHR = 0.5291772109; EV = 27.211386; NM = 1239.841984
L = open("chromo_his.xyz").read().splitlines(); n = int(L[0].split()[0])
els = [r.split()[0] for r in L[2:2+n]]
atomstr = "\n".join(L[2:2+n])

def cb_factory(tag):
    fh = open(f"opt_{tag}.dat", "w")
    def cb(envs):
        c = envs["coords"].reshape(-1, 3) * BOHR
        co = np.linalg.norm(c[CZ]-c[O]); nh = np.linalg.norm(c[ND1]-c[HD1])
        ho = np.linalg.norm(c[HD1]-c[O]); no = np.linalg.norm(c[ND1]-c[O])
        fh.write(f"cyc {envs['self'].cycle:3d}  CO={co:.4f}  NH={nh:.4f}  "
                 f"HO={ho:.4f}  NO={no:.4f}  E={envs['energy']:.6f}\n"); fh.flush()
    return cb

from pyscf.geomopt.geometric_solver import optimize
from pyscf import gto, dft, tddft

if mode == "cl":
    m = gto.M(atom=atomstr, basis=EBAS, charge=-1, spin=0, verbose=0, max_memory=8000)
    mf = dft.RKS(m); mf.xc = XC
    mopt = optimize(mf, maxsteps=80, callback=cb_factory("cl"))
    coords = mopt.atom_coords() * BOHR
elif mode in ("neo", "neoepc"):
    from pyscf import neo
    nm = neo.M(atom=atomstr, basis=EBAS, nuc_basis=NBAS, charge=-1, spin=0,
               quantum_nuc=[HD1], verbose=0, max_memory=8000)
    nmf = neo.KS(nm, xc=XC, epc=("17-2" if mode == "neoepc" else None))
    mopt = optimize(nmf, maxsteps=80, callback=cb_factory(mode))
    coords = mopt.atom_coords() * BOHR
else:
    raise SystemExit("mode must be cl|neo|neoepc")

# save optimized geometry
with open(f"geom_{mode}.xyz", "w") as f:
    f.write(f"{n}\n{mode} optimized; b3lyp/6-31g*\n")
    for el, c in zip(els, coords):
        f.write(f"{el} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n")
co = np.linalg.norm(coords[CZ]-coords[O]); ho = np.linalg.norm(coords[HD1]-coords[O])
no = np.linalg.norm(coords[ND1]-coords[O]); nh = np.linalg.norm(coords[ND1]-coords[HD1])
print(f"[{mode}] OPTIMIZED: C-O={co:.4f}  N-H={nh:.4f}  H..O={ho:.4f}  N..O={no:.4f}", flush=True)

# classical vertical TDDFT on the optimized geometry
gm = gto.M(atom="\n".join(f"{el} {c[0]} {c[1]} {c[2]}" for el, c in zip(els, coords)),
           basis=EBAS, charge=-1, spin=0, verbose=0, max_memory=8000)
mf2 = dft.RKS(gm); mf2.xc = XC; mf2.kernel()
td = tddft.TDDFT(mf2); td.nstates = 8; td.kernel()
osc = td.oscillator_strength(); e = np.asarray(td.e)*EV; i = int(np.argmax(osc))
print(f"[{mode}] TDDFT bright pi->pi*: {e[i]:.3f} eV = {NM/e[i]:.1f} nm (f={osc[i]:.3f})", flush=True)
print(f"[{mode}] RESULT  C-O={co:.4f} A   lambda={NM/e[i]:.1f} nm", flush=True)

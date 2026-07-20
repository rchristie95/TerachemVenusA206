#!/usr/bin/env python3
"""Build a Gaussian .cube of the STEOM-CCSD bright-state transition density
rho_0n(r) ~ phi_hole(r) * phi_particle(r) (dominant S1 NTO pair, mo92 x mo93),
preserving the QM-region atoms so PyMOL shows the chromophore + the +/- lobes."""
import numpy as np

def read_cube_full(fn):
    with open(fn) as f:
        c1, c2 = f.readline(), f.readline()
        t = f.readline().split(); natom = int(t[0]); origin = list(map(float, t[1:4]))
        axes = [f.readline() for _ in range(3)]
        na = abs(natom)
        atoms = [f.readline() for _ in range(na)]
        extra = f.readline() if natom < 0 else ""   # MO cube: orbital-index line
        n = [int(a.split()[0]) for a in axes]
        data = np.fromstring(" ".join(f.read().split()), sep=" ").reshape(n)
    return dict(c1=c1, c2=c2, natom=na, origin=origin, axes=axes, atoms=atoms, n=n, data=data)

base = "/home/robson/PetaChem/neo_model/orca_steom/steom_phenol_svpd.s1"
hole = read_cube_full(base + ".mo92a.cube")
part = read_cube_full(base + ".mo93a.cube")

rho = hole["data"] * part["data"]
rho = rho / np.abs(rho).max()        # normalize so isovalue ~0.02-0.06 works
print(f"grid {hole['n']}  rho range [{rho.min():+.3f}, {rho.max():+.3f}]  "
      f"+lobe vox {(rho>0.04).sum()}  -lobe vox {(rho<-0.04).sum()}")

out = "/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.cube"
with open(out, "w") as f:
    f.write("STEOM-CCSD S1 transition density (phi_hole * phi_particle)\n")
    f.write("rho_0n(r), normalized to max|rho|=1\n")
    f.write(f"{hole['natom']:5d}{hole['origin'][0]:12.6f}{hole['origin'][1]:12.6f}{hole['origin'][2]:12.6f}\n")
    for a in hole["axes"]:  f.write(a if a.endswith("\n") else a + "\n")
    for a in hole["atoms"]: f.write(a if a.endswith("\n") else a + "\n")
    flat = rho.flatten()                 # z-fastest (cube order, same as input reshape)
    for i in range(0, len(flat), 6):
        f.write("".join(f"{v:13.5E}" for v in flat[i:i+6]) + "\n")
print("wrote", out)

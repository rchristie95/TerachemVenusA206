#!/usr/bin/env python3
"""Build STEOM-CCSD S1 density cubes for the paper-style figures, from the NTO pair:
  transition density rho_0n  ~  phi_hole * phi_particle      (red +/ blue -)
  difference density rho_e-rho_g ~ |phi_particle|^2 - |phi_hole|^2  (green gain / violet loss)
Both on the same grid/atoms as the input cubes, each normalized to max|.|=1."""
import numpy as np

def read_cube_full(fn):
    with open(fn) as f:
        c1, c2 = f.readline(), f.readline()
        t = f.readline().split(); natom = int(t[0]); origin = list(map(float, t[1:4]))
        axes = [f.readline() for _ in range(3)]
        na = abs(natom)
        atoms = [f.readline() for _ in range(na)]
        if natom < 0: f.readline()                       # MO cube orbital-index line
        n = [int(a.split()[0]) for a in axes]
        data = np.fromstring(" ".join(f.read().split()), sep=" ").reshape(n)
    return dict(natom=na, origin=origin, axes=axes, atoms=atoms, n=n, data=data)

base = "/home/robson/PetaChem/neo_model/orca_steom/steom_phenol_svpd.s1"
hole = read_cube_full(base + ".mo92a.cube")
part = read_cube_full(base + ".mo93a.cube")

def write_cube(out, ref, field, title):
    field = field / np.abs(field).max()
    with open(out, "w") as f:
        f.write(title + "\n"); f.write("normalized to max|f|=1\n")
        f.write(f"{ref['natom']:5d}{ref['origin'][0]:12.6f}{ref['origin'][1]:12.6f}{ref['origin'][2]:12.6f}\n")
        for a in ref["axes"]:  f.write(a if a.endswith("\n") else a + "\n")
        for a in ref["atoms"]: f.write(a if a.endswith("\n") else a + "\n")
        flat = field.flatten()
        for i in range(0, len(flat), 6):
            f.write("".join(f"{v:13.5E}" for v in flat[i:i+6]) + "\n")
    print(f"wrote {out}  (+vox {(field>0.04).sum()}, -vox {(field<-0.04).sum()})")

trans = hole["data"] * part["data"]
diff  = part["data"]**2 - hole["data"]**2
write_cube("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens.cube", hole, trans,
           "STEOM-CCSD S1 transition density (phi_hole * phi_particle)")
write_cube("/home/robson/PetaChem/neo_model/orca_steom/steom_diffdens.cube", hole, diff,
           "STEOM-CCSD S1 difference density (|particle|^2 - |hole|^2)")

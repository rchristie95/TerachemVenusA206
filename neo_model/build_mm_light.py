#!/usr/bin/env python3
"""Lightweight MM field: whole residues within R A of the QM region (drops the
~28k bulk-solvent charges that cost time and contribute ~0). Residue-complete so
no split-residue artifacts; AMBER charges; boundary (<1.8 A of QM) dropped.
QM = CR2 + Tyr203 phenol (matches qm_cthrp.xyz). Run in TeraChem env."""
import sys, numpy as np
from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, PDBFile, NoCutoff
from openmm import unit
sys.path.insert(0, "/home/robson/PetaChem")
import qmmm_tddft_pipeline as P
from scipy.spatial import cKDTree

R = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
prm = AmberPrmtopFile("/home/robson/PetaChem/anionic_build/monomer_solv.prmtop")
topo = prm.topology
pos = np.array(PDBFile("/home/robson/PetaChem/anionic_build/monomer_min.pdb")
               .getPositions().value_in_unit(unit.angstrom))
charges, _ = P.get_atomic_charges_from_system(prm.createSystem(nonbondedMethod=NoCutoff))

PHENOL = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HD1", "HD2", "HE1", "HE2", "HH"}
qm = set()
for r in topo.residues():
    if r.name == "CR2": qm |= {a.index for a in r.atoms()}
    elif r.name == "TYR" and r.id == "202": qm |= {a.index for a in r.atoms() if a.name in PHENOL}
qm_xyz = pos[sorted(qm)]; tree = cKDTree(qm_xyz)
qm_charge = int(round(sum(charges[i] for i in qm)))

mm = []
for res in topo.residues():
    ridx = [a.index for a in res.atoms() if a.index not in qm]
    if not ridx: continue
    if tree.query(pos[ridx])[0].min() < R:
        mm += ridx
mm = sorted(set(mm))
mx = pos[mm]; mq = np.array([charges[i] for i in mm])
keep = tree.query(mx)[0] > 1.8                       # drop link-boundary charges
mq, mx = mq[keep], mx[keep]
print(f"R={R} A: {len(mq)} MM charges (vs ~30k full); net={mq.sum():+.3f}  "
      f"QM={qm_charge:+d}  QM+MM={mq.sum()+qm_charge:+.3f}")
with open("/home/robson/PetaChem/neo_model/mm_light_cthrp.dat", "w") as f:
    f.write(f"{len(mq)}\nlightweight MM (residues within {R} A) for CR2+Tyr-phenol QM\n")
    for q, (x, y, z) in zip(mq, mx):
        f.write(f"{q: .6f} {x: .6f} {y: .6f} {z: .6f}\n")
print("wrote mm_light_cthrp.dat")

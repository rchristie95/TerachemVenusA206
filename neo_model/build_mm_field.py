#!/usr/bin/env python3
"""Build the protein MM point-charge field for the anionic-chromophore QM region
(cr2_chromophore.xyz), for a pyscf ADC(2) QM/MM single point. QM = the chromophore
(charge -1); MM = every other protein/water atom (incl. His148) at its AMBER charge,
with charges within 1.8 A of any QM atom dropped (link-atom boundary). Run in the
TeraChem env (openmm). Writes neo_model/mm_field.dat (q x y z, Angstrom)."""
import sys, numpy as np
from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, PDBFile, NoCutoff
from openmm import unit
sys.path.insert(0, "/home/robson/PetaChem")
import qmmm_tddft_pipeline as P

prm = AmberPrmtopFile("/home/robson/PetaChem/anionic_build/monomer_solv.prmtop")
topo = prm.topology
pos = np.array(PDBFile("/home/robson/PetaChem/anionic_build/monomer_min.pdb")
               .getPositions().value_in_unit(unit.angstrom))
charges, _ = P.get_atomic_charges_from_system(prm.createSystem(nonbondedMethod=NoCutoff))

L = open("/home/robson/PetaChem/neo_model/cr2_chromophore.xyz").readlines(); nq = int(L[0].split()[0])
qm_xyz = np.array([[float(c) for c in L[2+i].split()[1:4]] for i in range(nq)])

cr2_idx = set(a.index for r in topo.residues() if r.name == "CR2" for a in r.atoms())
print(f"system atoms={topo.getNumAtoms()}  net charge={charges.sum():+.3f}")
print(f"CR2 residue atoms={len(cr2_idx)}  CR2 charge={sum(charges[i] for i in cr2_idx):+.3f}")

mm_idx = [a.index for a in topo.atoms() if a.index not in cr2_idx]
mm_xyz = pos[mm_idx]; mm_q = np.array([charges[i] for i in mm_idx])

# drop MM charges within 1.8 A of any QM atom (link-atom boundary); preserve total charge
from scipy.spatial import cKDTree
dmin, _ = cKDTree(qm_xyz).query(mm_xyz)
keep = dmin > 1.8
q_keep = mm_q[keep].copy()
q_keep += (mm_q.sum() - q_keep.sum()) / keep.sum()      # spread dropped charge -> exact total
xyz_keep = mm_xyz[keep]
print(f"MM candidates={len(mm_idx)}  dropped(<1.8A)={int((~keep).sum())}  kept={int(keep.sum())}")
print(f"MM total charge={q_keep.sum():+.3f}  =>  QM(-1)+MM={q_keep.sum()-1:+.3f} (vs system {charges.sum():+.3f})")
print(f"nearest MM->QM after exclusion: {cKDTree(qm_xyz).query(xyz_keep)[0].min():.2f} A")

with open("/home/robson/PetaChem/neo_model/mm_field.dat", "w") as f:
    f.write(f"{len(q_keep)}\nMM point charges (protein incl His148) for cr2_chromophore QM\n")
    for q, (x, y, z) in zip(q_keep, xyz_keep):
        f.write(f"{q: .6f} {x: .6f} {y: .6f} {z: .6f}\n")
print("wrote neo_model/mm_field.dat")

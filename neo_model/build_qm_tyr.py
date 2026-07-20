#!/usr/bin/env python3
"""Deterministic QM/MM builder for the cyanine wavelength loop. Geometry = the
MM-MINIMIZED monomer (monomer_min.pdb; deterministic, no MD/NVT). QM = a fixed
residue list; MM = every other atom at its AMBER charge, charges <1.8 A of QM
dropped (link boundary), total charge preserved. Run in the TeraChem env (openmm).
Usage: python build_qm_tyr.py <tag> <resname:resid> [<resname:resid> ...]
Writes neo_model/qm_<tag>.xyz and neo_model/mm_<tag>.dat ."""
import sys, numpy as np
from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, PDBFile, NoCutoff
from openmm import unit
sys.path.insert(0, "/home/robson/PetaChem")
import qmmm_tddft_pipeline as P

tag = sys.argv[1]
want = set(sys.argv[2:])                      # e.g. {"CR2:66","TYR:202","HIS:147"}
prm = AmberPrmtopFile("/home/robson/PetaChem/anionic_build/monomer_solv.prmtop")
topo = prm.topology
pos = np.array(PDBFile("/home/robson/PetaChem/anionic_build/monomer_min.pdb")
               .getPositions().value_in_unit(unit.angstrom))
charges, _ = P.get_atomic_charges_from_system(prm.createSystem(nonbondedMethod=NoCutoff))

def key(r): return f"{r.name}:{r.id}"
qm_res = [r for r in topo.residues() if key(r) in want]
print("QM residues:", [key(r) for r in qm_res], " (requested", want, ")")
assert len(qm_res) == len(want), "some requested residues not found"

qm_records = P.build_qm_atom_records(topo, pos, set(qm_res))
links, _ = P.build_link_atom_records(topo, pos, {r["global_index"] for r in qm_records})
qm = qm_records + links
qm_global = {r["global_index"] for r in qm_records}
qm_charge = int(round(sum(charges[i] for i in sorted(qm_global))))
print(f"QM real={len(qm_records)} link-H={len(links)} total={len(qm)} charge={qm_charge}")
P.write_xyz(qm, f"/home/robson/PetaChem/neo_model/qm_{tag}.xyz",
            f"QM={sorted(want)} (MM-min geom) charge={qm_charge}")

mm_idx = [a.index for a in topo.atoms() if a.index not in qm_global]
mm_xyz = pos[mm_idx]; mm_q = np.array([charges[i] for i in mm_idx])
from scipy.spatial import cKDTree
d, _ = cKDTree(np.array([r["coord"] for r in qm])).query(mm_xyz)
keep = d > 1.8
qk = mm_q[keep].copy(); qk += (mm_q.sum() - qk.sum()) / keep.sum(); xk = mm_xyz[keep]
print(f"MM kept={int(keep.sum())} dropped(<1.8A)={int((~keep).sum())} "
      f"MMcharge={qk.sum():+.3f}  QM+MM={qk.sum()+qm_charge:+.3f}")
with open(f"/home/robson/PetaChem/neo_model/mm_{tag}.dat", "w") as f:
    f.write(f"{len(qk)}\nMM field for QM={sorted(want)}\n")
    for q, (x, y, z) in zip(qk, xk):
        f.write(f"{q: .6f} {x: .6f} {y: .6f} {z: .6f}\n")
# report chromophore-O to Tyr203-ring stacking distance if Tyr present
print(f"wrote qm_{tag}.xyz ({len(qm)} atoms) + mm_{tag}.dat ({len(qk)} charges)")

#!/usr/bin/env python3
"""Build QM = full chromophore (CR2) + Tyr203 PHENOL RING only (the pi-stacker),
trimmed to keep ADC(2) feasible. Tyr203 backbone+CB dropped; the CG-CB bond is
capped with a link-H. Geometry = MM-minimized (deterministic). Gas-phase (no MM,
since its excitation effect is ~0). Run in the TeraChem env. Writes qm_cthrp.xyz."""
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

PHENOL = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HD1", "HD2", "HE1", "HE2", "HH"}
qm_global = set()
for r in topo.residues():
    if r.name == "CR2":
        qm_global |= {a.index for a in r.atoms()}
    elif r.name == "TYR" and r.id == "202":
        qm_global |= {a.index for a in r.atoms() if a.name in PHENOL}
print(f"QM real atoms: {len(qm_global)} (CR2 + Tyr203 phenol ring)")

qm_records = P.build_qm_atom_records(topo, pos, None, atom_indices=qm_global) \
    if False else None
# build records manually from indices (preserve order)
atoms = list(topo.atoms())
recs = []
for i in sorted(qm_global):
    a = atoms[i]
    el = a.element.symbol if a.element else "C"
    recs.append({"global_index": i, "symbol": el, "coord": pos[i]})
links, _ = P.build_link_atom_records(topo, pos, qm_global)
qm = recs + links
print(f"QM real={len(recs)} link-H={len(links)} total={len(qm)}")
P.write_xyz(qm, "/home/robson/PetaChem/neo_model/qm_cthrp.xyz",
            "QM = CR2 + Tyr203 phenol (MM-min geom), charge -1")
syms = [r["symbol"] for r in qm]
print("formula:", {s: syms.count(s) for s in sorted(set(syms))})
# stacking distance: min C(chromophore phenolate ring) - C(Tyr phenol ring)
print(f"wrote qm_cthrp.xyz ({len(qm)} atoms)")

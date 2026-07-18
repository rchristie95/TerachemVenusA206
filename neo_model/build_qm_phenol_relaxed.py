#!/usr/bin/env python3
"""QM = CR2 (QM/MM-relaxed, C-O 1.263 from the deterministic constrained opt) +
Tyr203 phenol (MM-min). Builds the full structure = monomer_min with CR2 atoms
swapped to their constrained-opt positions, then extracts QM + link atoms.
Run in TeraChem env. Writes qm_cthrp_relaxed.xyz."""
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
box = P.get_periodic_box_lengths_ang(topo); pbc = box is not None
atoms = list(topo.atoms())

# rebuild the 274-atom QM-region order (== qm_deprotonated/qm_opt order) to map -> global index
frozen = P.load_qm_selection("/home/robson/PetaChem/tc_simple_anionic/qm_selection.json")
selected, *_ = P.select_qm_residues(topo, pos, "CR2", 2.65, nearest_waters=5,
                                    include_resids=set(), frozen_keys=frozen,
                                    box_lengths_a=box, use_periodic=pbc)
recs = P.build_qm_atom_records(topo, pos, selected)
ql = open("/home/robson/PetaChem/tc_qmmm_opt_constrained/qm_opt.xyz").readlines()
qxyz = np.array([[float(c) for c in ql[2+i].split()[1:4]] for i in range(274)])

full = pos.copy(); ncr2 = 0
for i, r in enumerate(recs):
    if atoms[r["global_index"]].residue.name == "CR2":
        full[r["global_index"]] = qxyz[i]; ncr2 += 1
print(f"swapped {ncr2} CR2 atoms to QM/MM-relaxed positions")

PHENOL = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HD1", "HD2", "HE1", "HE2", "HH"}
qm = set()
for r in topo.residues():
    if r.name == "CR2": qm |= {a.index for a in r.atoms()}
    elif r.name == "TYR" and r.id == "202": qm |= {a.index for a in r.atoms() if a.name in PHENOL}
qmrec = [{"global_index": i, "symbol": atoms[i].element.symbol if atoms[i].element else "C",
          "coord": full[i]} for i in sorted(qm)]
links, _ = P.build_link_atom_records(topo, full, qm)
qmrec += links
P.write_xyz(qmrec, "/home/robson/PetaChem/neo_model/qm_cthrp_relaxed.xyz",
            "CR2(relaxed C-O 1.263) + Tyr203 phenol, charge -1")
# verify the chromophore phenolate C-O
x = np.array([r["coord"] for r in qmrec]); s = [r["symbol"] for r in qmrec]
for i in range(29):
    if s[i] == 'O':
        d = np.linalg.norm(x - x[i], axis=1)
        nb = [(j, s[j]) for j in range(len(x)) if j != i and d[j] < 1.45]
        if len(nb) == 1 and nb[0][1] == 'C':
            print(f"  phenolate-type O{i}-C{nb[0][0]} = {d[nb[0][0]]:.3f} A")
print(f"wrote qm_cthrp_relaxed.xyz ({len(qmrec)} atoms)")

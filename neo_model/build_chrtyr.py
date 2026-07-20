#!/usr/bin/env python3
"""Build a CLEAN in-protein QM region with NO His148 (His148 -> MM field), to remove
the His148->chromophore CT intruder that broke STEOM-CCSD on the 73-atom model.

QM = CR2 (QM/MM-relaxed, C-O 1.263, same swap as build_qm_phenol_relaxed.py)
   + FULL Tyr203 residue (all atoms of TYR 202), H-link capped.
MM = every other residue within R A of QM (INCLUDING His148), AMBER charges,
     charges within 1.8 A of any QM atom dropped (link boundary), net preserved.

Writes:
  orca_steom/geom_chrtyr.xyz   (QM geometry + link H)
  orca_steom/field_chrtyr.pc   (ORCA %pointcharges: line1=N, then 'q x y z' Angstrom)
Run in the TeraChem env (openmm).
"""
import sys, numpy as np
from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, PDBFile, NoCutoff
from openmm import unit
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/robson/PetaChem")
import qmmm_tddft_pipeline as P

OUT = "/home/robson/PetaChem/neo_model/orca_steom"
prm = AmberPrmtopFile("/home/robson/PetaChem/anionic_build/monomer_solv.prmtop")
topo = prm.topology
pos = np.array(PDBFile("/home/robson/PetaChem/anionic_build/monomer_min.pdb")
               .getPositions().value_in_unit(unit.angstrom))
charges, _ = P.get_atomic_charges_from_system(prm.createSystem(nonbondedMethod=NoCutoff))
atoms = list(topo.atoms())

# --- swap CR2 atoms to the QM/MM-relaxed positions (identical to build_qm_phenol_relaxed.py) ---
frozen = P.load_qm_selection("/home/robson/PetaChem/tc_simple_anionic/qm_selection.json")
box = P.get_periodic_box_lengths_ang(topo); pbc = box is not None
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

# --- QM region: CR2 (all) + FULL Tyr203 (all atoms of TYR 202).  His148 deliberately excluded. ---
qm = set()
for r in topo.residues():
    if r.name == "CR2":
        qm |= {a.index for a in r.atoms()}
    elif r.name == "TYR" and r.id == "202":
        qm |= {a.index for a in r.atoms()}
qm_charge = int(round(sum(charges[i] for i in qm)))
qmrec = [{"global_index": i, "symbol": atoms[i].element.symbol if atoms[i].element else "C",
          "coord": full[i]} for i in sorted(qm)]
links, _ = P.build_link_atom_records(topo, full, qm)
qmrec += links
P.write_xyz(qmrec, f"{OUT}/geom_chrtyr.xyz",
            f"CR2(relaxed C-O 1.263) + FULL Tyr203, His148 in MM, charge {qm_charge}")
print(f"QM: {len(qm)} heavy/real + {len(links)} link H = {len(qmrec)} atoms ; QM charge {qm_charge:+d}")

# --- MM field: whole residues within R A of QM, excl QM, drop <1.8 A; preserve net charge ---
R = 12.0
qm_xyz = full[sorted(qm)]; tree = cKDTree(qm_xyz)
mm = []
for res in topo.residues():
    ridx = [a.index for a in res.atoms() if a.index not in qm]
    if ridx and tree.query(full[ridx])[0].min() < R:
        mm += ridx
mm = sorted(set(mm))
mx = full[mm]; mq = np.array([charges[i] for i in mm])
keep = tree.query(mx)[0] > 1.8
mq_k = mq[keep].copy(); mx_k = mx[keep]
mq_k += (mq.sum() - mq_k.sum()) / keep.sum()             # spread dropped charge -> exact MM total

# verify His148 (topo id "147", HID/HIE/HIP) made it into the field
his_idx = [a.index for r in topo.residues() if r.name in ("HID","HIE","HIP","HIS") and r.id == "147"
           for a in r.atoms()]
his_in_mm = [i for i in his_idx if i in set(mm)]
his_dropped = [i for i in his_idx if i in set(mm) and tree.query(full[[i]])[0][0] <= 1.8]
print(f"MM: {len(mq_k)} charges within {R} A (dropped <1.8A: {int((~keep).sum())}); "
      f"MM net {mq_k.sum():+.3f} ; QM+MM {mq_k.sum()+qm_charge:+.3f}")
print(f"His148: residue atoms={len(his_idx)}, in MM field={len(his_in_mm)}, "
      f"dropped(<1.8A of QM)={len(his_dropped)}, "
      f"nearest His->QM={tree.query(full[his_idx])[0].min():.2f} A" if his_idx else "His148 NOT FOUND")

with open(f"{OUT}/field_chrtyr.pc", "w") as f:
    f.write(f"{len(mq_k)}\n")
    for q, (x, y, z) in zip(mq_k, mx_k):
        f.write(f"{q: .6f} {x: .6f} {y: .6f} {z: .6f}\n")
print(f"wrote {OUT}/geom_chrtyr.xyz and field_chrtyr.pc")

# --- sanity: phenolate C-O bond length in the QM geometry ---
sym = [r["symbol"] for r in qmrec]; xyz = np.array([r["coord"] for r in qmrec])
for i, s in enumerate(sym):
    if s == "O":
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        nb = [(j, d[j]) for j in range(len(sym)) if sym[j] == "C" and 0.1 < d[j] < 1.6]
        nb.sort(key=lambda t: t[1])
        if nb and nb[0][1] < 1.45:
            print(f"  phenolate-type O(idx{i})-C(idx{nb[0][0]}) = {nb[0][1]:.3f} A")

#!/usr/bin/env python3
"""TDDFT-screen candidate QM regions / QM-MM link schemes for STEOM stability.

Rationale: STEOM-CCSD crashes when near-degenerate IP/EA satellite states coalesce into
complex eigenpairs. In TDDFT these show up as DARK charge-transfer states sitting below the
bright pi->pi*. TDDFT is cheap + has no EOM solver, so we use it to count dark-below-bright
states for each candidate region; fewest = most STEOM-stable. CR2 is always QM (relaxed
C-O 1.263); residues are added at chosen 'extent' (phenol ring / sidechain / full residue)
= the link scheme (where the bond is cut, then H-capped). Run in TeraChem env.

For each spec, writes tc_screen_<name>/{geometry.xyz, mm_charges.dat, tddft.in}.
"""
import sys, os, numpy as np
from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, PDBFile, NoCutoff
from openmm import unit
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/robson/PetaChem")
import qmmm_tddft_pipeline as P

BASE = "/home/robson/PetaChem"
prm = AmberPrmtopFile(f"{BASE}/anionic_build/monomer_solv.prmtop")
topo = prm.topology
pos = np.array(PDBFile(f"{BASE}/anionic_build/monomer_min.pdb").getPositions().value_in_unit(unit.angstrom))
charges, _ = P.get_atomic_charges_from_system(prm.createSystem(nonbondedMethod=NoCutoff))
atoms = list(topo.atoms())

# CR2 -> relaxed positions (C-O 1.263)
frozen = P.load_qm_selection(f"{BASE}/tc_simple_anionic/qm_selection.json")
box = P.get_periodic_box_lengths_ang(topo); pbc = box is not None
selected, *_ = P.select_qm_residues(topo, pos, "CR2", 2.65, nearest_waters=5,
                                    include_resids=set(), frozen_keys=frozen,
                                    box_lengths_a=box, use_periodic=pbc)
recs = P.build_qm_atom_records(topo, pos, selected)
ql = open(f"{BASE}/tc_qmmm_opt_constrained/qm_opt.xyz").readlines()
qxyz = np.array([[float(c) for c in ql[2+i].split()[1:4]] for i in range(274)])
full = pos.copy()
for i, r in enumerate(recs):
    if atoms[r["global_index"]].residue.name == "CR2":
        full[r["global_index"]] = qxyz[i]

BACKBONE = {"N","H","H1","H2","H3","CA","HA","HA2","HA3","C","O","OXT"}
PHENOL   = {"CG","CD1","CD2","CE1","CE2","CZ","OH","HD1","HD2","HE1","HE2","HH"}

def res_atoms(pred, resid, extent):
    out = set()
    for r in topo.residues():
        if pred(r.name) and r.id == resid:
            for a in r.atoms():
                if extent == "full": out.add(a.index)
                elif extent == "sidechain" and a.name not in BACKBONE: out.add(a.index)
                elif extent == "phenol" and a.name in PHENOL: out.add(a.index)
    return out

cr2 = set(a.index for r in topo.residues() if r.name == "CR2" for a in r.atoms())
isTYR = lambda n: n == "TYR"
isHIS = lambda n: n.startswith("HI")
isSER = lambda n: n == "SER"

# Tyr203 = topo id "202", His148 = topo id "147", Ser205 = topo id "204" (AMBER off-by-one)
specs = {
  "cr2only":   [],
  "phenol":    [(isTYR,"202","phenol")],                                  # = the 44-atom (STEOM OK)
  "tyrsc":     [(isTYR,"202","sidechain")],                               # +CB vs phenol ring
  "tyrfull":   [(isTYR,"202","full")],                                    # = the 54-atom (EA intruder)
  "hisonly":   [(isHIS,"147","sidechain")],                              # isolate His148 CT
  "phenolhis": [(isTYR,"202","phenol"), (isHIS,"147","sidechain")],       # phenol + His sidechain
  "phenolser": [(isTYR,"202","phenol"), (isSER,"204","sidechain")],       # phenol + Ser205 (H-bond)
}

print(f"{'spec':<11} {'QM atoms':>9} {'charge':>7} {'MM chg':>7}  His148")
for name, adds in specs.items():
    qm = set(cr2)
    for pred, rid, ext in adds:
        qm |= res_atoms(pred, rid, ext)
    qmc = int(round(sum(charges[i] for i in qm)))
    qmrec = [{"global_index": i, "symbol": atoms[i].element.symbol if atoms[i].element else "C",
              "coord": full[i]} for i in sorted(qm)]
    links, _ = P.build_link_atom_records(topo, full, qm)
    qmrec += links
    d = f"{BASE}/tc_screen_{name}"; os.makedirs(d, exist_ok=True)
    P.write_xyz(qmrec, f"{d}/geometry.xyz", f"screen {name}: charge {qmc}")
    qm_xyz = full[sorted(qm)]; tree = cKDTree(qm_xyz)
    mm = []
    for res in topo.residues():
        ridx = [a.index for a in res.atoms() if a.index not in qm]
        if ridx and tree.query(full[ridx])[0].min() < 12.0: mm += ridx
    mm = sorted(set(mm)); mx = full[mm]; mq = np.array([charges[i] for i in mm])
    keep = tree.query(mx)[0] > 1.8; mqk = mq[keep].copy(); mxk = mx[keep]
    mqk += (mq.sum() - mqk.sum()) / keep.sum()
    his_qm = any(atoms[i].residue.name.startswith("HI") for i in qm)
    with open(f"{d}/mm_charges.dat", "w") as f:
        f.write(f"{len(mqk)}\nMM Point Charges\n")
        for q, (x, y, z) in zip(mqk, mxk): f.write(f"{q:.8f} {x:.8f} {y:.8f} {z:.8f}\n")
    with open(f"{d}/tddft.in", "w") as f:
        f.write(f"coordinates geometry.xyz\nrun energy\ncis yes\nbasis 6-311g**\nmethod wb97xd3\n"
                f"charge {qmc}\nspinmult 1\npointcharges mm_charges.dat\nscrdir scr\n"
                f"cisnumstates 6\ncismaxiter 200\ncismax 500\nscf diis+a\nthreall 1.0e-13\ngpus 1\nend\n")
    print(f"{name:<11} {len(qmrec):>9} {qmc:>7} {len(mqk):>7}  {'QM' if his_qm else 'MM/absent'}")
print("done")

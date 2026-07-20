#!/usr/bin/env python3
"""
Build the QM/MM TDDFT setup (qm_deprotonated.xyz, mm_charges.dat,
qm_setup_settings.in) from the AMBER (anionic published-FF) relaxed monomer,
reusing terachem_full_pipeline's Stage-1 helpers but:
  - taking topology/charges from the AMBER prmtop (published CR2 RESP charges),
  - SKIPPING CR2 deprotonation (the published CR2 is already the anionic phenolate),
  - using the same QM selection / embedding parameters as the original Stage 1.

Output dir (default tc_simple_anionic) is consumed by Stage 2 via
  terachem_full_pipeline.py --skip-simple --tddft-args "--input-dir tc_simple_anionic ..."
"""
import sys
from pathlib import Path
import numpy as np

from openmm import Platform
Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile, NoCutoff
from openmm import unit

import terachem_full_pipeline as P

PRMTOP = "anionic_build/monomer_solv.prmtop"
INPCRD = "anionic_build/monomer_solv.inpcrd"
MINPDB = "anionic_build/monomer_min.pdb"
OUT = Path("tc_simple_anionic"); OUT.mkdir(exist_ok=True)

prm = AmberPrmtopFile(PRMTOP)
inp = AmberInpcrdFile(INPCRD)
topology = prm.topology
if inp.boxVectors is not None:
    topology.setPeriodicBoxVectors(inp.boxVectors)

# minimized positions (atom order matches prmtop topology)
positions_ang = np.array(PDBFile(MINPDB).getPositions().value_in_unit(unit.angstrom))

# charges from the AMBER system (published CR2 RESP + ff14SB)
system = prm.createSystem(nonbondedMethod=NoCutoff)
atom_charges, charge_source = P.get_atomic_charges_from_system(system)

box_lengths_a = P.get_periodic_box_lengths_ang(topology)
use_pbc = box_lengths_a is not None

# 1) QM residue selection: CR2 + protein within 2.65 A + 5 nearest waters.
#    AUDIT FIX: the YFP-defining Tyr203 pi-stacks with the chromophore and *is*
#    the YFP red-shift, but in the anionic-relaxed geometry it sits at 2.92 A
#    (just outside the 2.65 A cutoff) and was dropped, blue-shifting the
#    spectrum to 421 nm. Force-include the pi-stacking Tyr (closest TYR to CR2
#    beyond the cutoff, within 4 A) so the YFP chromophore is described correctly.
cr2_idx = [a.index for r in topology.residues() if r.name == "CR2" for a in r.atoms()]
cr2_xyz = positions_ang[cr2_idx]
force_ids = set()
best = None
for res in topology.residues():
    if res.name != "TYR":
        continue
    ridx = [a.index for a in res.atoms()]
    dmin = float(np.min(np.linalg.norm(positions_ang[ridx][:, None, :] - cr2_xyz[None, :, :], axis=2)))
    if dmin < 4.0 and (best is None or dmin < best[0]):
        # the closest TYR overall is the pi-stacker (Tyr203)
        pass
    if 2.65 <= dmin < 4.0:
        # candidate just outside cutoff
        if best is None or dmin < best[0]:
            best = (dmin, res.id)
if best is not None:
    force_ids.add(best[1])
    print(f"[*] Force-including pi-stacking Tyr203 (res {best[1]}, {best[0]:.2f} A from CR2)")

# Reproducible QM region: reuse a frozen selection if present, else compute it
# (distance + forced Tyr203) and persist it. This guarantees the identical QM
# residue set across every NVT frame despite the stochastic dynamics.
SELFILE = OUT / "qm_selection.json"
frozen = P.load_qm_selection(SELFILE) if SELFILE.exists() else None
if frozen:
    print(f"[*] Using frozen QM selection from {SELFILE} ({len(frozen)} residues)")
selected, cr2_res, n_w, avail_w, far_w = P.select_qm_residues(
    topology, positions_ang, "CR2", 2.65, nearest_waters=5,
    include_resids=force_ids,
    box_lengths_a=box_lengths_a, use_periodic=use_pbc,
    frozen_keys=frozen)
if not frozen:
    P.save_qm_selection(selected, SELFILE)
    print(f"[*] Saved QM selection to {SELFILE} for reproducible reuse")
print(f"[*] QM residues: {len(selected)} (CR2 + neighbours + {n_w} waters)")
print("    " + ", ".join(sorted(f"{r.name}{r.id}" for r in selected)))

# 2) QM atom records — NO deprotonation (CR2 already anionic)
qm_records = P.build_qm_atom_records(topology, positions_ang, selected)
qm_global = {r["global_index"] for r in qm_records}

# 3) link atoms on severed bonds
links, cut_bonds = P.build_link_atom_records(topology, positions_ang, qm_global)
qm_with_links = qm_records + links
print(f"[*] QM real atoms: {len(qm_records)}  link-H: {len(links)}")

# 4) QM charge = rounded sum of MM charges over QM atoms (anionic CR2 included)
qm_sum = float(np.sum(atom_charges[sorted(qm_global)]))
tc_charge = int(np.rint(qm_sum))
print(f"[*] QM charge: {tc_charge} (raw MM sum {qm_sum:+.4f})")

# 5) MM embedding (PCM active -> exclude non-QM waters), matching Stage-1 defaults
WATER = P.WATER_RESIDUE_NAMES
nonqm_water = {a.index for a in topology.atoms()
               if a.residue.name in WATER and a.index not in qm_global}
mm_idx, _, uncapped, capd = P.select_mm_embedding_indices(
    topology, positions_ang, qm_global, set(nonqm_water),
    cutoff_a=0.0, min_distance_a=1.2, exclusion_hops=2, max_point_charges=0,
    box_lengths_a=box_lengths_a, use_periodic=use_pbc)
qm_ref = np.array([r["coord"] for r in qm_with_links], dtype=float)
act_idx, act_q, stats = P.apply_mm_short_range_repulsion(
    list(mm_idx), positions_ang, atom_charges, qm_ref, 1.8, 2.8,
    box_lengths_a=box_lengths_a, use_periodic=use_pbc, topology=topology,
    preserve_residue_charge=True, preserve_total_charge=True)
print(f"[*] MM point charges: {len(act_idx)} (net {float(np.sum(act_q)):+.3f} e)")

# 6) write Stage-2 inputs
P.write_xyz(qm_with_links, OUT / "qm_deprotonated.xyz", "QM region (anionic CR2) for TD-DFT")
P.write_mm_pointcharges(OUT / "mm_charges.dat", positions_ang, atom_charges, act_idx, charge_values=act_q)
P.write_qm_setup_settings(OUT / "qm_setup_settings.in", "wb97xd3", "6-311g**",
                          tc_charge, 1, "cosmo", 78.39, "iswig", 1.40)

# 7) protein-only relaxed monomer (chain A) for Stage-3 PyMOL alignment
prot = []
for atom, pos in zip(topology.atoms(), positions_ang):
    rn = atom.residue.name
    if rn in WATER or rn in {"NA","CL","Na+","Cl-"}: continue
    el = atom.element.symbol if atom.element else "C"
    name = atom.name
    nm = name if len(name) == 4 else " " + name
    prot.append(f"ATOM  {len(prot)+1:5d} {nm:<4s} {rn[:3]:>3s} A{atom.residue.id:>4}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {el:>2s}\n")
(OUT / "monomer_relaxed.pdb").write_text("".join(prot) + "TER\nEND\n")
print(f"[*] wrote {OUT}/ : qm_deprotonated.xyz, mm_charges.dat, qm_setup_settings.in, monomer_relaxed.pdb")
print(f"    protein atoms in monomer_relaxed.pdb: {len(prot)}")

#!/usr/bin/env python3
"""Multi-model protein-only PDB for render_video(), unwrapped across PBC.

render_video() in run_nvt.py reads a multi-model PDB (count_pdb_models), not a
DCD, which is why ORCA_ONLY_PIPELINE.md passes --openmm-trajectory-file *.pdb.

Topology is solvated_protonated.pdb: OpenMM built the System from it, so its
atom order IS the DCD order by construction. protein_topology.pdb is a
post-imaged subset and does not share that order.

PBC: the DCD stores raw wrapped coordinates, so residues are unwrapped
SEQUENTIALLY along the chain (each residue minimum-imaged against the previous
one). A single global anchor would fail here -- the protein spans ~60 A in a
91.86 A box, i.e. more than the half-box a one-anchor minimum image assumes.
"""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, 'tddft_ensemble_decoherence/scripts')
from prepare_production_frame import read_dcd_frame

DCD  = Path('tc_tandem_v4_correctseq/tandem_v4_1ns.dcd')
TOP  = 'tc_tandem_v4_correctseq/solvated_protonated.pdb'
OUT  = Path('tandem_v4_correct_seq/tandem_v4_traj_300.pdb')
NFRAMES, TOTAL = 300, 1200
SOLVENT = {'HOH','WAT','TIP3','NA','CL','K','SOD','CLA','NA+','CL-'}

recs = [l.rstrip('\n') for l in open(TOP) if l.startswith(('ATOM','HETATM'))]
keep, resid_of = [], []
for i, l in enumerate(recs):
    if l[17:20].strip() not in SOLVENT:
        keep.append(i); resid_of.append(int(l[22:26]))
keep = np.array(keep); resid_of = np.array(resid_of)
uniq = sorted(set(resid_of.tolist()))
groups = [np.where(resid_of == r)[0] for r in uniq]
print(f'protein atoms kept {len(keep)} in {len(groups)} residues')

def unwrap(P, box):
    """Sequential minimum-image unwrap, residue by residue along the chain."""
    out = P.copy()
    prev = out[groups[0]].mean(0)
    for g in groups:
        c = out[g].mean(0)
        shift = np.round((c - prev) / box) * box
        out[g] -= shift
        prev = out[g].mean(0)
    return out

idx = np.linspace(0, TOTAL - 1, NFRAMES).astype(int)
with OUT.open('w') as fh:
    for m, f in enumerate(idx, 1):
        xyz, box, _ = read_dcd_frame(DCD, int(f))
        P = unwrap(xyz[keep], np.asarray(box, dtype=float))
        P -= P.mean(0)                      # centre for a stable camera
        fh.write(f'MODEL     {m:4d}\n')
        for (li, (x, y, z)) in zip(keep, P):
            fh.write(f'{recs[li][:30]}{x:8.3f}{y:8.3f}{z:8.3f}{recs[li][54:]}\n')
        fh.write('ENDMDL\n')
fh_extent = None
print(f'wrote {OUT} : {NFRAMES} models from frames {idx[0]}..{idx[-1]}')

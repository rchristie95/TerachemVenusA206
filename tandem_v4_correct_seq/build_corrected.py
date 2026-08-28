#!/usr/bin/env python3
"""Restore the two LGMDELYK barrel C-terminal tails absent from tandem_dimer_v2.pdb.

v2 already carries the correct linker Gly (v1 had a spurious Cys). The only
remaining deviations from Addgene #250938 (pdVenus-TD, 511 aa) are the two
8-residue C-terminal tails LGMDELYK, which are disordered in 1MYW and therefore
absent from the crystal barrels the tandem was built from.

Alignment of the built 491-residue chain against the deposit gave exactly:
    insert LGMDELYK after built residue 229  (barrel 1 C-term, MID-CHAIN)
    C246G                                    (already correct in v2)
    insert LGMDELYK after built residue 491  (barrel 2 C-term, TERMINAL)
plus two X<->GYG rows that are the mature CR2 chromophores and are expected.

491 + 16 = 507 residues, which is the deposit's 511 with each CR2 standing in
for a GYG tripeptide.
"""
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import sys

TAIL = ['LEU','GLY','MET','ASP','GLU','LEU','TYR','LYS']

fixer = PDBFixer(filename='/home/robson/PetaChem/tandem_dimer_v2.pdb')
chains = list(fixer.topology.chains())
res = list(chains[0].residues())
print(f'input: {len(chains)} chain(s), {len(res)} residues')
print(f'  residue at index 229 = {res[229].name} (linker start, tail goes BEFORE it)')
print(f'  final residue        = {res[-1].name}')

fixer.findMissingResidues()
print(f'  PDBFixer auto-detected missing residues: {fixer.missingResidues}')

# insert BEFORE index 229 (barrel-1 C-term) and at the very end (barrel-2 C-term)
fixer.missingResidues = {(0, 229): TAIL, (0, len(res)): TAIL}
fixer.findNonstandardResidues()
print(f'  nonstandard (left untouched): {[(r.name) for r,_ in fixer.nonstandardResidues]}')
fixer.findMissingAtoms()
fixer.addMissingAtoms()

out = 'tandem_dimer_v4_511aa.pdb'
with open(out,'w') as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=False)
n = len(list(list(fixer.topology.chains())[0].residues()))
print(f'wrote {out}: {n} residues (expect 507)')

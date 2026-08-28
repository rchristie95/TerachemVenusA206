#!/usr/bin/env python3
"""Replace the two PDBFixer-built LGMDELYK tails with AlphaFold3 conformations.

Only the 8-residue tails move. Each tail is superposed using the 15 barrel
residues immediately preceding it, so the A206 crystal register that the
coupling depends on is untouched -- AF3's inter-domain register is known to be
wrong (11.7 A from crystal) and is deliberately NOT imported.

Residue mapping (v4 has CR2 as one residue where AF3 has the GYG precursor):
    v4   1-65    <-> AF3   1-65
    v4  66 (CR2) <-> AF3  66-68 (GYG)
    v4  67-335   <-> AF3  69-337
    v4 336 (CR2) <-> AF3 338-340
    v4 337-507   <-> AF3 341-511
"""
import numpy as np, sys

CIF = sys.argv[1] if len(sys.argv) > 1 else \
    '/tmp/claude-1000/-home-robson-PetaChem/05591632-d63e-4446-a0e3-108a439ad793/scratchpad/af511/fold_dvenus_td511_s1_model_0.cif'
SRC, OUT = 'tandem_dimer_v4_511aa.pdb', 'tandem_dimer_v4_af3tails.pdb'

def v4_to_af3(r):
    if r <= 65:  return r
    if r == 66:  return None
    if r <= 335: return r + 2
    if r == 336: return None
    return r + 4

# ---- read AF3 ----
af = {}
with open(CIF) as fh:
    inloop = False
    for ln in fh:
        if ln.startswith('_atom_site.'): inloop = True; continue
        if inloop:
            p = ln.split()
            if len(p) < 18: 
                if af: break
                continue
            if p[0] != 'ATOM': continue
            af.setdefault(int(p[15]), {})[p[3]] = np.array([float(p[10]), float(p[11]), float(p[12])])
print(f'AF3 residues read: {len(af)}')

# ---- read v4 ----
lines = [l for l in open(SRC) if l.startswith(('ATOM','HETATM'))]
def rnum(l): return int(l[22:26])
def aname(l): return l[12:16].strip()

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P-pc).T @ (Q-qc)
    U,S,Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T@U.T))
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return R, qc - R@pc

BB = ['N','CA','C','O']
out = {}
for label, tail, anchor in [('barrel-1', range(230,238), range(215,230)),
                            ('barrel-2', range(500,508), range(485,500))]:
    P, Q = [], []
    for r in anchor:
        a = v4_to_af3(r)
        if a is None or a not in af: continue
        v4at = {aname(l): np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                for l in lines if rnum(l)==r}
        for nm in BB:
            if nm in v4at and nm in af[a]:
                P.append(af[a][nm]); Q.append(v4at[nm])
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    rmsd = np.sqrt((((P@R.T)+t - Q)**2).sum(1).mean())
    print(f'{label}: anchored on {len(P)} backbone atoms, superposition RMSD {rmsd:.2f} A')
    for r in tail:
        a = v4_to_af3(r)
        for nm, xyz in af[a].items():
            out[(r,nm)] = R@xyz + t

# ---- write ----
n=0
with open(OUT,'w') as fh:
    for l in lines:
        k = (rnum(l), aname(l))
        if k in out:
            x,y,z = out[k]; n+=1
            fh.write(f'{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}')
        else:
            fh.write(l)
    fh.write('END\n')
print(f'wrote {OUT}: {n} tail atoms replaced with AF3 conformations')

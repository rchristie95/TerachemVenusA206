#!/usr/bin/env python3
"""Huang-Rhys factors from the TeraChem Hessian + excited-state gradient.

TeraChem has no ESD-equivalent module, so the last step is done here. In the
vertical-gradient (linear coupling) approximation, projecting the excited-state
gradient at the Franck-Condon point onto the ground-state mass-weighted normal
modes gives, in atomic units,

    S_i = g_i^2 / (2 w_i^3)

with g_i the mass-weighted gradient along mode i. The reorganisation energy is
lambda = sum_i S_i w_i, which is the physical sanity check on the result.

ALL parameters here are computed: nothing is taken from experiment or assumed.
Basis is 6-31G* throughout (see the tc/*.in headers for why the chain must be
single-basis, and note this differs from the 6-311+G** of the site energies).
"""
import numpy as np, json
from pathlib import Path

HARTREE_CM = 219474.6313702
AMU_ME = 1822.888486209
F = Path('tc/scr_freq'); G = Path('tc/scr_esgrad')

# --- frequencies -------------------------------------------------------------
L = (F/'Frequencies.dat').read_text().split('\n')
freq = []
for i, l in enumerate(L):
    if '---' in l and i > 0:
        try: freq += [float(x) for x in L[i-1].split()]
        except ValueError: pass
freq = np.array(freq)

# --- mass-weighted normal modes ----------------------------------------------
modes, cur = [], []
for l in (F/'Mass.weighted.modes.dat').read_text().split('\n'):
    if l.startswith('=== Mode'):
        if cur: modes.append(cur)
        cur = []
    else:
        try: cur.append(float(l))
        except ValueError: pass
if cur: modes.append(cur)
Lmw = np.array(modes)                      # (nmode, 3N)
Lmw /= np.linalg.norm(Lmw, axis=1, keepdims=True)

# --- masses and excited-state gradient ---------------------------------------
mass = np.array([float(l.split()[1]) for l in (F/'Atom.mass.dat').read_text().split('\n')[1:] if l.strip()])
gl = (G/'grad.xyz').read_text().split('\n')
grad = np.array([[float(x) for x in l.split()[1:4]] for l in gl[2:2+len(mass)]])

# --- mass-weight the gradient, project, and form S ---------------------------
gmw = (grad / np.sqrt(mass * AMU_ME)[:, None]).ravel()
gi = Lmw @ gmw
w_au = freq / HARTREE_CM
ok = freq > 1.0
S = np.zeros_like(freq)
S[ok] = gi[ok]**2 / (2.0 * w_au[ok]**3)

lam = float((S * freq).sum())
print(f"{len(freq)} modes, {int((freq<0).sum())} imaginary")
print(f"total reorganisation energy  lambda = sum S_i w_i = {lam:.1f} cm^-1")
print(f"total Huang-Rhys             sum S_i             = {S.sum():.3f}\n")
o = np.argsort(S)[::-1][:12]
print(f"{'omega (cm^-1)':>14} {'S':>9} {'S*omega':>10}")
for i in o:
    print(f"{freq[i]:14.1f} {S[i]:9.4f} {S[i]*freq[i]:10.1f}")

# effective single mode for the one-mode Holstein model used in the CD code
hi = (freq > 1000) & (freq < 1800)
S_eff = float(S[hi].sum())
w_eff = float((S[hi]*freq[hi]).sum()/S[hi].sum()) if S[hi].sum() > 0 else float('nan')
print(f"\nHIGH-FREQUENCY WINDOW 1000-1800 cm^-1 (what the one-mode CD model represents)")
print(f"  S_eff = sum S_i             = {S_eff:.3f}")
print(f"  w_eff = <w>_S weighted      = {w_eff:.1f} cm^-1")
print(f"\n  the CD model currently ASSUMES  S = 0.35, omega = 1450 cm^-1")
json.dump({'lambda_cm': lam, 'S_total': float(S.sum()), 'S_eff_highfreq': S_eff,
           'omega_eff_cm': w_eff, 'basis': '6-31G*',
           'modes': [{'omega_cm': float(freq[i]), 'S': float(S[i])} for i in np.argsort(freq)]},
          open('results_huang_rhys.json','w'), indent=2)
print("\nwrote results_huang_rhys.json")

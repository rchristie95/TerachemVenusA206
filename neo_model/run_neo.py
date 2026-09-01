#!/usr/bin/env python3
"""
Proton-quantum effect on the FULL anionic chromophore pi->pi* (CT) excitation,
via the His148-ND1-H...O(phenolate) bridge.

Model: full anionic CR2 chromophore (phenolate ring + methine bridge +
imidazolinone) + His148 imidazole, from the QM/MM-relaxed protein geometry.

Same heavy-atom geometry + same electronic functional/basis both ways:
  (A) classical proton : standard DFT + TDDFT
  (B) quantum  proton  : NEO-DFT (HD1 quantum) + NEO-TDDFT
NEO-TDDFT (davidson) does not support epc, so epc is OFF in both the NEO
ground state and response (mean-field quantum-proton estimate). NEO-TDDFT
mixes proton-vibrational + electronic roots; the electronic pi->pi* is matched
to (A) by energy. (B)-(A) = mean-field proton-quantum shift on the excitation.
"""
import json, numpy as np
XYZ = "chromo_his.xyz"
idx = json.load(open("idx_full.json"))
HD1, O, CZ, ND1 = idx["HD1"], idx["O"], idx["CZ"], idx["ND1"]
XC, EBAS, NBAS = "b3lyp", "6-31g*", "pb4d"
NST_A, NST_B = 8, 16
EV, NM = 27.211386, 1239.841984

L = open(XYZ).read().splitlines(); n = int(L[0].split()[0])
atomstr = "\n".join(L[2:2+n])
xyz = np.array([[float(c) for c in r.split()[1:4]] for r in L[2:2+n]])
print(f"# {n} atoms  C-O={np.linalg.norm(xyz[CZ]-xyz[O]):.3f}  "
      f"N-H={np.linalg.norm(xyz[ND1]-xyz[HD1]):.3f}  H..O={np.linalg.norm(xyz[HD1]-xyz[O]):.3f}", flush=True)

# ---------- (A) classical proton ----------
print("\n=== (A) classical-proton DFT + TDDFT ===", flush=True)
from pyscf import gto, dft, tddft
m = gto.M(atom=atomstr, basis=EBAS, charge=-1, spin=0, verbose=0, max_memory=8000)
mf = dft.RKS(m); mf.xc = XC; mf.kernel()
td = tddft.TDDFT(mf); td.nstates = NST_A; td.kernel()
osc = td.oscillator_strength(); eA = np.asarray(td.e) * EV
iA = int(np.argmax(osc)); E_bright = eA[iA]
print(f"  E(DFT)={mf.e_tot:.6f}")
for k in range(len(eA)):
    print(f"    S{k+1}: {eA[k]:.3f} eV  {NM/eA[k]:6.1f} nm  f={osc[k]:.3f}"
          + ("   <== bright pi->pi*" if k == iA else ""))
print(f"  CLASSICAL bright: {E_bright:.3f} eV = {NM/E_bright:.1f} nm  (f={osc[iA]:.3f})", flush=True)

# ---------- (B) quantum proton (NEO, no epc) ----------
print("\n=== (B) quantum-proton NEO-DFT + NEO-TDDFT (HD1 quantum, no epc) ===", flush=True)
from pyscf import neo
nmol = neo.M(atom=atomstr, basis=EBAS, nuc_basis=NBAS, charge=-1, spin=0,
             quantum_nuc=[HD1], verbose=0, max_memory=8000)
nmf = neo.KS(nmol, xc=XC)          # epc left as None
nmf.scf()
print(f"  E(NEO-DFT)={nmf.e_tot:.6f}  converged={nmf.converged}")
# proton expectation position (trace-normalised)
try:
    for t, comp in nmf.components.items():
        if t == "e": continue
        cmol = comp.mol
        dm = comp.make_rdm1()
        if isinstance(dm, np.ndarray) and dm.ndim == 3: dm = dm.sum(0)
        S = cmol.intor("int1e_ovlp"); norm = float(np.einsum("ij,ji->", dm, S))
        r = cmol.intor("int1e_r")    # lab-frame, origin 0
        # int1e_r is in BOHR (pyscf internal units); xyz[] is in ANGSTROM.
        BOHR_TO_ANGSTROM = 0.5291772109
        pos = np.einsum("xij,ji->x", r, dm).real / norm * BOHR_TO_ANGSTROM
        print(f"  quantum-proton norm={norm:.3f}  <r>=[{pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f}]")
        print(f"    classical HD1=[{xyz[HD1][0]:.3f} {xyz[HD1][1]:.3f} {xyz[HD1][2]:.3f}]"
              f"  shift={np.linalg.norm(pos-xyz[HD1]):.3f} A toward "
              f"({'phenolate O' if np.linalg.norm(pos-xyz[O])<np.linalg.norm(xyz[HD1]-xyz[O]) else 'His N'})")
        print(f"    <r>..O={np.linalg.norm(pos-xyz[O]):.3f}  (classical H..O={np.linalg.norm(xyz[HD1]-xyz[O]):.3f})")
        break
except Exception as e:
    print("  (proton <r> unavailable:", repr(e), ")")

td2 = neo.CTDDFT(nmf); td2.nstates = NST_B
td2.kernel()
eB = np.asarray(td2.e) * EV
print("  NEO-TDDFT roots (proton-vibrational + electronic, mixed):")
for k in range(len(eB)):
    print(f"    R{k+1}: {eB[k]:.3f} eV  {NM/max(eB[k],1e-3):7.1f} nm")
j = int(np.argmin(np.abs(eB - E_bright)))
print(f"  matched electronic root R{j+1}: {eB[j]:.3f} eV = {NM/eB[j]:.1f} nm")
print(f"\n=== PROTON-QUANTUM SHIFT (B-A) = {eB[j]-E_bright:+.3f} eV  "
      f"({NM/eB[j]-NM/E_bright:+.1f} nm) ===")
print(f"    classical {NM/E_bright:.1f} nm -> quantum {NM/eB[j]:.1f} nm")

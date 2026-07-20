#!/usr/bin/env python3
"""
terachem_tddft_analysis.py

Performs TD-DFT analysis on the optimized geometry from the hybrid pipeline.
1. Calculates first 15 excited states.
2. Identifies the brightest state in the 450-600 nm range.
3. Generates Transition and Difference densities for that specific state.
"""

import os
import sys
import re
import shutil
import subprocess
from pathlib import Path
import numpy as np

# Configuration
# WORKDIR = Path("tc_tddft_analysis_big")
# INPUT_DIR = Path("tc_hybrid_amoeba_big")
WORKDIR = Path("tc_tddft_analysis_lbfgs")
INPUT_DIR = Path("tc_hybrid_lbfgs")

GEOM_FILE = "final_absorption.xyz"
CHARGES_FILE = "mm_charges.dat"

TC_METHOD = "wb97xd3"
TC_BASIS = "6-311g**"
# TC_METHOD = "hf"
# TC_BASIS = "3-21g"
TC_CHARGE = -1
TC_SPIN = 1
NUM_STATES = 20
TC_PATH = os.environ.get("TC_PATH", "terachem")

def setup_workspace():
    if not WORKDIR.exists():
        WORKDIR.mkdir()
    
    # Copy necessary files
    src_geom = INPUT_DIR / GEOM_FILE
    if not src_geom.exists():
        src_geom = INPUT_DIR / "qm_step.xyz"
    if not src_geom.exists():
        src_geom = INPUT_DIR / "qm_1.xyz"
    
    if not src_geom.exists():
        print(f"[!] Geometry file not found in {INPUT_DIR}")
        sys.exit(1)
        
    shutil.copy(src_geom, WORKDIR / "geometry.xyz")
    
    src_charges = INPUT_DIR / CHARGES_FILE
    if src_charges.exists():
        shutil.copy(src_charges, WORKDIR / "mm_charges.dat")
    else:
        print("[!] Warning: MM charges file not found. Running Gas Phase?")

    return WORKDIR / "geometry.xyz", WORKDIR / "mm_charges.dat"

def run_tddft_energy(geom_path, charges_path):
    print("[*] Running TD-DFT Energy Calculation...")
    inp_file = WORKDIR / "energy.in"
    scr_dir = WORKDIR / "scr_energy"
    out_file = WORKDIR / "energy.out"
    
    with open(inp_file, 'w') as f:
        f.write(f"coordinates {geom_path.name}\n")
        f.write(f"run energy\n")
        f.write("cis yes\n") 
        f.write(f"basis {TC_BASIS}\n")
        f.write(f"method {TC_METHOD}\n")
        f.write(f"charge {TC_CHARGE}\n")
        f.write(f"spinmult {TC_SPIN}\n")
        if charges_path.exists():
            f.write(f"pointcharges {charges_path.name}\n")
        f.write(f"scrdir {scr_dir.name}\n")
        f.write(f"cisnumstates {NUM_STATES}\n")
        f.write("cismaxiter 200\n")
        f.write("cismax 500\n")
        f.write("scf diis+a\n")
        f.write("threall 1.0e-13\n")
        f.write("end\n")

    cmd = [TC_PATH, inp_file.name]
    print(f"    - Executing: {' '.join(cmd)}")
    with open(out_file, 'w') as log:
        subprocess.run(cmd, cwd=WORKDIR, stdout=log, stderr=subprocess.STDOUT)
        
    return out_file

def parse_excitation_energies(out_file):
    print("[*] Parsing Excited States...")
    states = {} 
    
    if not out_file.exists(): return []
    text = out_file.read_text(errors="replace")
    
    # 1. Parse Energies and Oscillator Strengths
    in_table = False
    for line in text.splitlines():
        if "Final Excited State Results" in line: in_table = True; continue
        if in_table:
            if "---" in line or not line.strip(): continue
            if "Printing MM field" in line: break
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                try:
                    r = int(parts[0])
                    ev = float(parts[2])
                    osc = float(parts[3])
                    nm = 1239.84193 / ev if ev > 0 else 0.0
                    states[r] = {'root': r, 'nm': nm, 'osc': osc, 'ev': ev}
                except: continue

    # 2. Parse Transition Dipole Moments (magnitude |T|)
    in_dipole = False
    for line in text.splitlines():
        if "Transition dipole moments:" in line and "between" not in line:
            in_dipole = True; continue
        if "Transition dipole moments between" in line:
            in_dipole = False; continue
            
        if in_dipole:
            parts = line.split()
            # Look for lines starting with integer (Root index)
            if len(parts) >= 5 and parts[0].isdigit():
                try:
                    r = int(parts[0])
                    if r in states:
                        states[r]['mu_mag'] = float(parts[4])
                except: continue

    state_list = sorted(states.values(), key=lambda x: x['root'])

    print(f"    - Found {len(state_list)} excited states.")
    print(f"    - {'Root':<5} {'Wavelength':<10} {'Osc. Str.':<10} {'|mu| (a.u.)'}")
    print("    " + "-"*45)
    for s in state_list:
        mu_val = s.get('mu_mag', 0.0)
        print(f"    - {s['root']:<5} {s['nm']:.2f} nm     {s['osc']:.4f}     {mu_val:.4f}")
        
    return state_list

def select_brightest_state(states, nm_min=450, nm_max=600):
    candidates = [s for s in states if nm_min <= s['nm'] <= nm_max]
    if not candidates:
        if states: candidates = states
        else: return None
    best = max(candidates, key=lambda x: x['osc'])
    print(f"\n[*] Selected Target State: Root {best['root']}")
    print(f"    - Wavelength: {best['nm']:.2f} nm")
    print(f"    - Osc. Strength: {best['osc']:.4f}")
    return best['root']

def generate_densities(geom_path, charges_path, root):
    print(f"\n[*] Generating Densities for Root {root}...")
    inp_file = WORKDIR / "plot.in"
    scr_dir = WORKDIR / "scr_plot"
    out_file = WORKDIR / "plot.out"
    target_root = max(int(root), 1)
    # Solve at least as many CIS states as the energy run (NUM_STATES) so the
    # Davidson eigensolver converges to the SAME state set/ordering. A plot run
    # with cisnumstates=target_root can lock onto a DIFFERENT state in a dense
    # manifold and plot the wrong transition density -> the renorm factor leaves
    # ~1.0 and the downstream coupling J is computed from the wrong state.
    cis_states = max(target_root, int(NUM_STATES))

    with open(inp_file, 'w') as f:
        f.write(f"coordinates {geom_path.name}\n")
        f.write("run energy\n")
        f.write("cis yes\n")
        f.write(f"basis {TC_BASIS}\nmethod {TC_METHOD}\ncharge {TC_CHARGE}\nspinmult {TC_SPIN}\n")
        if charges_path.exists(): f.write(f"pointcharges {charges_path.name}\n")
        f.write(f"scrdir {scr_dir.name}\ncisnumstates {cis_states}\n")
        f.write("cismaxiter 200\ncismax 500\n")
        f.write("scf diis+a\nthreall 1.0e-13\n")
        f.write("cisdiffdensity yes\ncistransdensity yes\n")
        f.write(f"cistarget {target_root}\n")
        f.write("end\n")

    with open(out_file, 'w') as log:
        subprocess.run([TC_PATH, inp_file.name], cwd=WORKDIR, stdout=log, stderr=subprocess.STDOUT)
    
    diff_src = scr_dir / f"diffdens_{target_root}.dx"
    trans_src = scr_dir / f"transdens_{target_root}.dx"
    if diff_src.exists(): shutil.copy(diff_src, WORKDIR / f"abs_diffdens_{target_root}.dx")
    if trans_src.exists(): shutil.copy(trans_src, WORKDIR / f"abs_transdens_{target_root}.dx")

def main():
    geom, charges = setup_workspace()
    if (WORKDIR / "energy.out").exists():
        print("[*] Using existing energy.out...")
        out_file = WORKDIR / "energy.out"
    else:
        out_file = run_tddft_energy(geom, charges)
    states = parse_excitation_energies(out_file)
    target_root = select_brightest_state(states)
    if target_root:
        if not (WORKDIR / f"abs_transdens_{target_root}.dx").exists():
            generate_densities(geom, charges, target_root)
        else:
            print(f"[*] Using existing densities for Root {target_root}...")
    print(f"\n[Done] Results in {WORKDIR}")

if __name__ == "__main__": main()

#!/usr/bin/env python3
"""Transition-dipole basis ladder on the production frames.

|mu| sets every absolute CD intensity in this work, because the couplet
amplitude scales as |mu|^4 (J ~ mu^2 and the chirality triple product ~ mu^2).
Production runs at def2-SVP and gives 9.40 D; two independent experimental
routes indicate 7.6-7.9 D. The 23% double-to-triple-zeta shift quoted in the
manuscript comes from a formaldehyde series, not from this chromophore in this
embedding, so the ladder below measures it directly.

    def2-SVP    already on disk (production)
    def2-SVPD   adds diffuse functions at fixed double zeta
    def2-TZVPD  adds zeta at fixed diffuse
    def2-QZVPD  convergence check, a few frames only

SVP -> SVPD isolates the diffuse-function effect and SVPD -> TZVPD the zeta
effect, which the single SVP -> TZVP comparison would confound.

Runs inside the EXISTING prepared frame directories, reusing each frame's own
geometry.xyz and mm_charges.dat, and writes to its own td_mu_<basis>.in and
scr_mu_<basis>/ so nothing from the production sweep is touched. Stages run in
priority order and results are written incrementally, so a partial run is still
usable.

State identity: unlike the TDX dark barrel, site A and site B each have one
dominant bright state carrying ~85% of the oscillator strength in the 7-root
window, so the brightest root is a safe criterion here. Every root is recorded
anyway so the choice is auditable.
"""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from parse_terachem_tddft import parse as parse_tc                 # noqa: E402

ROOT = HERE.parent
TERACHEM = Path("/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem")
WINDOW = "v5camp"
# Measured cost at 44 atoms: SVPD ~384 s/site. TZVPD and QZVPD scale steeply,
# so the ensemble-level question (how much do diffuse functions move |mu|?) is
# answered on many frames at SVPD, and the convergence question (does |mu| come
# down with zeta?) on a few frames through the full ladder. Stages run in this
# order and write incrementally, so a partial night still answers the first.
# QZVPD reduced to one frame: at ~2400 basis functions it is a convergence
# indicator, not a statistic, and may not fit the card at all.
STAGES = [("svpd", 24), ("tzvpd", 8), ("qzvpd", 1)]
OUT_DIR = ROOT / "results"


def stratified(n: int) -> list[int]:
    """n frames spread evenly across the |Delta| distribution."""
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/{WINDOW}_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: abs(float(r["site_A_energy_cm-1"])
                                - float(r["site_B_energy_cm-1"])))
    idx = np.linspace(0, len(rows) - 1, n).round().astype(int)
    return [int(rows[i]["frame_index"]) for i in idx]


def run_one(frame: int, site: str, basis: str, gpu: str) -> dict | None:
    d = ROOT / "frames" / WINDOW / f"frame_{frame:04d}" / f"site_{site}"
    out = d / f"td_mu_{basis}.out"
    inp = d / f"td_mu_{basis}.in"
    if not out.exists() or "Job finished:" not in out.read_text(errors="replace"):
        shutil.copy2(ROOT / f"input_templates/tddft_mu_{basis}.in", inp)
        t0 = time.time()
        with out.open("w") as fh:
            rc = subprocess.run([str(TERACHEM), inp.name], cwd=d, stdout=fh,
                                stderr=subprocess.STDOUT,
                                env={**__import__("os").environ,
                                     "CUDA_VISIBLE_DEVICES": gpu}).returncode
        runtime = time.time() - t0
        shutil.rmtree(d / f"scr_mu_{basis}", ignore_errors=True)
        if rc:
            return {"frame": frame, "site": site, "basis": basis,
                    "status": f"terachem rc={rc}", "runtime_s": runtime}
    else:
        runtime = float("nan")
    parsed = parse_tc(out)
    if parsed["status"] != "converged":
        return {"frame": frame, "site": site, "basis": basis,
                "status": parsed["status"], "runtime_s": runtime}
    roots = sorted(parsed["roots"], key=lambda r: r["root"])
    best = max(roots, key=lambda r: r["oscillator_strength"])
    mu = float(np.linalg.norm(best["transition_dipole_au"])) * 2.5417464519
    total_f = sum(r["oscillator_strength"] for r in roots)
    return {"frame": frame, "site": site, "basis": basis, "status": "complete",
            "root": best["root"], "energy_eV": best["energy_eV"],
            "f": best["oscillator_strength"], "mu_D": mu,
            "f_fraction_of_window": best["oscillator_strength"] / total_f,
            "runtime_s": runtime,
            "all_roots": [{"root": r["root"], "E_eV": r["energy_eV"],
                           "f": r["oscillator_strength"]} for r in roots]}


def main() -> None:
    gpu = sys.argv[1] if len(sys.argv) > 1 else "0"
    tag = sys.argv[2] if len(sys.argv) > 2 else "a"
    nworkers = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    offset = "abc".index(tag)
    OUT = OUT_DIR / f"basis_ladder_{tag}.json"
    results = json.loads(OUT.read_text()) if OUT.exists() else []
    done = {(r["frame"], r["site"], r["basis"]) for r in results
            if r.get("status") == "complete"}
    wanted = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    for basis, n in STAGES:
        if wanted and basis not in wanted:
            continue
        jobs = [(f, st) for f in stratified(n) for st in ("A", "B")]
        jobs = jobs[offset::nworkers]
        print(f"[stage] def2-{basis}: worker {tag} takes {len(jobs)} of "
              f"{2*n} site-frames", flush=True)
        for frame, site in jobs:
                if (frame, site, basis) in done:
                    continue
                r = run_one(frame, site, basis, gpu)
                if r is None:
                    continue
                results.append(r)
                OUT.write_text(json.dumps(results, indent=2) + "\n")
                if r["status"] == "complete":
                    print(f"  {basis} frame {frame} site {site}: root {r['root']} "
                          f"{r['energy_eV']:.4f} eV  f={r['f']:.3f} "
                          f"|mu|={r['mu_D']:.3f} D  ({r['runtime_s']:.0f} s)",
                          flush=True)
                else:
                    print(f"  {basis} frame {frame} site {site}: {r['status']}",
                          flush=True)
        ok = [r for r in results if r["basis"] == basis and r["status"] == "complete"]
        if ok:
            m = np.array([r["mu_D"] for r in ok])
            print(f"[stage] def2-{basis} done: <|mu|> = {m.mean():.3f} +/- "
                  f"{m.std(ddof=1):.3f} D over {m.size} site-frames", flush=True)


if __name__ == "__main__":
    main()

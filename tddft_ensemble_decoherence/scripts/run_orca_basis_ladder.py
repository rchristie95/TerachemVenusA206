#!/usr/bin/env python3
"""The basis-ladder rungs TeraChem cannot run, in ORCA.

TeraChem's basis format is S/P/D only, so def2-TZVPD (26 F-shells on this
44-atom region) dies silently right after basis setup. That is a hard limit, not
a resource or guess problem -- both were chased. def2-SVPD has no F-shells,
which is why it ran.

The manuscript now asserts that the linear-response transition dipole is
basis-converged, on the strength of two numbers: 0.3% from def2-SVP to def2-SVPD
across 48 TeraChem site-frames, and 0.2% from def2-SVPD to def2-TZVPD in
dipole_gauge/GAUGE_VERDICT.md. The second is a single geometry with CAM-B3LYP,
not the production functional and not an ensemble. This runs BOTH rungs in ORCA,
same code and same wB97X-D3, so the SVPD -> TZVPD ratio is internal rather than
cross-code, and checks that assertion on production frames.

Reuses each frame's existing geometry.xyz and mm_charges.dat, converted to
ORCA's .pc format the same way run_steom_v5.write_inputs does (TeraChem writes
N / comment / q x y z; ORCA wants N / q x y z). TDA is set to match production.
"""
from __future__ import annotations

import csv, json, os, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
ORCA_DIR = REPO / "orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
ORCA = ORCA_DIR / "orca"
OUT = ROOT / "results" / "orca_basis_ladder.json"
DEBYE_PER_AU = 2.5417464519

TEMPLATE = """! wB97X-D3 {basis} RIJCOSX def2/J TightSCF
# Basis rung TeraChem cannot run: def2-TZVPD needs F functions and TeraChem's
# basis format is S/P/D only. Same geometry and point-charge field the
# production TeraChem job was given for this frame and site.
%maxcore {maxcore}
%pal nprocs {nprocs} end
%pointcharges "field.pc"
%tddft
  nroots 7
  tda true
  maxdim 10
end
* xyzfile -1 1 geometry.xyz
"""


def stratified(n):
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/v5camp_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: abs(float(r["site_A_energy_cm-1"])
                                - float(r["site_B_energy_cm-1"])))
    idx = np.linspace(0, len(rows) - 1, n).round().astype(int)
    return [int(rows[i]["frame_index"]) for i in idx]


def prepare(frame, site, basis, nprocs, maxcore):
    src = ROOT / "frames/v5camp" / f"frame_{frame:04d}" / f"site_{site}"
    d = ROOT / "orca_ladder" / f"frame_{frame:04d}" / f"site_{site}_{basis}"
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / "geometry.xyz", d / "geometry.xyz")
    lines = (src / "mm_charges.dat").read_text().splitlines()
    n = int(lines[0].split()[0])
    body = lines[2:2 + n]
    if len(body) != n:
        raise RuntimeError(f"{src}: header says {n}, found {len(body)}")
    (d / "field.pc").write_text(f"{n}\n" + "\n".join(body) + "\n")
    (d / "orca.inp").write_text(TEMPLATE.format(
        basis=basis, nprocs=nprocs, maxcore=maxcore))
    return d


def parse(path):
    """Lowest bright root: energy in eV, oscillator strength, |mu| in Debye."""
    text = path.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    rows, block = [], False
    for line in text.splitlines():
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in line:
            block, hdr = True, 0
            continue
        if block:
            p = line.split()
            if len(p) >= 7 and p[0].replace("-", "").isdigit() is False and "->" in line:
                try:
                    ev = float(p[3]); f = float(p[6])
                    rows.append((ev, f))
                except (ValueError, IndexError):
                    pass
            elif rows and not line.strip():
                break
    if not rows:
        return None
    ev, f = max(rows, key=lambda r: r[1])
    mu = np.sqrt(1.5 * f / (ev / 27.211386)) * DEBYE_PER_AU
    return {"energy_eV": ev, "f": f, "mu_D": float(mu), "n_roots": len(rows)}


def main():
    frames = stratified(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
    bases = (sys.argv[2].split(",") if len(sys.argv) > 2
             else ["def2-SVPD", "def2-TZVPD"])
    nprocs, maxcore = 16, 3000
    results = json.loads(OUT.read_text()) if OUT.exists() else []
    done = {(r["frame"], r["site"], r["basis"]) for r in results if r.get("mu_D")}
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = ":".join([
        str(ORCA_DIR), "/home/robson/anaconda3/envs/openmpi416/lib",
        env.get("LD_LIBRARY_PATH", "")]).strip(":")
    env["PATH"] = "/home/robson/anaconda3/envs/openmpi416/bin:" + env["PATH"]

    for basis in bases:
        for frame in frames:
            for site in ("A", "B"):
                if (frame, site, basis) in done:
                    continue
                d = prepare(frame, site, basis, nprocs, maxcore)
                out = d / "orca.out"
                print(f"[start] {basis} frame {frame} site {site}", flush=True)
                t0 = time.time()
                with out.open("w") as fh:
                    subprocess.run([str(ORCA), "orca.inp"], cwd=d, stdout=fh,
                                   stderr=subprocess.STDOUT, env=env)
                dt = time.time() - t0
                for pat in ("*.tmp*", "*.gbw", "*.densities", "*.cis"):
                    for f in d.glob(pat):
                        f.unlink(missing_ok=True)
                r = parse(out) or {}
                r.update(frame=frame, site=site, basis=basis, runtime_h=dt / 3600)
                results.append(r)
                OUT.write_text(json.dumps(results, indent=2) + "\n")
                if r.get("mu_D"):
                    print(f"[done]  {basis} frame {frame} site {site}: "
                          f"{r['energy_eV']:.4f} eV  f={r['f']:.4f}  "
                          f"|mu|={r['mu_D']:.3f} D  ({dt/3600:.2f} h)", flush=True)
                else:
                    print(f"[FAIL]  {basis} frame {frame} site {site} "
                          f"({dt/3600:.2f} h)", flush=True)
        ok = [r for r in results if r["basis"] == basis and r.get("mu_D")]
        if ok:
            m = np.array([r["mu_D"] for r in ok])
            print(f"[stage] {basis}: <|mu|> = {m.mean():.3f} +/- "
                  f"{m.std(ddof=1) if m.size > 1 else 0:.3f} D over {m.size}",
                  flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""STEOM-DLPNO-CCSD references for the v5 campaign.

Four production frames, uniformly strided every 500 DCD frames from the start of
production (900, 1400, 1900, 2400), two sites each. These are the absolute
energy/dipole scale that the 500 def2-SVP TDDFT frames are calibrated against,
so they must see EXACTLY the geometry and embedding field that TDDFT sees.

That is guaranteed structurally rather than by inspection: the frames are
prepared by the same prepare_production_frame.py invocation run_ensemble.py
builds for the v5camp window, into a separate output directory. The two were
checked byte-for-byte identical on frame 900 before this script was written.
The separate directory is required -- run_ensemble.py treats a non-empty frame
directory without preparation.json as a corrupt frame and refuses it, so the
STEOM stream must not write into frames/v5camp/.

ORCA is CPU-bound and TeraChem holds the GPU, so this stream overlaps the TDDFT
ensemble rather than adding to its wall clock.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(REPO / "weekend_dipole"))

ORCA_DIR = REPO / "orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
ORCA = ORCA_DIR / "orca"
TOPOLOGY = REPO / "tc_tandem_v5_campaign/solvated_protonated.pdb"
TRAJECTORY = REPO / "tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd"
CACHE = ROOT / "cache/embedding_charges_v5.npz"
CR2_PRMTOP = REPO / "anionic_build/monomer_solv.prmtop"
OUT = ROOT / "steom_v5"

# Uniform stride 500 over the production block (DCD frames 900-2899), as
# specified. Every one of these is also in the TDDFT manifest (900 + 4k), so
# each STEOM point has a same-geometry TDDFT partner.
FRAMES = [900, 1400, 1900, 2400]
SITES = ["A", "B"]

# Same %mdci block as weekend_dipole/steom_def2_svp.inp, which converged in
# 2 h 09 m. TCutPNOSingles/NDav/MaxIter are the documented recipe for getting
# parallel EOM to converge at all; do not drop them.
TEMPLATE = """! STEOM-DLPNO-CCSD def2-SVP AutoAux RIJCOSX TightSCF TightPNO
# v5 campaign STEOM reference, frame {frame} site {site}.
# Geometry and point charges are the identical files TeraChem is given for the
# same frame under frames/v5camp/frame_{frame:04d}/site_{site}/.
%maxcore {maxcore}
%pal nprocs {nprocs} end
%pointcharges "{pc}"
%mdci
  NRoots 5
  DoRootWise true
  FollowCIS true
  RootHoming true
  CheckEachRoot true
  DTol 1e-5
  DoStoreSTEOM true
  AddL2Term True
  DoDbFilter true
  NRootsCISNAT 10
  TCutPNOSingles 1e-12
  NDav 400
  MaxIter 400
end
* xyzfile -1 1 {xyz}
"""


def log(msg: str) -> None:
    print(f"[{time.strftime('%F %H:%M:%S')}] {msg}", flush=True)


def prepare(frame: int) -> Path:
    """Run the v5camp preparation for one frame into steom_v5/frame_XXXX/."""
    d = OUT / f"frame_{frame:04d}"
    if (d / "preparation.json").exists():
        return d
    if d.exists():
        shutil.rmtree(d)
    cmd = [
        sys.executable, str(ROOT / "scripts/prepare_production_frame.py"),
        "--topology", str(TOPOLOGY),
        "--trajectory", str(TRAJECTORY),
        "--frame-index", str(frame),
        "--site-b-offset", "8",
        "--output-dir", str(d),
        "--embedding-cache", str(CACHE),
        "--retain-partner-cr2-charges",
        "--amber-cr2-prmtop", str(CR2_PRMTOP),
        "--full-system-embedding", "--link-only-exclusion",
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return d


def write_inputs(frame: int, site: str, nprocs: int, maxcore: int) -> Path:
    """Convert the TeraChem frame files into an ORCA job in place.

    TeraChem's point-charge file is `N / comment / q x y z`; ORCA's .pc is the
    same without the comment line. Nothing is reordered or rescaled, so the two
    codes see the same field to the printed precision.
    """
    d = OUT / f"frame_{frame:04d}" / f"site_{site}"
    lines = (d / "mm_charges.dat").read_text().splitlines()
    n = int(lines[0].split()[0])
    body = lines[2:2 + n]
    if len(body) != n:
        raise RuntimeError(f"{d}: header says {n} charges, found {len(body)}")
    (d / "field.pc").write_text(f"{n}\n" + "\n".join(body) + "\n")
    (d / "steom.inp").write_text(TEMPLATE.format(
        frame=frame, site=site, nprocs=nprocs, maxcore=maxcore,
        pc="field.pc", xyz="geometry.xyz"))
    return d


def free_ram_gb() -> int:
    """Available RAM. `free -g` column 7 -- NOT disk.

    Kept deliberately separate from free_disk_gb(). Conflating the two cost an
    hour on 2026-08-31: the completion pass was handed a 300 GB threshold that
    was reasoning about scratch, compared it against 62 GB of total RAM, and sat
    in its hold loop waiting for a number this machine cannot reach.
    """
    out = subprocess.run(["free", "-g"], capture_output=True, text=True).stdout
    return int(out.splitlines()[1].split()[6])


def free_disk_gb() -> int:
    """Free space on the filesystem holding the scratch."""
    st = os.statvfs(OUT)
    return int(st.f_bavail * st.f_frsize / 1024 ** 3)


def clean_scratch(d: Path) -> None:
    """DLPNO scratch has no restart value and is large. Never touch densities."""
    for p in list(d.glob("steom*.tmp*")) + list(d.glob("steom.bas[0-9]*")) + \
             list(d.glob("steom.hostnames")):
        if "densit" in p.name.lower():
            continue
        try:
            p.unlink()
        except OSError:
            pass


def run_one(job: tuple[int, str], nprocs: int, maxcore: int, hours: float,
            min_free_gb: int, min_disk_gb: int) -> dict:
    frame, site = job
    d = OUT / f"frame_{frame:04d}" / f"site_{site}"
    tag = f"frame {frame} site {site}"
    if (d / "steom.out").exists() and "ORCA TERMINATED NORMALLY" in \
            (d / "steom.out").read_text(errors="replace"):
        log(f"skip  {tag} (already complete)")
        return {"frame": frame, "site": site, "status": "complete"}

    waited = 0
    while waited < 3600:
        ram, disk = free_ram_gb(), free_disk_gb()
        if ram >= min_free_gb and disk >= min_disk_gb:
            break
        log(f"hold  {tag}: {ram} GB RAM (need {min_free_gb}), "
            f"{disk} GB disk (need {min_disk_gb})")
        time.sleep(120)
        waited += 120

    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = ":".join([
        str(ORCA_DIR), "/home/robson/anaconda3/envs/openmpi416/lib",
        env.get("LD_LIBRARY_PATH", "")]).strip(":")
    env["PATH"] = "/home/robson/anaconda3/envs/openmpi416/bin:" + env["PATH"]

    log(f"start {tag} ({nprocs} ranks, {free_ram_gb()} GB RAM, {free_disk_gb()} GB disk)")
    t0 = time.time()
    with (d / "steom.out").open("w") as fh:
        rc = subprocess.run(["timeout", f"{hours}h", str(ORCA), "steom.inp"],
                            cwd=d, stdout=fh, stderr=subprocess.STDOUT,
                            env=env).returncode
    dt = time.time() - t0
    clean_scratch(d)
    text = (d / "steom.out").read_text(errors="replace")
    ok = "ORCA TERMINATED NORMALLY" in text
    status = "complete" if ok else ("timeout" if rc == 124 else f"failed rc={rc}")
    log(f"{'OK   ' if ok else 'FAIL '} {tag} in {dt / 3600:.2f} h ({status})")
    result = {"frame": frame, "site": site, "status": status,
              "runtime_hours": dt / 3600}
    if ok:
        try:
            from parse_gauge import blocks_after_leftright
            length, _ = blocks_after_leftright(text.splitlines())
            if length:
                bright = max(length, key=lambda r: r["f"])
                result.update(energy_eV=bright["eV"], f=bright["f"],
                              mu_D=(bright["D2"] ** 0.5) * 2.5417464,
                              root=bright["root"])
                log(f"      {tag}: E={bright['eV']:.4f} eV  f={bright['f']:.4f}"
                    f"  |mu|={result['mu_D']:.3f} D")
        except Exception as exc:      # parsing must never lose a finished job
            result["parse_error"] = repr(exc)
    (d / "steom_result.json").write_text(json.dumps(result, indent=2))
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    # Defaults sized to the measured footprint, not to the core count. A
    # DLPNO-STEOM run of this class writes 250-300 GB of PNO/PAO scratch
    # (weekend_dipole/run_ddterm_chain.sh:59) and that is set by the QM basis
    # and the 134 correlated electrons, NOT by the point-charge count -- so it
    # cannot be reduced by trimming the embedding. On a 437 GB disk that permits
    # exactly one job at a time. This defaulted to 4 on 2026-08-30; the burn hit
    # 699 GB/h once the first job reached the EOM dressing and a frame had to be
    # killed to recover.
    ap.add_argument("--concurrency", type=int, default=1)
    ap.add_argument("--nprocs", type=int, default=10)
    ap.add_argument("--maxcore", type=int, default=3000)
    ap.add_argument("--hours", type=float, default=10.0)
    ap.add_argument("--min-free-gb", type=int, default=12,
                    help="Minimum available RAM in GB before a job starts.")
    ap.add_argument("--min-disk-gb", type=int, default=150,
                    help="Minimum free DISK in GB before a job starts; one "
                         "job peaks near 86 GB of scratch at def2-SVP.")
    ap.add_argument("--frames", type=int, nargs="*", default=FRAMES)
    args = ap.parse_args()

    OUT.mkdir(exist_ok=True)
    log(f"===== v5 STEOM stream: frames {args.frames} x sites {SITES} =====")
    jobs = []
    for frame in args.frames:
        prepare(frame)
        for site in SITES:
            write_inputs(frame, site, args.nprocs, args.maxcore)
            jobs.append((frame, site))
    log(f"prepared {len(jobs)} jobs, running {args.concurrency} at a time")

    with ThreadPoolExecutor(max_workers=args.concurrency) as pool:
        results = list(pool.map(
            lambda j: run_one(j, args.nprocs, args.maxcore, args.hours,
                              args.min_free_gb, args.min_disk_gb), jobs))

    (OUT / "steom_v5_summary.json").write_text(json.dumps(results, indent=2))
    ok = sum(r["status"] == "complete" for r in results)
    log(f"===== v5 STEOM stream done: {ok}/{len(results)} complete =====")
    for r in results:
        if r["status"] == "complete" and "energy_eV" in r:
            log(f"  frame {r['frame']} site {r['site']}: "
                f"{r['energy_eV']:.4f} eV  |mu|={r['mu_D']:.3f} D")
    sys.exit(0 if ok == len(results) else 1)


if __name__ == "__main__":
    main()

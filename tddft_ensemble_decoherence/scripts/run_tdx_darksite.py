#!/usr/bin/env python3
"""Second QM region for the dVenus-TDX campaign: Delta and J against barrel 2.

WHY. The TDX site-A sweep (500 frames, complete) measures one chromophore. The
question it cannot answer is what the OTHER barrel does, and therefore what
Delta and J are for this construct -- which is what decides whether the inferred
CD spectrum has a bisignate couplet at all.

The barrel-2 region is the same physical slot as site A, asymmetric only in
content: the deposited construct deletes the chromophore tyrosine, so the mature
CR2 is the un-matured Gly-Gly, while the pi-stacked tyrosine Venus wraps around
the chromophore is still there (residue 473, not 464: collapsing CR2 -- one PDB
residue standing for three amino acids -- into two glycines shifts everything
downstream by one).

    site A       CR2(29)    + phenol(12) + 3 link H = 44 atoms at -1 e
    site Bdark   GlyGly(14) + phenol(12) + 3 link H = 29 atoms at  0 e

Two things this pass must get right, both verified on frame 100 before launch:

  * the live chromophore is a PARTNER for this site, sitting in its embedding
    field 26 A away. load_or_build_charges zeroes CR2 by name, so without
    --retain-partner-cr2-charges the field sums to +0.849200 e with 29 charges
    pinned at exactly 0.0 -- the dark site polarised by a field missing -1 e.
    With the flag it sums to -0.150800 e and nothing is pinned. This is the same
    error that ran on five v5camp calibration frames.

  * the region is NEUTRAL. The input templates hard-code "charge -1", the mature
    anionic chromophore; process_site takes the charge from the preparation
    summary instead and records it in the method profile.

State selection differs by necessity: there is no CR2 density to track against,
so the brightest computed root is taken. That is a different criterion from site
A's density-overlap tracking and is labelled as such in every output.

Site A is never touched: this pass writes only site_Bdark/ inside the existing
frame directories, and the TDC reads site A's already-stored root densities.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_ensemble as R                                           # noqa: E402

ROOT = R.ROOT
REPO = R.REPO
WINDOW = "v7tdx"
SITE = "Bdark"
TOPOLOGY = REPO / "tc_tandem_tdx_dark2/solvated_protonated.pdb"
TRAJECTORY = REPO / "tc_tandem_tdx_dark2/tandem_tdx_dark2_1ns.dcd"
CACHE = ROOT / "cache/embedding_charges_v7tdx.npz"
CR2_PRMTOP = REPO / "anionic_build/monomer_solv.prmtop"


def prepare_dark(frame: int, frame_dir: Path) -> dict:
    site_dir = frame_dir / f"site_{SITE}"
    marker = frame_dir / "preparation_darksite.json"
    if marker.exists() and (site_dir / "geometry.xyz").exists():
        return json.loads(marker.read_text())
    staging = frame_dir / "_darkprep"
    if staging.exists():
        subprocess.run(["rm", "-rf", str(staging)], check=True)
    R.run_checked([
        sys.executable, str(ROOT / "scripts/prepare_production_frame.py"),
        "--topology", str(TOPOLOGY), "--trajectory", str(TRAJECTORY),
        "--frame-index", str(frame), "--sites", SITE,
        "--output-dir", str(staging), "--embedding-cache", str(CACHE),
        "--full-system-embedding", "--link-only-exclusion",
        "--retain-partner-cr2-charges", "--amber-cr2-prmtop", str(CR2_PRMTOP),
    ])
    summary = json.loads((staging / "preparation.json").read_text())
    site = summary["sites"][SITE]
    if site["qm_net_charge_e"] != 0:
        raise RuntimeError(f"frame {frame}: dark region is not neutral")
    if site["partner_cr2_residues_restored"] != 1:
        raise RuntimeError(
            f"frame {frame}: the live CR2 was not restored into the field "
            f"({site['partner_cr2_residues_restored']} restored)")
    if site_dir.exists():
        subprocess.run(["rm", "-rf", str(site_dir)], check=True)
    (staging / f"site_{SITE}").rename(site_dir)
    marker.write_text(json.dumps(summary, indent=2) + "\n")
    subprocess.run(["rm", "-rf", str(staging)], check=True)
    return summary


def main() -> None:
    ap = __import__("argparse").ArgumentParser()
    ap.add_argument("--manifest-suffix", default="")
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--roots", type=int, default=7)
    ap.add_argument("--profile", default="svp")
    ap.add_argument("--site-dir-suffix", default="",
                    help="write into site_Bdark<suffix>/ instead of site_Bdark/, "
                         "so a wider root window can be run without disturbing "
                         "the completed 7-root sweep")
    ap.add_argument("--every", type=int, default=1,
                    help="stride over the manifest, for diagnostic subsets")
    args = ap.parse_args()

    path = ROOT / "manifests" / f"{WINDOW}dark{args.manifest_suffix}_frames.csv"
    columns, rows = R.load_manifest(path)
    for n, row in enumerate(rows):
        if n % args.every:
            continue
        if row["convergence_status"] == "complete":
            continue
        frame = int(row["frame_index"])
        frame_dir = ROOT / "frames" / WINDOW / f"frame_{frame:04d}"
        site_a_summary = frame_dir / "frame_summary.json"
        if not site_a_summary.exists():
            raise RuntimeError(f"frame {frame}: site A is not complete")
        print(f"[start] {WINDOW} dark frame={frame} time={row['time_ps']} ps", flush=True)
        t0 = time.time()
        try:
            prep = prepare_dark(frame, frame_dir)
            if args.site_dir_suffix:
                alt = frame_dir / f"site_{SITE}{args.site_dir_suffix}"
                alt.mkdir(exist_ok=True)
                for name in ("geometry.xyz", "mm_charges.dat"):
                    target = alt / name
                    if not target.exists():
                        import shutil as _sh
                        _sh.copy2(frame_dir / f"site_{SITE}" / name, target)
            site_b = R.process_site(
                frame_dir / f"site_{SITE}{args.site_dir_suffix}", WINDOW,
                args.gpu, None, args.roots,
                args.profile, R.REFERENCE, qm_charge=0, selection="brightest")
            a = json.loads(site_a_summary.read_text())["site_A"]
            b = site_b["root"]
            # The suffix MUST reach this path: without it a wider-root rerun
            # silently reuses the narrow-root TDC and reports the old J.
            tdc_path = frame_dir / f"tdc_dark{args.site_dir_suffix}/summary.json"
            if not tdc_path.exists():
                tdc_path.parent.mkdir(parents=True, exist_ok=True)
                R.run_checked([
                    sys.executable, str(ROOT / "scripts/mutual_tdc.py"),
                    "--site-a", str(frame_dir / f"site_A/root{a['root']}_density.npz"),
                    "--site-b", str(site_b["compact_path"]),
                    "--out", str(tdc_path), "--epsilon", "1.77",
                    # This frame's OWN cell, not run_ensemble.BOX_A, which is a
                    # hard-coded 90.7522 A matching neither trajectory (v5 runs
                    # at 90.6101, TDX at 90.5625). At a 26 A separation in a
                    # 90 A box the nearest image is unambiguous either way, so
                    # this is provenance rather than a numerical correction.
                    "--box-lengths", *[str(v) for v in prep["box_A"]],
                ])
            tdc = json.loads(tdc_path.read_text())
            delta = a["energy_cm-1"] - b["energy_cm-1"]
            row.update({
                "site_B_energy_eV": f"{b['energy_eV']:.10f}",
                "site_B_energy_cm-1": f"{b['energy_cm-1']:.6f}",
                "site_B_wavelength_nm": f"{b['wavelength_nm']:.6f}",
                "site_B_oscillator_strength": f"{b['oscillator_strength']:.8f}",
                "site_B_selected_root": str(b["root"]),
                "site_B_assignment_score": f"{site_b['assignment']['selected_assignment_score']:.8f}",
                "J_vacuum_cm-1": f"{tdc['J_vacuum_cm-1']:.8f}",
                "J_screened_cm-1": f"{tdc['J_screened_cm-1']:.8f}",
                "J_pda_vacuum_cm-1": f"{tdc['J_pda_vacuum_cm-1']:.8f}",
                "J_pda_screened_cm-1": f"{tdc['J_pda_screened_cm-1']:.8f}",
                "separation_A": f"{tdc['separation_A']:.8f}",
                "dipole_angle_deg": f"{tdc['dipole_angle_deg']:.8f}",
                "grid_normalization_error_B": f"{site_b['grid_normalization_error']:.8e}",
                "convergence_status": "complete",
                "failure_class": "",
                "runtime_seconds": f"{time.time() - t0:.6f}",
                "gpu": args.gpu,
                "notes": (
                    "dVenus-TDX barrel 2 (Gly-Gly + Tyr473, 29 atoms, neutral); "
                    "brightest-root selection, NOT density-overlap tracking; "
                    "live CR2 restored at -1 e in the field"
                ),
            })
            print(f"[done] dark frame={frame} root={b['root']} "
                  f"EB={b['energy_eV']:.4f} eV f={b['oscillator_strength']:.5f} "
                  f"Delta={delta:+.0f} cm-1 J={tdc['J_screened_cm-1']:+.4f} cm-1",
                  flush=True)
        except Exception as error:                      # noqa: BLE001
            row["convergence_status"] = "failed"
            row["failure_class"] = type(error).__name__
            row["notes"] = repr(error)[:400]
            print(f"[failed] dark frame={frame}: {error!r}", flush=True)
        R.save_manifest(path, columns, rows)


if __name__ == "__main__":
    main()

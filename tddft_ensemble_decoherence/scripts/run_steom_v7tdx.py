#!/usr/bin/env python3
"""One DLPNO-STEOM-CCSD reference for the dVenus-TDX campaign.

The v5 campaign ran eight (four frames x two sites) to measure the STEOM-minus-
TDDFT shift with an error bar. TDX gets one, at the first sampled frame, because
the shift is not what this campaign is for: what the TDX ensemble measures is the
band SHAPE -- the static disorder and the reorganisation energy -- and shape is
independent of the origin. The single job is a spot check that the TDX
chromophore behaves like the TD one at the correlated level, not a calibration
with a quoted uncertainty, and it must not be reported as one.

Everything except the construct is the v5 recipe: same def2-SVP STEOM input,
same %mdci block (TCutPNOSingles/NDav/MaxIter are the documented recipe for
getting parallel EOM to converge at all), same geometry and point-charge files
that TeraChem was given for the same frame.

Differences forced by the construct, all of them deliberate:
  * one site, not two -- barrel 2 has no chromophore;
  * no --site-b-offset and no --retain-partner-cr2-charges -- there is no
    partner CR2 whose RESP charge could be restored, and removing that -1 e from
    the field is the physical content of TDX;
  * its own embedding cache, since the topology has 508 residues not 511.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_steom_v5 as v5                                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
TOPOLOGY = REPO / "tc_tandem_tdx_dark2/solvated_protonated.pdb"
TRAJECTORY = REPO / "tc_tandem_tdx_dark2/tandem_tdx_dark2_1ns.dcd"
CACHE = ROOT / "cache/embedding_charges_v7tdx.npz"
OUT = ROOT / "steom_v7tdx"
FRAME = 100          # the first sampled frame, DCD index; 2 ps into production
SITE = "A"


def prepare(frame: int) -> Path:
    d = OUT / f"frame_{frame:04d}"
    if (d / "preparation.json").exists():
        return d
    if d.exists():
        shutil.rmtree(d)
    subprocess.run([
        sys.executable, str(ROOT / "scripts/prepare_production_frame.py"),
        "--topology", str(TOPOLOGY),
        "--trajectory", str(TRAJECTORY),
        "--frame-index", str(frame),
        "--sites", SITE,
        "--output-dir", str(d),
        "--embedding-cache", str(CACHE),
        "--full-system-embedding", "--link-only-exclusion",
    ], check=True, capture_output=True, text=True)
    return d


def main() -> None:
    v5.OUT = OUT
    OUT.mkdir(parents=True, exist_ok=True)
    v5.log(f"===== TDX STEOM: frame {FRAME} site {SITE} (one reference) =====")
    prepare(FRAME)
    # v5.write_inputs does two things that must not be reimplemented here: it
    # converts TeraChem's mm_charges.dat (N / comment / q x y z) into ORCA's .pc
    # (N / q x y z), and it writes steom.inp from the shared template. Writing
    # the input without the conversion is what failed the first attempt --
    # ORCA aborted in Startup with CANNOT OPEN FILE on the missing field.pc.
    v5.write_inputs(FRAME, SITE, nprocs=5, maxcore=2000)
    result = v5.run_one((FRAME, SITE), nprocs=5, maxcore=2000, hours=24.0,
                        min_free_gb=12, min_disk_gb=150)
    (OUT / "steom_v7tdx_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    v5.log(f"wrote {OUT / 'steom_v7tdx_summary.json'}: {result}")


if __name__ == "__main__":
    main()

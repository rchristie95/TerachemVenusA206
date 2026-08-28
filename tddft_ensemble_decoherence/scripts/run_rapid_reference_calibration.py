#!/usr/bin/env python3
"""Calibrate the rapid TDDFT profile to the retained same-geometry STEOM state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from run_ensemble import REFERENCE, ROOT, process_site, sha256


REPO = ROOT.parent
STEOM_REFERENCE_EV = 2.366633


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu", default="0")
    parser.add_argument("--roots", type=int, choices=(4, 7), default=4)
    parser.add_argument(
        "--profile", choices=("rapid", "rapid_cam"), default="rapid_cam"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "rapid_100ns50",
        help="Campaign result directory (default: rapid_100ns50)",
    )
    args = parser.parse_args()

    campaign_dir = args.output_dir.resolve()
    outdir = campaign_dir / "reference_calibration"
    outdir.mkdir(parents=True, exist_ok=True)
    for name in ("geometry.xyz", "mm_charges.dat"):
        source = REPO / "tc_tddft_44" / name
        destination = outdir / name
        if destination.exists() and sha256(destination) != sha256(source):
            raise RuntimeError(f"Refusing to overwrite mismatched retained file: {destination}")
        if not destination.exists():
            shutil.copy2(source, destination)

    site = process_site(
        outdir,
        "long100ns50",
        args.gpu,
        None,
        args.roots,
        args.profile,
        REFERENCE,
    )
    selected = site["root"]
    delta = STEOM_REFERENCE_EV - float(selected["energy_eV"])
    output_text = (outdir / "td_density.out").read_text(errors="replace")
    version = None
    for line in output_text.splitlines():
        if "TeraChem" in line and ("version" in line.lower() or "v1." in line.lower()):
            version = line.strip()
            break
    gpu = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=index,name,driver_version",
            "--format=csv,noheader",
        ],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    compact = Path(site["compact_path"])
    result = {
        "status": "complete",
        "profile": args.profile,
        "method": "CAM-B3LYP" if args.profile == "rapid_cam" else "B3LYP",
        "basis": "6-31G*",
        "purpose": "rapid prototyping; not a replacement for the validated wB97X-D3/6-311G** benchmark",
        "roots": args.roots,
        "selected_root": int(selected["root"]),
        "selection_basis": site["assignment"]["selection_basis"],
        "assignment_score_against_retained_reference": site["assignment"][
            "selected_assignment_score"
        ],
        "tddft_reference_energy_eV": float(selected["energy_eV"]),
        "tddft_reference_wavelength_nm": float(selected["wavelength_nm"]),
        "tddft_reference_oscillator_strength": float(selected["oscillator_strength"]),
        "steom_reference_energy_eV": STEOM_REFERENCE_EV,
        "steom_reference_wavelength_nm": 523.90,
        "steom_reference_oscillator_strength": 0.885045,
        "steom_constant_shift_eV": delta,
        "constant_shift_effect": (
            "Changes absolute site energies/wavelengths only; cancels exactly from "
            "the A-B gap and cannot change a gap-fluctuation decoherence time."
        ),
        "state_tracking_reference": str(compact.resolve()),
        "state_tracking_reference_sha256": sha256(compact),
        "geometry": str((outdir / "geometry.xyz").resolve()),
        "geometry_sha256": sha256(outdir / "geometry.xyz"),
        "point_charges": str((outdir / "mm_charges.dat").resolve()),
        "point_charges_sha256": sha256(outdir / "mm_charges.dat"),
        "input_template": str(
            (
                ROOT
                / (
                    f"input_templates/tddft_density_rapid_cam_{args.roots}roots.in"
                    if args.profile == "rapid_cam"
                    else f"input_templates/tddft_density_rapid_{args.roots}roots.in"
                )
            ).resolve()
        ),
        "input_template_sha256": sha256(
            ROOT
            / (
                f"input_templates/tddft_density_rapid_cam_{args.roots}roots.in"
                if args.profile == "rapid_cam"
                else f"input_templates/tddft_density_rapid_{args.roots}roots.in"
            )
        ),
        "terachem_version_line": version,
        "gpu": gpu,
        "runtime_seconds": site["runtime_seconds"],
    }
    target = campaign_dir / "calibration.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Parse compact, machine-readable diagnostics from a TeraChem TDDFT output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


FINAL_ROW = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+"
    r"(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+)\s+->\s+(\d+)\s*:\s*([A-Z])\s+->\s+([A-Z])"
)
DIPOLE_ROW = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s*$"
)
CI_ROW = re.compile(
    r"^\s*(\d+)\s+->\s+(\d+)\s*:\s*([A-Z])\s+->\s+([A-Z])\s*:\s*(-?\d+\.\d+)"
)


def parse(path: Path) -> dict:
    text = path.read_text(errors="replace")
    completed = "Job finished:" in text
    requested = re.findall(r"^\s*cisnumstates\s+(\d+)\s*$", text, flags=re.MULTILINE)
    requested_roots = int(requested[-1]) if requested else None
    final_energy = re.findall(r"FINAL ENERGY:\s*(-?\d+\.\d+)\s+a\.u\.", text)
    mm_count = re.findall(r"Number of MM point charges:\s*(\d+)", text)
    gpu = re.findall(r"Device\s+\d+:\s+([^,]+),", text)

    roots: dict[int, dict] = {}
    final_section = text.split("Final Excited State Results:")[-1]
    for line in final_section.splitlines():
        match = FINAL_ROW.match(line)
        if not match:
            continue
        root = int(match.group(1))
        roots[root] = {
            "root": root,
            "total_energy_au": float(match.group(2)),
            "energy_eV": float(match.group(3)),
            "energy_cm-1": float(match.group(3)) * 8065.544005,
            "wavelength_nm": 1239.841984 / float(match.group(3)),
            "oscillator_strength": float(match.group(4)),
            "s2": float(match.group(5)),
            "largest_ci_coefficient": float(match.group(6)),
            "largest_excitation": {
                "occupied": int(match.group(7)),
                "virtual": int(match.group(8)),
                "occupied_type": match.group(9),
                "virtual_type": match.group(10),
            },
            "transition_dipole_au": None,
            "transition_dipole_magnitude_au": None,
            "ci_coefficients": [],
        }

    marker = "Transition dipole moments:\n"
    if marker in text:
        dipole_section = text.split(marker, 1)[1]
        for line in dipole_section.splitlines():
            match = DIPOLE_ROW.match(line)
            if not match:
                continue
            root = int(match.group(1))
            if root in roots:
                roots[root]["transition_dipole_au"] = [
                    float(match.group(2)),
                    float(match.group(3)),
                    float(match.group(4)),
                ]
                roots[root]["transition_dipole_magnitude_au"] = float(match.group(5))
            if requested_roots is not None and root == requested_roots:
                break

    for root, item in roots.items():
        block_match = re.search(
            rf"Root {root}: Largest CI coefficients:\n(.*?)(?:\n\s*\n|\Z)",
            text,
            flags=re.DOTALL,
        )
        if not block_match:
            continue
        for line in block_match.group(1).splitlines():
            match = CI_ROW.match(line)
            if match:
                item["ci_coefficients"].append(
                    {
                        "occupied": int(match.group(1)),
                        "virtual": int(match.group(2)),
                        "occupied_type": match.group(3),
                        "virtual_type": match.group(4),
                        "coefficient": float(match.group(5)),
                    }
                )

    error_patterns = {
        "scf_failure": r"SCF.*(?:failed|not converged)|failed to converge SCF",
        "excited_state_failure": r"Davidson.*(?:failed|not converged)|CIS.*not converged",
        "gpu_memory_failure": r"out of memory|CUDA_ERROR_OUT_OF_MEMORY",
        "numerical_failure": r"\bnan\b|numerical failure",
    }
    detected_errors = [
        label for label, pattern in error_patterns.items() if re.search(pattern, text, re.IGNORECASE)
    ]
    converged = (
        completed
        and requested_roots is not None
        and len(roots) == requested_roots
        and all(item["transition_dipole_au"] is not None for item in roots.values())
        and not detected_errors
    )
    return {
        "status": "converged" if converged else "failed_or_incomplete",
        "completed": completed,
        "requested_roots": requested_roots,
        "parsed_roots": len(roots),
        "scf_final_energy_au": float(final_energy[-1]) if final_energy else None,
        "mm_point_charges": int(mm_count[-1]) if mm_count else None,
        "gpu_model": gpu[-1].strip() if gpu else None,
        "detected_errors": detected_errors,
        "roots": [roots[index] for index in sorted(roots)],
        "source_output": str(path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = parse(args.output)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.json:
        args.json.write_text(rendered)
    print(rendered, end="")
    if result["status"] != "converged":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

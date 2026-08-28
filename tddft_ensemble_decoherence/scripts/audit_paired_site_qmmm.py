#!/usr/bin/env python3
"""Audit the paired 50-frame, 100-site TDDFT QM/MM campaign.

Each retained dimer/MM snapshot must contain two independent TeraChem jobs:
one 44-atom chromophore-plus-stacked-Tyr QM region for site A and one for site
B.  This auditor never launches TeraChem.  It writes a site-level restart
manifest that distinguishes the 100 retained calculations from archived retry
attempts and classifies every missing or inconsistent site explicitly.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WINDOW = "short100fs50_fixedmm"
EXPECTED_SITES = {
    "A": {"cr2_key": ["A", "66"], "stacked_tyr_key": ["A", "202"]},
    "B": {"cr2_key": ["A", "328"], "stacked_tyr_key": ["A", "464"]},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Missing CSV header: {path}")
        return reader.fieldnames, list(reader)


def selected_site(summary: dict, site: str) -> dict:
    if "global_assignment" in summary:
        return summary["global_assignment"][f"site_{site}"]
    return summary[f"site_{site}"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", default=DEFAULT_WINDOW)
    parser.add_argument("--expected-frames", type=int, default=50)
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=ROOT / "rapid_100fs50_fixedmm",
    )
    args = parser.parse_args()

    result_dir = args.result_dir.resolve()
    frame_root = ROOT / "frames" / args.window
    source_manifest = result_dir / "frame_manifest.csv"
    if not source_manifest.exists():
        source_manifest = ROOT / "manifests" / f"{args.window}_frames.csv"
    _, frame_rows = load_rows(source_manifest)
    if len(frame_rows) != args.expected_frames:
        raise RuntimeError(
            f"Expected {args.expected_frames} frame rows, found {len(frame_rows)}"
        )

    site_rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    paired_complete_frames = 0
    topology_paths: set[str] = set()
    trajectory_paths: set[str] = set()

    for frame_row in frame_rows:
        frame = int(frame_row["frame_index"])
        frame_dir = frame_root / f"frame_{frame:04d}"
        preparation_path = frame_dir / "preparation.json"
        summary_path = frame_dir / "frame_summary_global.json"
        if not summary_path.exists():
            summary_path = frame_dir / "frame_summary.json"
        preparation = (
            json.loads(preparation_path.read_text())
            if preparation_path.exists()
            else {}
        )
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        if preparation:
            topology_paths.add(preparation["topology"])
            trajectory_paths.add(preparation["trajectory"])

        pair_complete = True
        for site in ("A", "B"):
            problems: list[str] = []
            site_dir = frame_dir / f"site_{site}"
            output_path = site_dir / "td_density.out"
            density_summary_path = site_dir / "density_summary.json"
            method_path = site_dir / "method_profile.json"
            runtime_path = site_dir / "runtime_seconds.txt"
            prep_site = preparation.get("sites", {}).get(site, {})

            if frame_row.get("convergence_status") != "complete":
                problems.append(
                    f"frame_manifest_{frame_row.get('convergence_status') or 'unclassified'}"
                )
            if not preparation_path.exists():
                problems.append("missing_preparation")
            if prep_site.get("qm_atoms") != 44:
                problems.append("qm_atom_count_not_44")
            if prep_site.get("physical_qm_atoms") != 41:
                problems.append("physical_qm_atom_count_not_41")
            if prep_site.get("cr2_key") != EXPECTED_SITES[site]["cr2_key"]:
                problems.append("wrong_chromophore_mapping")
            if (
                prep_site.get("stacked_tyr_key")
                != EXPECTED_SITES[site]["stacked_tyr_key"]
            ):
                problems.append("wrong_stacked_tyr_mapping")
            caps = prep_site.get("link_caps", [])
            if len(caps) != 3:
                problems.append("link_cap_count_not_3")
            elif any(
                abs(float(cap.get("cap_direction_dot_parent_to_partner", 0.0)) - 1.0)
                > 1.0e-10
                for cap in caps
            ):
                problems.append("invalid_link_cap_direction")
            if prep_site.get("excluded_link_atom_indices_one_based") != [42, 43, 44]:
                problems.append("wrong_density_link_atom_exclusion")

            method = json.loads(method_path.read_text()) if method_path.exists() else {}
            density_summary = (
                json.loads(density_summary_path.read_text())
                if density_summary_path.exists()
                else {}
            )
            if not method_path.exists():
                problems.append("missing_method_profile")
            if not output_path.exists():
                problems.append("missing_tddft_output")
            elif "Job finished:" not in output_path.read_text(errors="replace"):
                problems.append("unfinished_tddft_output")
            if density_summary.get("status") != "converged":
                problems.append("unconverged_tddft_summary")
            if not runtime_path.exists():
                problems.append("missing_runtime")

            chosen: dict = {}
            if summary:
                try:
                    chosen = selected_site(summary, site)
                except KeyError:
                    problems.append("missing_selected_state")
            else:
                problems.append("missing_frame_summary")
            compact_path = (
                Path(chosen["compact_path"]) if chosen.get("compact_path") else None
            )
            if compact_path is None or not compact_path.exists():
                problems.append("missing_selected_compact_density")

            status = "complete" if not problems else "failed"
            pair_complete &= status == "complete"
            entry: dict[str, object] = {
                "window": args.window,
                "frame_index": frame,
                "time_fs": frame_row["time_fs"],
                "coordinate_sha256": frame_row["coordinate_sha256"],
                "site": site,
                "qm_partition": (
                    "29-atom anionic CR2 + 12-atom stacked Tyr phenol + "
                    "3 link H; one site QM at a time in the dimer/MM snapshot"
                ),
                "qm_atoms": prep_site.get("qm_atoms", ""),
                "physical_qm_atoms": prep_site.get("physical_qm_atoms", ""),
                "link_atoms": len(caps),
                "excluded_link_atom_indices_one_based": json.dumps(
                    prep_site.get("excluded_link_atom_indices_one_based", [])
                ),
                "cr2_key": ":".join(prep_site.get("cr2_key", [])),
                "stacked_tyr_key": ":".join(prep_site.get("stacked_tyr_key", [])),
                "mm_point_charges": prep_site.get("mm_point_charges", ""),
                "mm_net_charge_e": prep_site.get("mm_net_charge_e", ""),
                "method_profile": method.get("profile", ""),
                "roots_calculated": method.get("roots", ""),
                "selected_root": chosen.get("root", ""),
                "selected_energy_eV": chosen.get("energy_eV", ""),
                "selected_oscillator_strength": chosen.get(
                    "oscillator_strength", ""
                ),
                "assignment_score": chosen.get("reference_spatial_overlap", ""),
                "runtime_seconds": (
                    float(runtime_path.read_text()) if runtime_path.exists() else ""
                ),
                "gpu": frame_row.get("gpu", ""),
                "tddft_output": str(output_path.resolve()),
                "tddft_output_sha256": (
                    sha256(output_path) if output_path.exists() else ""
                ),
                "selected_compact_density": (
                    str(compact_path.resolve()) if compact_path is not None else ""
                ),
                "selected_compact_density_sha256": (
                    sha256(compact_path)
                    if compact_path is not None and compact_path.exists()
                    else ""
                ),
                "status": status,
                "failure_class": ";".join(problems),
            }
            site_rows.append(entry)
            if problems:
                failures.append(
                    {"frame_index": frame, "site": site, "problems": problems}
                )
        paired_complete_frames += int(pair_complete)

    columns = list(site_rows[0])
    site_manifest = result_dir / "paired_site_qmmm_manifest.csv"
    temporary_manifest = site_manifest.with_suffix(".csv.tmp")
    with temporary_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(site_rows)
    temporary_manifest.replace(site_manifest)

    complete_by_site = {
        site: sum(row["status"] == "complete" and row["site"] == site for row in site_rows)
        for site in ("A", "B")
    }
    profiles = Counter(
        f"{row['method_profile']}:{row['roots_calculated']}roots"
        for row in site_rows
        if row["status"] == "complete"
    )
    total_runtime = sum(
        float(row["runtime_seconds"])
        for row in site_rows
        if row["runtime_seconds"] != ""
    )
    audit = {
        "status": "complete" if not failures else "failed",
        "policy": (
            "Two separate chromophore-plus-stacked-Tyr QM/MM TDDFTs are "
            "required for every retained dimer/MM frame. Archived retry "
            "directories are preserved but excluded from the retained count. "
            "This audit never launches or replaces a calculation."
        ),
        "window": args.window,
        "source_frame_manifest": str(source_manifest.resolve()),
        "site_manifest": str(site_manifest.resolve()),
        "expected_frames": args.expected_frames,
        "expected_site_calculations": 2 * args.expected_frames,
        "retained_site_rows": len(site_rows),
        "complete_site_calculations": sum(
            row["status"] == "complete" for row in site_rows
        ),
        "complete_by_site": complete_by_site,
        "paired_complete_frames": paired_complete_frames,
        "failed_site_calculations": len(failures),
        "failures": failures,
        "method_profiles": dict(sorted(profiles.items())),
        "summed_terachem_runtime_seconds": total_runtime,
        "topologies": sorted(topology_paths),
        "trajectories": sorted(trajectory_paths),
        "qm_region": {
            "per_site_atoms": 44,
            "physical_atoms": 41,
            "composition": (
                "29-atom anionic CR2 chromophore + 12-atom stacked Tyr phenol "
                "+ 3 link hydrogens"
            ),
            "link_atom_density_policy": (
                "Atoms 42-44 are retained in QM electronic-structure jobs and "
                "excluded from compact transition densities and TDC."
            ),
        },
        "environment": (
            "Both site jobs use the same instantaneous full-dimer NVT frame "
            "with the validated fixed frame-0 local MM point-charge identities; "
            "only one chromophore-plus-Tyr region is QM in each job."
        ),
        "restart_semantics": (
            "The production worker skips a site whose retained td_density.out "
            "contains 'Job finished:' and resumes only a missing/incomplete site; "
            "a complete frame manifest prevents duplicate paired jobs."
        ),
    }
    audit_path = result_dir / "paired_site_qmmm_audit.json"
    temporary_audit = audit_path.with_suffix(".json.tmp")
    temporary_audit.write_text(json.dumps(audit, indent=2) + "\n")
    temporary_audit.replace(audit_path)
    print(json.dumps(audit, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())

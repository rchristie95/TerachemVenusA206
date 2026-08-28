#!/usr/bin/env python3
"""Write compact, machine-readable provenance for the completed 100-fs run."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]
RESULT = REPO / "tddft_ensemble_decoherence/wb97xd3_100fs100"
FRAME_ROOT = REPO / "tddft_ensemble_decoherence/frames/short100fs100_wb97xd3"
TRAJECTORY = (
    REPO
    / "tc_tandem_nvt_100fs_100frames_20260724/tandem_nvt_100fs_100frames.dcd"
)
TOPOLOGY = (
    REPO / "tc_tandem_nvt_100fs_100frames_20260724/solvated_protonated.pdb"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command(*args: str) -> str:
    return subprocess.check_output(args, cwd=REPO, text=True).strip()


def output_record(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main() -> None:
    analysis = json.loads((RESULT / "summary.json").read_text())
    tracking = json.loads((RESULT / "global_state_tracking.json").read_text())
    site_audit = json.loads((RESULT / "paired_site_qmmm_audit.json").read_text())
    benchmark = json.loads((RESULT / "adaptive_parallel_benchmark.json").read_text())
    calibration = json.loads((RESULT / "calibration.json").read_text())
    movie = json.loads(
        (
            RESULT
            / "tddft_density_wb97xd3_100fs_50s_fixed_camera_isosurfaces.provenance.json"
        ).read_text()
    )

    frame_summaries = [
        json.loads(
            (
                FRAME_ROOT / f"frame_{frame:04d}/frame_summary_global.json"
            ).read_text()
        )
        for frame in range(100)
    ]
    tdc = [row["tdc"] for row in frame_summaries]
    symmetry_relative = np.asarray(
        [abs(float(row["symmetry_relative_error"])) for row in tdc]
    )
    cpu_gpu_relative = np.asarray(
        [abs(float(row["cpu_gpu_relative_error"])) for row in tdc]
    )

    git_status = command("git", "status", "--short", "--branch").splitlines()
    gpu = command(
        "nvidia-smi",
        "--query-gpu=name,driver_version,memory.total",
        "--format=csv,noheader,nounits",
    )
    tc_output = (
        FRAME_ROOT / "frame_0000/site_A/td_density.out"
    ).read_text(errors="replace")
    tc_version = next(
        line.strip(" |*")
        for line in tc_output.splitlines()
        if "TeraChem v" in line
    )

    inputs = {
        name: output_record(REPO / relative)
        for name, relative in {
            "tddft_template": "tddft_ensemble_decoherence/input_templates/tddft_density_7roots.in",
            "state_tracker": "tddft_ensemble_decoherence/scripts/global_state_tracking_100fs.py",
            "analysis_script": "tddft_ensemble_decoherence/scripts/summarize_100fs50.py",
            "fixed_embedding_selection": "tddft_ensemble_decoherence/wb97xd3_100fs100/fixed_embedding_selection.json",
            "source_checkpoint": "tc_tandem_nvt_rerun_20260722_retry2/nvt.chk",
            "final_checkpoint": "tc_tandem_nvt_100fs_100frames_20260724/nvt.chk",
        }.items()
    }
    outputs = {
        name: output_record(RESULT / relative)
        for name, relative in {
            "frame_manifest": "frame_manifest.csv",
            "paired_site_manifest": "paired_site_qmmm_manifest.csv",
            "global_state_tracking": "global_state_tracking.json",
            "summary": "summary.json",
            "numerical_report": "NUMERICAL_REPORT.md",
            "histograms": "energy_coupling_histograms.png",
            "coherence_plot": "correlation_and_coherence.png",
            "movie": "tddft_density_wb97xd3_100fs_50s_fixed_camera_isosurfaces.mp4",
            "movie_montage": "video_isosurfaces_first_middle_last.png",
            "legacy_sphere_movie": "tddft_density_wb97xd3_100fs_50s_fixed_camera_visible.mp4",
        }.items()
    }

    provenance = {
        "status": "complete",
        "scope": (
            "100 consecutive full-dimer NVT frames spanning 1--100 fs; "
            "two independently embedded chromophore-plus-Tyr site TDDFTs per frame"
        ),
        "repository": {
            "path": str(REPO),
            "branch": command("git", "branch", "--show-current"),
            "commit": command("git", "rev-parse", "HEAD"),
            "status_short_branch": git_status,
            "note": "Pre-existing and campaign worktree changes were preserved; no commit or push was made.",
        },
        "software": {
            "terachem": tc_version,
            "orca_retained_steom": "6.1.1 RELEASE",
            "cuda_gpu_driver_memory_mib": gpu,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "trajectory": {
            "path": str(TRAJECTORY.resolve()),
            "sha256": sha256(TRAJECTORY),
            "frames": 100,
            "time_fs": [1.0, 100.0],
            "sampling_interval_fs": 1.0,
            "atoms": 69456,
            "periodic_box_A": [90.7521749130615] * 3,
            "topology": str(TOPOLOGY.resolve()),
            "topology_sha256": sha256(TOPOLOGY),
            "interpolated_frames": 0,
        },
        "electronic_structure": {
            "method": "wB97X-D3",
            "basis": "6-311G**",
            "charge": -1,
            "multiplicity": 1,
            "roots_computed": 7,
            "target": "bright anionic chromophore pi-pi*, dominant 93 -> 95 configuration",
            "qm_atoms_per_site": 44,
            "physical_density_atoms": 41,
            "link_atoms_excluded_one_based": [42, 43, 44],
            "embedding": (
                "validated fixed frame-0 local MM point-charge identities in "
                "the instantaneous full-dimer NVT geometry"
            ),
        },
        "campaign": {
            "frames_complete": site_audit["paired_complete_frames"],
            "site_jobs_complete": site_audit["complete_site_calculations"],
            "site_jobs_failed": site_audit["failed_site_calculations"],
            "failure_fraction": 0.0,
            "retry_count": sum(
                int(row.get("retry_count", 0))
                for row in benchmark.get("retry_audit", [])
                if isinstance(row, dict)
            ),
            "selected_concurrent_workers": benchmark["selected_workers"],
            "summed_terachem_runtime_seconds": site_audit[
                "summed_terachem_runtime_seconds"
            ],
        },
        "state_tracking": {
            "required_configuration": "93 -> 95",
            "site_A_root_counts": analysis["selected_root_counts"]["site_A"],
            "site_B_root_counts": analysis["selected_root_counts"]["site_B"],
            "changes_from_static_assignment": {
                site: len(tracking["sites"][site]["changed_frames"])
                for site in ("A", "B")
            },
            "minimum_adjacent_full_density_overlap": {
                site: tracking["sites"][site][
                    "minimum_adjacent_full_spatial_overlap"
                ]
                for site in ("A", "B")
            },
            "tracking_sensitivity_consensus": "12/12 paths at every frame for both sites",
            "exact_all_root_density_edges": tracking[
                "exact_spatial_transition_edges"
            ],
        },
        "tdc": {
            "dielectric_convention": "J_screened = J_vacuum / epsilon_opt",
            "epsilon_opt": 1.77,
            "mean_vacuum_cm-1": analysis["distributions"]["J_vacuum_cm-1"][
                "mean"
            ],
            "mean_screened_cm-1": analysis["distributions"][
                "J_screened_cm-1"
            ]["mean"],
            "mean_pda_screened_cm-1": analysis["distributions"][
                "J_PDA_screened_cm-1"
            ]["mean"],
            "maximum_AB_BA_relative_error": float(symmetry_relative.max()),
            "maximum_cpu_gpu_subset_relative_error": float(
                cpu_gpu_relative.max()
            ),
            "cpu_subset_points_per_frame": 1000,
            "coupling_energy_shift_or_scaling": "none",
        },
        "calibration": calibration,
        "decoherence": {
            "cumulant_T2star_raw_TDDFT_fs": analysis["cumulant"][
                "T2star_raw_TDDFT_fs"
            ],
            "cumulant_T2star_shifted_STEOM_fs": analysis["cumulant"][
                "T2star_constant_shifted_STEOM_fs"
            ],
            "constant_shift_gap_invariance_max_abs_eV": analysis[
                "constant_shift_gap_invariance_max_abs_eV"
            ],
            "bath_positive_correlation_time_fs": analysis["cumulant"][
                "bath_positive_correlation_time_fs"
            ],
            "direct_full_signed_J_fs": analysis[
                "direct_two_site_propagation"
            ]["T2_e_minus_1_full_signed_J_fs"],
            "direct_gap_only_fs": analysis["direct_two_site_propagation"][
                "T2_e_minus_1_gap_only_fs"
            ],
            "stride_sensitivity": analysis["cumulant"]["stride_sensitivity"],
            "segment_sensitivity": analysis["cumulant"][
                "segment_sensitivity"
            ],
            "bootstrap": analysis["cumulant"]["moving_block_bootstrap"],
            "limitation": (
                "The 100-fs record is sufficient for a prototype ultrafast "
                "estimate but not a stationary 1-ps block analysis or slow "
                "inhomogeneous broadening."
            ),
        },
        "movie": movie,
        "inputs": inputs,
        "outputs": outputs,
    }
    (RESULT / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()

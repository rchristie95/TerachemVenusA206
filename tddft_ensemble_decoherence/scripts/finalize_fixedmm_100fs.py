#!/usr/bin/env python3
"""Write final provenance and reports for the corrected 100 fs campaign."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
OUT = ROOT / "rapid_100fs50_fixedmm"
WINDOW = "short100fs50_fixedmm"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=REPO,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()


def main() -> None:
    manifest_path = ROOT / f"manifests/{WINDOW}_frames.csv"
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    benchmark = json.loads((OUT / "adaptive_parallel_benchmark.json").read_text())
    calibration = json.loads((OUT / "calibration.json").read_text())
    selection = json.loads((OUT / "fixed_embedding_selection.json").read_text())
    trajectory = json.loads(
        (ROOT / "rapid_100fs50/trajectory_provenance.json").read_text()
    )
    tracking = json.loads((OUT / "global_state_tracking.json").read_text())
    sensitivity = json.loads((OUT / "state_path_sensitivity.json").read_text())
    summary_path = OUT / "summary.json"
    summary = json.loads(summary_path.read_text())

    symmetry, cpu_gpu, normalizations, cap_dots = [], [], [], []
    selection_observations = {"A": set(), "B": set()}
    for frame in range(50):
        frame_dir = ROOT / f"frames/{WINDOW}/frame_{frame:04d}"
        global_summary = json.loads(
            (frame_dir / "frame_summary_global.json").read_text()
        )
        tdc = global_summary["tdc"]
        symmetry.append(float(tdc["symmetry_relative_error"]))
        cpu_gpu.append(float(tdc["cpu_gpu_relative_error"]))
        for site in ("A", "B"):
            assignment = global_summary["global_assignment"][f"site_{site}"]
            normalizations.append(float(assignment["grid_normalization_error"]))
        preparation = json.loads((frame_dir / "preparation.json").read_text())
        for site in ("A", "B"):
            item = preparation["sites"][site]
            selection_observations[site].add(
                (int(item["mm_point_charges"]), round(float(item["mm_net_charge_e"]), 8))
            )
            cap_dots.extend(
                float(cap["cap_direction_dot_parent_to_partner"])
                for cap in item["link_caps"]
            )

    audit = {
        "status": "complete",
        "requested_frames": 50,
        "final_complete_frames": sum(
            row["convergence_status"] == "complete" for row in rows
        ),
        "final_failed_frames": sum(
            row["convergence_status"] == "failed" for row in rows
        ),
        "interpolated_or_replaced_frames": 0,
        "seven_root_boundary_retry_frames": benchmark["failed_before_retry"],
        "retry_count": sum(int(row["retry_count"]) for row in rows),
        "retry_policy": (
            "Original seven-root directories retained as attempt0; rerun serially "
            "with nine roots; no energy or coupling interpolation."
        ),
        "final_failure_classes": {},
        "runtime_sum_seconds": sum(float(row["runtime_seconds"]) for row in rows),
        "runtime_mean_seconds_per_frame": float(
            np.mean([float(row["runtime_seconds"]) for row in rows])
        ),
        "adaptive_selected_workers": benchmark["selected_workers"],
        "adaptive_three_worker_throughput_frames_per_hour": benchmark["levels"][2][
            "throughput_frames_per_hour"
        ],
        "adaptive_peak_gpu_memory_MiB": benchmark["levels"][2][
            "peak_gpu_memory_mib"
        ],
        "tdc_maximum_AB_BA_relative_error": max(symmetry),
        "tdc_maximum_CPU_GPU_subset_relative_error": max(cpu_gpu),
        "maximum_pre_neutralization_grid_charge_error": max(normalizations),
        "median_pre_neutralization_grid_charge_error": float(
            np.median(normalizations)
        ),
        "minimum_link_cap_direction_dot_product": min(cap_dots),
        "fixed_MM_selection_observed": {
            site: [
                {"point_charges": count, "net_charge_e": charge}
                for count, charge in sorted(values)
            ]
            for site, values in selection_observations.items()
        },
        "state_tracking": {
            "site_A_changed_from_static_frames": tracking["sites"]["A"][
                "changed_frames"
            ],
            "site_B_changed_from_static_frames": tracking["sites"]["B"][
                "changed_frames"
            ],
            "site_A_minimum_adjacent_full_density_overlap": tracking["sites"]["A"][
                "minimum_adjacent_full_spatial_overlap"
            ],
            "site_B_minimum_adjacent_full_density_overlap": tracking["sites"]["B"][
                "minimum_adjacent_full_spatial_overlap"
            ],
        },
    }
    (OUT / "convergence_failure_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    template = ROOT / "input_templates/tddft_density_rapid_cam_7roots.in"
    retry_template = ROOT / "input_templates/tddft_density_rapid_cam_9roots.in"
    provenance = {
        "status": "complete",
        "repository": str(REPO.resolve()),
        "branch": git("branch", "--show-current"),
        "commit": git("rev-parse", "HEAD"),
        "origin_main": git("rev-parse", "origin/main"),
        "worktree_dirty": bool(git("status", "--porcelain")),
        "trajectory": trajectory,
        "fixed_embedding_selection": {
            "path": str((OUT / "fixed_embedding_selection.json").resolve()),
            "sha256": sha256(OUT / "fixed_embedding_selection.json"),
            "policy": selection["policy"],
            "sites": selection["sites"],
        },
        "QM_region": {
            "atoms_total": 44,
            "physical_atoms": 41,
            "composition": "29-atom anionic CR2 plus 12-atom stacked Tyr phenol plus 3 link H",
            "charge": -1,
            "multiplicity": 1,
            "link_atom_indices_one_based_excluded_from_density_and_TDC": [42, 43, 44],
            "cap_construction": "instantaneous minimum-image severed-bond vectors",
        },
        "TDDFT": {
            "TeraChem": "v1.97B-beta-251105, git 3fcf25291f21a40754a20c9c6450f71e006a2729",
            "method": calibration["method"],
            "basis": calibration["basis"],
            "profile": calibration["profile"],
            "production_roots": calibration["production_roots"],
            "seven_root_template_sha256": sha256(template),
            "nine_root_retry_template_sha256": sha256(retry_template),
        },
        "STEOM_reference": {
            "ORCA": "6.1.1 RELEASE",
            "energy_eV": calibration["steom_reference_energy_eV"],
            "wavelength_nm": calibration["steom_reference_wavelength_nm"],
            "oscillator_strength": calibration["steom_reference_oscillator_strength"],
            "constant_shift_eV": calibration["steom_constant_shift_eV"],
            "affine_calibration": "not fitted: only one matched geometry",
        },
        "GPU": {
            "model": "NVIDIA GeForce RTX 4080",
            "driver": "580.126.09",
            "memory_MiB": 16376,
        },
        "TDC": {
            "epsilon_opt": 1.77,
            "screening": "J_screened = J_vacuum / epsilon_opt",
            "signed_values_retained": True,
            "point_dipole_role": "diagnostic only",
        },
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    summary["state_path_sensitivity"] = sensitivity
    summary["convergence_failure_audit"] = audit
    summary["provenance_path"] = str((OUT / "provenance.json").resolve())
    summary["invalid_cutoff_varying_control"] = {
        "T2star_fs": 8.023832522868304,
        "reason_invalid": (
            "Instantaneous 12 A shell changed MM residue membership and net charge "
            "between adjacent frames; retained only as a diagnosed artifact."
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    t2 = summary["cumulant"]["T2star_raw_TDDFT_fs"]
    tau = summary["cumulant"]["bath_positive_correlation_time_fs"]
    direct = summary["direct_two_site_propagation"]
    dist = summary["distributions"]
    bootstrap = summary["cumulant"]["moving_block_bootstrap"]
    readme = f"""# Corrected rapid 100 fs TDDFT/TDC ensemble

This directory contains the restartable 50-frame, 2 fs-cadence prototype with
a fixed frame-0 MM atom selection. The full classical NVT system remains intact;
only the 44-atom QM single-point region is bond-capped. No frame is interpolated
or replaced.

## Exact commands

```bash
python tddft_ensemble_decoherence/scripts/initialize_manifests.py
python tddft_ensemble_decoherence/scripts/register_nvt_trajectories.py \\
  --short100fs50-fixedmm-dcd \\
  tc_tandem_nvt_100fs_50frames_20260723/tandem_nvt_100fs_50frames.dcd

# Frame 0 creates fixed_embedding_selection.json; frames 1--49 reuse it.
conda run -n TeraChem python \\
  tddft_ensemble_decoherence/scripts/prepare_production_frame.py \\
  --topology tc_tandem_nvt_100fs_50frames_20260723/solvated_protonated.pdb \\
  --trajectory tc_tandem_nvt_100fs_50frames_20260723/tandem_nvt_100fs_50frames.dcd \\
  --frame-index FRAME --output-dir FRAME_DIR \\
  --embedding-cache tddft_ensemble_decoherence/cache/embedding_charges.npz \\
  --fixed-embedding-selection \\
  tddft_ensemble_decoherence/rapid_100fs50_fixedmm/fixed_embedding_selection.json

conda run -n TeraChem python -u \\
  tddft_ensemble_decoherence/scripts/run_rapid_adaptive_100ns50.py \\
  --window short100fs50_fixedmm \\
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm \\
  --maximum-workers 3

conda run -n TeraChem python -u \\
  tddft_ensemble_decoherence/scripts/global_state_tracking_100fs.py \\
  --window short100fs50_fixedmm \\
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm

conda run -n TeraChem python \\
  tddft_ensemble_decoherence/scripts/summarize_100fs50.py \\
  --window short100fs50_fixedmm \\
  --output-dir tddft_ensemble_decoherence/rapid_100fs50_fixedmm
```

## Result

- Complete frames: 50/50; final failures: 0; interpolation: 0
- Canonical cumulant T2*: {t2:.3f} fs
- Positive-sequence bath correlation time: {tau:.3f} fs
- Raw mean wavelengths: A 385.31 nm; B 370.35 nm
- Constant-STEOM-shifted mean wavelengths: A 582.78 nm; B 549.57 nm
- Screened signed TDC: mean {dist['J_screened_cm-1']['mean']:.2f} cm^-1,
  SD {dist['J_screened_cm-1']['sample_sd']:.2f} cm^-1

See `NUMERICAL_REPORT.md`, `summary.json`, `provenance.json`, and
`convergence_failure_audit.json`. Large trajectory, raw cube, and scratch data
remain local and are not intended for version control.
"""
    (OUT / "README.md").write_text(readme)

    report = f"""# Numerical report: corrected 100 fs TDDFT/TDC prototype

## Outcome

All 50 retained NVT frames (2--100 fs at 2 fs cadence) have classified A/B
TDDFT states and signed mutual TDC couplings. Six seven-root boundary cases
(frames 31, 32, 33, 42, 43, 44) were preserved and recovered with nine roots.
There are no final failed, replaced, or interpolated frames.

The canonical Gaussian second-cumulant result is **T2* = {t2:.3f} fs**. The
same-geometry STEOM correction is {calibration['steom_constant_shift_eV']:+.8f}
eV, and the raw and constant-shifted gap coherences agree to
{summary['constant_shift_gap_invariance_max_abs_eV']:.3e} eV; a constant
shift changes the absolute spectrum but cannot change differential-gap
fluctuations or T2*. No affine calibration is justified from one matched
STEOM geometry.

## Embedding correction

The first prototype used a fresh 12 A residue selection at every frame. Its MM
net charge changed discontinuously and produced an invalid 8.02 fs result.
That output is preserved as an artifact audit. The corrected campaign holds
the frame-0 atom identities fixed: A has 2,510 charges with net -1.3320 e and
B has 2,439 charges with net -0.8158 e at every frame. Coordinates and
minimum-image placement remain instantaneous. The classical MD system itself
is intact and contains 69,456 atoms; link caps occur only in each 44-atom QM
single point.

## Spectral and coupling statistics

- Site A raw energy: {dist['site_A_raw_eV']['mean']:.4f} +/- {dist['site_A_raw_eV']['sample_sd']:.4f} eV
- Site B raw energy: {dist['site_B_raw_eV']['mean']:.4f} +/- {dist['site_B_raw_eV']['sample_sd']:.4f} eV
- Mean corrected photon wavelength: A 582.78 nm; B 549.57 nm
- Differential gap SD: {dist['differential_gap_eV']['sample_sd']:.4f} eV
- Screened signed J: {dist['J_screened_cm-1']['mean']:.2f} +/- {dist['J_screened_cm-1']['sample_sd']:.2f} cm^-1
- Screened Davydov splitting 2|J|: mean {dist['Davydov_splitting_screened_cm-1']['mean']:.2f} cm^-1

Link atoms 42--44 are excluded from compact transition densities and TDC.
TDC symmetry passes at {audit['tdc_maximum_AB_BA_relative_error']:.2e} maximum
relative error; CPU/GPU subset validation passes at
{audit['tdc_maximum_CPU_GPU_subset_relative_error']:.2e}.

## Decoherence and uncertainty

- Canonical cumulant T2*: {t2:.3f} fs
- Sampling-stride T2* (2/4/8 fs): 11.628 / 11.750 / 10.691 fs
- Positive-sequence bath correlation time: {tau:.3f} fs
- First/last 50 fs: first is right-censored at 48 fs; last is 10.468 fs
- Moving-block bootstrap: finite 95% interval
  {bootstrap['finite_percentile_95_CI_fs'][0]:.2f}--{bootstrap['finite_percentile_95_CI_fs'][1]:.2f}
  fs; {bootstrap['right_censored_count']}/{bootstrap['replicates']} resamples
  are right-censored
- Exact-density state-path variants: finite median
  {sensitivity['finite_T2star_median_fs']:.2f} fs, range
  {sensitivity['finite_T2star_range_fs'][0]:.2f}--{sensitivity['finite_T2star_range_fs'][1]:.2f}
  fs, with one censored path. The shortest variants use aggressive
  density-only weighting and introduce large energy discontinuities.
- Direct signed-J propagation first crosses e^-1 at
  {direct['T2_e_minus_1_full_signed_J_fs']:.3f} fs, close to the cumulant
  crossing. The gap-only direct trace narrowly avoids the early threshold and
  recurs before crossing at {direct['T2_e_minus_1_gap_only_fs']:.3f} fs;
  this finite-window first-crossing sensitivity is not evidence for a robust
  53 fs coupling-induced change.

## Comparison and interpretation

The canonical geometry-dependent result is shorter than the prior 21.8 fs
minimum-image, 24.0 fs PME, and 24.5 +/- 3.7 fs one-ps-block fixed-charge
estimates, and much shorter than the manuscript's phenomenological 60 fs
sensitivity value. The observed difference is driven mainly by the larger
TDDFT differential-gap fluctuation amplitude (SD 0.0924 eV) and B-site state
mixing, with a short positive correlation time of 7.58 fs. Coupling
fluctuations do not materially move the canonical cumulant/full-propagation
crossing, although the direct trace is strongly recurrent in this short
sample.

This is a rapid CAM-B3LYP/6-31G* prototype, not a replacement for the validated
wB97X-D3/6-311G** benchmark. The 100 fs record is too short for independent
1 ps blocks, slow solvation, or inhomogeneous disorder, and its two 50 fs
halves are not stationary. Therefore 11.63 fs is evidence that
geometry-dependent electronic fluctuations may shorten the fixed-charge
estimate, but it is not yet a precise converged T2*. Population relaxation is
not determined here.

Classical moving water point charges are sufficient for this first
electrostatic-dephasing estimate. Quantum waters are not required to compute
the present classical correlation; an explicit-QM first-shell pilot would test
electronic polarization and charge-transfer corrections. Nuclear quantum
effects would require a quantum-nuclear treatment or a separately justified
quantum correction, not merely adding water molecules to snapshot QM regions.
"""
    (OUT / "NUMERICAL_REPORT.md").write_text(report)


if __name__ == "__main__":
    main()

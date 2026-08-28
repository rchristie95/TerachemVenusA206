#!/usr/bin/env python3
"""Globally track the bright diabatic state and refresh affected TDC values.

The first-pass selector compares every root to one static reference.  Near
state crossings that can jump between adiabatic roots.  This script retains
those assignments, then uses a reference-anchored Viterbi path with adjacent
physical-atom transition-charge overlap, transition-dipole orientation,
energy continuity, reference spatial overlap, oscillator strength, and
orbital character.  A broad weight sensitivity grid is recorded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import sys

import numpy as np

from run_ensemble import (
    DELTA_STEOM_EV,
    ROOT,
    load_manifest,
    run_checked,
    save_manifest,
)
from select_state_by_overlap import overlap as full_spatial_overlap


WINDOW = "short100fs50"
MANIFEST = ROOT / "manifests/short100fs50_frames.csv"
OUTDIR = ROOT / "rapid_100fs50"
FRAME_COUNT = 50
BOX_A = [90.7521749130615] * 3
# The first hybrid pass audits the selected path with the full spatial
# transition density.  Use exact all-root transition matrices around every
# interval where that audit was below 0.80, including one-frame shoulders.
# Elsewhere the selected full-density overlap was already high and the atomic
# partition is retained as the inexpensive proxy.
EXACT_SPATIAL_EDGES = {
    "A": set(range(8, 13)) | set(range(25, 32)) | set(range(41, 49)),
    # Site B changes adiabatic root repeatedly even with the fixed MM field;
    # use the full physical transition density at every edge.
    "B": set(range(49)),
}
_EDGE_MATRIX_CACHE: dict[tuple[str, int], np.ndarray] = {}


def load_site(site: str) -> list[list[dict]]:
    frames = []
    for frame in range(FRAME_COUNT):
        site_dir = ROOT / f"frames/{WINDOW}/frame_{frame:04d}/site_{site}"
        assignment = json.loads((site_dir / "state_assignment.json").read_text())
        parsed = json.loads((site_dir / "density_summary.json").read_text())
        overlaps = {
            int(item["root"]): item for item in assignment["candidates"]
        }
        states = []
        for item in parsed["roots"]:
            root = int(item["root"])
            compact_path = site_dir / f"root{root}_density.npz"
            compact = np.load(compact_path)
            states.append(
                {
                    "root": root,
                    "energy_eV": float(item["energy_eV"]),
                    "energy_cm-1": float(item["energy_cm-1"]),
                    "wavelength_nm": float(item["wavelength_nm"]),
                    "oscillator_strength": float(item["oscillator_strength"]),
                    "transition_dipole_au": item["transition_dipole_au"],
                    "transition_dipole_magnitude_au": item[
                        "transition_dipole_magnitude_au"
                    ],
                    "largest_excitation": item["largest_excitation"],
                    "reference_spatial_overlap": float(
                        overlaps[root]["spatial_overlap_absolute"]
                    ),
                    "atom_transition_q": np.asarray(
                        compact["atom_transition_q"], dtype=float
                    ),
                    "grid_normalization_error": abs(
                        float(compact["net_charge_before_neutralization"])
                    ),
                    "compact_path": compact_path,
                }
            )
        frames.append(states)
    return frames


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return abs(float(np.dot(left, right) / denominator)) if denominator else 0.0


def edge_density_overlap(
    site: str,
    edge: int,
    previous_index: int,
    current_index: int,
    frames: list[list[dict]],
) -> float:
    """Return an exact or atomic-proxy adjacent transition-density overlap."""
    if edge not in EXACT_SPATIAL_EDGES[site]:
        return cosine(
            frames[edge][previous_index]["atom_transition_q"],
            frames[edge + 1][current_index]["atom_transition_q"],
        )
    key = (site, edge)
    if key not in _EDGE_MATRIX_CACHE:
        cache_dir = OUTDIR / "adjacent_density_transition_matrices"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"site_{site}_edge_{edge:04d}_{edge + 1:04d}.npy"
        expected_shape = (len(frames[edge]), len(frames[edge + 1]))
        if cache_path.exists():
            matrix = np.load(cache_path)
            if matrix.shape != expected_shape:
                raise RuntimeError(
                    f"Invalid cached edge matrix {cache_path}: "
                    f"{matrix.shape} != {expected_shape}"
                )
        else:
            matrix = np.empty(expected_shape, dtype=float)
            for left_index, left in enumerate(frames[edge]):
                left_density = np.load(left["compact_path"])
                for right_index, right in enumerate(frames[edge + 1]):
                    right_density = np.load(right["compact_path"])
                    matrix[left_index, right_index] = float(
                        full_spatial_overlap(
                            left_density,
                            right_density,
                            sigma=0.12,
                            cutoff=0.35,
                        )["spatial_overlap_absolute"]
                    )
            temporary = cache_path.with_suffix(".npy.tmp")
            with temporary.open("wb") as handle:
                np.save(handle, matrix)
            temporary.replace(cache_path)
            print(
                f"[exact-density-matrix] site={site} edge={edge}->{edge + 1} "
                f"shape={matrix.shape}",
                flush=True,
            )
        _EDGE_MATRIX_CACHE[key] = matrix
    return float(_EDGE_MATRIX_CACHE[key][previous_index, current_index])


def viterbi(
    site: str,
    frames: list[list[dict]],
    *,
    density_weight: float,
    dipole_weight: float,
    energy_weight: float,
    reference_weight: float,
    oscillator_weight: float,
    occupied_weight: float,
    energy_scale_eV: float,
) -> tuple[list[dict], float]:
    def emission(state: dict) -> float:
        # The retained bright anionic chromophore pi-pi* reference is the
        # dominant 93 -> 95 configuration.  Every frame in both sites contains
        # exactly one such candidate.  Treating all occupied-93 excitations as
        # equivalent allowed intensity borrowing from the distinct 93 -> 96
        # state to create unphysical 0.15--0.28 eV one-fs jumps.
        excitation = state["largest_excitation"]
        if not (
            excitation["occupied"] == 93
            and excitation["virtual"] == 95
        ):
            return -1.0e6
        return (
            reference_weight * state["reference_spatial_overlap"]
            + oscillator_weight * np.sqrt(max(state["oscillator_strength"], 0.0))
            + occupied_weight
        )

    scores = [np.asarray([emission(state) for state in frames[0]])]
    parents = []
    for frame in range(1, len(frames)):
        matrix = np.empty((len(frames[frame - 1]), len(frames[frame])))
        for previous_index, previous in enumerate(frames[frame - 1]):
            previous_mu = np.asarray(previous["transition_dipole_au"], dtype=float)
            for current_index, current in enumerate(frames[frame]):
                current_mu = np.asarray(current["transition_dipole_au"], dtype=float)
                matrix[previous_index, current_index] = (
                    scores[-1][previous_index]
                    + density_weight
                    * edge_density_overlap(
                        site,
                        frame - 1,
                        previous_index,
                        current_index,
                        frames,
                    )
                    + dipole_weight * cosine(previous_mu, current_mu)
                    - energy_weight
                    * abs(previous["energy_eV"] - current["energy_eV"])
                    / energy_scale_eV
                    + emission(current)
                )
        parents.append(np.argmax(matrix, axis=0))
        scores.append(np.max(matrix, axis=0))
    index = int(np.argmax(scores[-1]))
    indices = [index]
    for parent in reversed(parents):
        index = int(parent[index])
        indices.append(index)
    indices.reverse()
    return [frames[frame][indices[frame]] for frame in range(len(frames))], float(
        np.max(scores[-1])
    )


def sensitivity(
    site: str, frames: list[list[dict]]
) -> tuple[list[dict], list[dict]]:
    paths = []
    records = []
    for energy_weight in (0.25, 0.5, 1.0, 2.0):
        for density_weight in (1.0, 2.0, 4.0):
            path, score = viterbi(
                site,
                frames,
                density_weight=density_weight,
                dipole_weight=0.3,
                energy_weight=energy_weight,
                reference_weight=0.6,
                oscillator_weight=0.3,
                occupied_weight=0.15,
                energy_scale_eV=0.2,
            )
            paths.append(path)
            records.append(
                {
                    "density_weight": density_weight,
                    "dipole_weight": 0.3,
                    "energy_weight": energy_weight,
                    "reference_weight": 0.6,
                    "oscillator_weight": 0.3,
                    "occupied_93_weight": 0.15,
                    "energy_scale_eV": 0.2,
                    "path_score": score,
                    "roots": [state["root"] for state in path],
                    "maximum_adjacent_energy_jump_eV": max(
                        abs(path[index]["energy_eV"] - path[index - 1]["energy_eV"])
                        for index in range(1, len(path))
                    ),
                }
            )
    canonical, _ = viterbi(
        site,
        frames,
        density_weight=2.0,
        dipole_weight=0.3,
        energy_weight=0.5,
        reference_weight=0.6,
        oscillator_weight=0.3,
        occupied_weight=0.15,
        energy_scale_eV=0.2,
    )
    for frame, state in enumerate(canonical):
        roots = [path[frame]["root"] for path in paths]
        state["sensitivity_consensus_count"] = roots.count(state["root"])
        state["sensitivity_path_count"] = len(paths)
        state["alternative_roots"] = {
            str(root): roots.count(root) for root in sorted(set(roots))
        }
    return canonical, records


def selected_spatial_path_overlap(path: list[dict]) -> list[float]:
    values = [1.0]
    for frame in range(1, len(path)):
        left = np.load(path[frame - 1]["compact_path"])
        right = np.load(path[frame]["compact_path"])
        values.append(
            float(
                full_spatial_overlap(
                    left,
                    right,
                    sigma=0.12,
                    cutoff=0.35,
                )["spatial_overlap_absolute"]
            )
        )
    return values


def main() -> None:
    global WINDOW, MANIFEST, OUTDIR, FRAME_COUNT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--window",
        choices=(
            "short100fs50",
            "short100fs50_fixedmm",
            "short100fs100_wb97xd3",
        ),
        default="short100fs50",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    WINDOW = args.window
    FRAME_COUNT = 100 if WINDOW == "short100fs100_wb97xd3" else 50
    if WINDOW == "short100fs100_wb97xd3":
        # The initial inexpensive proxy pass identified every edge whose
        # selected full-density overlap was below 0.95.  Re-track with complete
        # 7x7 physical-density overlap matrices at those edges and their
        # one-frame shoulders.  This removes reliance on the atomic partition
        # precisely where root identity is least certain.
        EXACT_SPATIAL_EDGES["A"] = {
            1, 2, 3,
            12, 13, 14,
            36, 37, 38,
            41, 42, 43,
            72, 73, 74,
            81, 82, 83,
        }
        EXACT_SPATIAL_EDGES["B"] = {
            13, 14, 15,
            26, 27, 28,
            90, 91, 92,
        }
    MANIFEST = ROOT / "manifests" / f"{WINDOW}_frames.csv"
    OUTDIR = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else ROOT
        / (
            "rapid_100fs50_fixedmm"
            if WINDOW == "short100fs50_fixedmm"
            else "wb97xd3_100fs100"
            if WINDOW == "short100fs100_wb97xd3"
            else "rapid_100fs50"
        )
    )
    calibration = json.loads((OUTDIR / "calibration.json").read_text())
    delta = float(calibration["steom_constant_shift_eV"])
    columns, rows = load_manifest(MANIFEST)
    if any(row["convergence_status"] != "complete" for row in rows):
        raise RuntimeError(
            f"Global state tracking requires {FRAME_COUNT} classified complete frames"
        )
    retained_manifest = OUTDIR / "frame_manifest_static_reference_assignments.csv"
    if not retained_manifest.exists():
        shutil.copy2(MANIFEST, retained_manifest)

    paths = {}
    audit = {
        "status": "running",
        "method": (
            "reference-anchored Viterbi/parallel transport with exact "
            "all-root spatial-density transition matrices at ambiguous edges"
        ),
        "canonical_weights": {
            "adjacent_atomic_transition_density_cosine": 2.0,
            "adjacent_transition_dipole_absolute_cosine": 0.3,
            "adjacent_energy_change_per_0.2_eV": -0.5,
            "static_reference_spatial_overlap": 0.6,
            "sqrt_oscillator_strength": 0.3,
            "required_dominant_configuration": "occupied 93 -> virtual 95",
            "target_configuration_bonus": 0.15,
        },
        "phase_convention": "all compact densities phase-aligned to the retained reference",
        "link_atoms_excluded": [42, 43, 44],
        "exact_spatial_transition_edges": {
            site: sorted(edges) for site, edges in EXACT_SPATIAL_EDGES.items()
        },
        "sensitivity": {},
        "sites": {},
    }
    with retained_manifest.open(newline="") as handle:
        retained_rows = list(csv.DictReader(handle))
    if len(retained_rows) != FRAME_COUNT:
        raise RuntimeError("Retained static-assignment manifest is incomplete")
    for site in ("A", "B"):
        candidates = load_site(site)
        path, variants = sensitivity(site, candidates)
        paths[site] = path
        exact_overlap = selected_spatial_path_overlap(path)
        static_roots = [
            int(row[f"site_{site}_selected_root"]) for row in retained_rows
        ]
        changes = [
            frame
            for frame in range(FRAME_COUNT)
            if static_roots[frame] != path[frame]["root"]
        ]
        audit["sensitivity"][site] = variants
        audit["sites"][site] = {
            "static_roots": static_roots,
            "global_roots": [state["root"] for state in path],
            "changed_frames": changes,
            "selected_adjacent_full_spatial_overlap": exact_overlap,
            "minimum_adjacent_full_spatial_overlap": min(exact_overlap[1:]),
            "median_adjacent_full_spatial_overlap": float(
                np.median(exact_overlap[1:])
            ),
            "maximum_adjacent_energy_jump_eV": max(
                abs(path[index]["energy_eV"] - path[index - 1]["energy_eV"])
                for index in range(1, FRAME_COUNT)
            ),
            "frames": [
                {
                    key: value
                    for key, value in state.items()
                    if key
                    not in {
                        "atom_transition_q",
                        "compact_path",
                    }
                }
                | {"compact_path": str(state["compact_path"].resolve())}
                for state in path
            ],
        }

    for frame in range(FRAME_COUNT):
        frame_dir = ROOT / f"frames/{WINDOW}/frame_{frame:04d}"
        a, b = paths["A"][frame], paths["B"][frame]
        roots_changed = (
            int(retained_rows[frame]["site_A_selected_root"]) != a["root"]
            or int(retained_rows[frame]["site_B_selected_root"]) != b["root"]
        )
        tdc_path = (
            frame_dir
            / f"tdc_global_A{a['root']}_B{b['root']}"
            / "summary.json"
        )
        if roots_changed and not tdc_path.exists():
            run_checked(
                [
                    sys.executable,
                    str(ROOT / "scripts/mutual_tdc.py"),
                    "--site-a",
                    str(a["compact_path"]),
                    "--site-b",
                    str(b["compact_path"]),
                    "--out",
                    str(tdc_path),
                    "--epsilon",
                    "1.77",
                    "--box-lengths",
                    *[str(value) for value in BOX_A],
                ]
            )
        tdc_source = tdc_path if roots_changed else frame_dir / "tdc/summary.json"
        tdc = json.loads(tdc_source.read_text())
        if tdc["status"] != "complete":
            raise RuntimeError(f"Global TDC validation failed at frame {frame}")
        global_assignment = {
            "frame_index": frame,
            "selection_method": "global reference-anchored Viterbi/parallel transport",
            "site_A": {
                key: value
                for key, value in a.items()
                if key not in {"atom_transition_q", "compact_path"}
            }
            | {"compact_path": str(a["compact_path"].resolve())},
            "site_B": {
                key: value
                for key, value in b.items()
                if key not in {"atom_transition_q", "compact_path"}
            }
            | {"compact_path": str(b["compact_path"].resolve())},
            "roots_changed_from_static_reference_assignment": roots_changed,
            "tdc_source": str(tdc_source.resolve()),
        }
        assignment_path = frame_dir / "state_assignment_global.json"
        assignment_path.write_text(json.dumps(global_assignment, indent=2) + "\n")
        summary_path = frame_dir / "frame_summary_global.json"
        summary_path.write_text(
            json.dumps(
                {
                    "window": WINDOW,
                    "frame_index": frame,
                    "global_assignment": global_assignment,
                    "tdc": tdc,
                    "steom_constant_shift_eV": delta,
                    "original_frame_summary": rows[frame]["compact_output_path"],
                },
                indent=2,
            )
            + "\n"
        )
        rows[frame].update(
            {
                "site_A_energy_eV": f"{a['energy_eV']:.10f}",
                "site_A_energy_cm-1": f"{a['energy_cm-1']:.6f}",
                "site_A_wavelength_nm": f"{a['wavelength_nm']:.6f}",
                "site_B_energy_eV": f"{b['energy_eV']:.10f}",
                "site_B_energy_cm-1": f"{b['energy_cm-1']:.6f}",
                "site_B_wavelength_nm": f"{b['wavelength_nm']:.6f}",
                "site_A_corrected_eV": f"{a['energy_eV'] + delta:.10f}",
                "site_B_corrected_eV": f"{b['energy_eV'] + delta:.10f}",
                "site_A_oscillator_strength": f"{a['oscillator_strength']:.8f}",
                "site_B_oscillator_strength": f"{b['oscillator_strength']:.8f}",
                "site_A_selected_root": str(a["root"]),
                "site_B_selected_root": str(b["root"]),
                "site_A_assignment_score": f"{a['reference_spatial_overlap']:.8f}",
                "site_B_assignment_score": f"{b['reference_spatial_overlap']:.8f}",
                "J_vacuum_cm-1": f"{tdc['J_vacuum_cm-1']:.8f}",
                "J_screened_cm-1": f"{tdc['J_screened_cm-1']:.8f}",
                "J_pda_vacuum_cm-1": f"{tdc['J_pda_vacuum_cm-1']:.8f}",
                "J_pda_screened_cm-1": f"{tdc['J_pda_screened_cm-1']:.8f}",
                "separation_A": f"{tdc['separation_A']:.8f}",
                "dipole_angle_deg": f"{tdc['dipole_angle_deg']:.8f}",
                "grid_normalization_error_A": f"{a['grid_normalization_error']:.8e}",
                "grid_normalization_error_B": f"{b['grid_normalization_error']:.8e}",
                "compact_output_path": str(summary_path.resolve()),
                "notes": (
                    "global density/configuration/brightness/dipole/energy state "
                    "path; signed TDC; original static assignment retained; no interpolation"
                ),
            }
        )
        save_manifest(MANIFEST, columns, rows)
        print(
            f"[global] frame={frame} roots={a['root']}/{b['root']} "
            f"changed={roots_changed} J={tdc['J_screened_cm-1']:.3f}",
            flush=True,
        )
    audit["status"] = "complete"
    audit["manifest"] = str(MANIFEST.resolve())
    audit["retained_static_manifest"] = str(retained_manifest.resolve())
    (OUTDIR / "global_state_tracking.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print("[global-state-tracking-complete]", flush=True)


if __name__ == "__main__":
    main()

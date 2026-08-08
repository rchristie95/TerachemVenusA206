#!/usr/bin/env python3
"""Validate compact published data without ORCA, OpenMM, PyMOL, or a GPU."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def file_sha256(path: Path) -> str:
    """Return the byte-exact SHA-256 used for generated ORCA outputs."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_last_electric_spectrum(path: Path) -> list[dict[str, object]]:
    """Parse the final ORCA electric-dipole absorption table.

    A STEOM output can contain more than one spectrum.  The published state is
    taken from the final electric-dipole table, so this parser deliberately
    returns the last complete non-empty block rather than the first match.
    """
    heading = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    end_heading = "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
    blocks: list[list[dict[str, object]]] = []
    current: list[dict[str, object]] | None = None
    row_pattern = re.compile(
        rf"^\s*(\S+)\s+->\s+(\S+)\s+({FLOAT})\s+({FLOAT})\s+"
        rf"({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+"
        rf"({FLOAT})\s*$"
    )
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if heading in line:
            if current:
                blocks.append(current)
            current = []
            continue
        if current is not None and end_heading in line:
            if current:
                blocks.append(current)
            current = None
            continue
        if current is None:
            continue
        match = row_pattern.match(line)
        if not match:
            continue
        target = match.group(2)
        root_match = re.match(r"(\d+)-", target)
        current.append(
            {
                "transition": f"{match.group(1)} -> {target}",
                "root": int(root_match.group(1)) if root_match else target,
                "energy_eV": float(match.group(3)),
                "wavenumber_cm-1": float(match.group(4)),
                "wavelength_nm": float(match.group(5)),
                "oscillator_strength": float(match.group(6)),
                "D2_au2": float(match.group(7)),
                "transition_dipole_au": [
                    float(match.group(8)),
                    float(match.group(9)),
                    float(match.group(10)),
                ],
            }
        )
    if current:
        blocks.append(current)
    if not blocks:
        raise AssertionError(f"No electric-dipole absorption spectrum in {path}")
    return blocks[-1]


def parse_final_scf_energy(path: Path) -> float:
    pattern = re.compile(rf"FINAL SINGLE POINT ENERGY\s+({FLOAT})")
    matches = pattern.findall(path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        raise AssertionError(f"No final single-point energy in {path}")
    return float(matches[-1])


def parse_final_cc_total_energy(path: Path) -> float:
    pattern = re.compile(rf"^E\(TOT\)\s+\.{{3}}\s+({FLOAT})\s*$", re.MULTILINE)
    matches = pattern.findall(path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        raise AssertionError(f"No final coupled-cluster total energy in {path}")
    return float(matches[-1])


def parse_t1_diagnostic(path: Path) -> float:
    pattern = re.compile(rf"^T1 diagnostic\s+\.{{3}}\s+({FLOAT})\s*$", re.MULTILINE)
    matches = pattern.findall(path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        raise AssertionError(f"No T1 diagnostic in {path}")
    return float(matches[-1])


def validate_checksums() -> None:
    for raw in (ROOT / "reference/SHA256SUMS").read_text().splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        expected, rel = raw.split(maxsplit=1)
        path = ROOT / rel
        if not path.exists():
            print(f"SKIP external/generated: {rel}")
            continue
        data = path.read_bytes()
        # Git may materialise text files as CRLF on Windows. SHA256SUMS records
        # the repository-canonical LF bytes so validation is cross-platform;
        # binary assets must remain byte-exact (normalising arbitrary binary
        # bytes can silently change a hash whenever b"\r\n" occurs in them).
        if path.suffix.lower() in {
            ".csv",
            ".dat",
            ".in",
            ".inp",
            ".inpcrd",
            ".json",
            ".pc",
            ".pdb",
            ".prmtop",
            ".txt",
            ".xyz",
        }:
            data = data.replace(b"\r\n", b"\n")
        actual = hashlib.sha256(data).hexdigest()
        assert actual == expected, (
            f"checksum mismatch: {rel}; expected {expected}, obtained {actual}"
        )


def validate_tandem_statistics() -> None:
    ref = json.loads((ROOT / "reference/orca_validation.json").read_text())
    production = ROOT / "coupling_nvt_production_cr2_1000_20260721/coupling_samples.csv"
    with production.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = [float(row["J_cm"]) for row in rows]
    pda_values = [float(row["J_pda_cm"]) for row in rows]
    separations = [float(row["separation_A"]) for row in rows]
    angles = [float(row["angle_deg"]) for row in rows]
    expected = ref["tandem_ensemble"]
    assert len(values) == expected["frames"] == 1000
    assert abs(statistics.fmean(values) - expected["mean_J_cm-1"]) < 1e-12
    # coupling_ensemble.py reports the sample standard deviation (ddof=1).
    assert abs(statistics.stdev(values) - expected["std_J_cm-1"]) < 1e-12
    assert min(values) == expected["min_J_cm-1"]
    assert max(values) == expected["max_J_cm-1"]
    assert abs(statistics.fmean(pda_values) - expected["mean_PDA_cm-1"]) < 1e-12
    assert abs(statistics.stdev(pda_values) - expected["std_PDA_cm-1"]) < 1e-12
    assert abs(
        statistics.fmean(values) / statistics.fmean(pda_values)
        - expected["mean_TDC_over_mean_PDA"]
    ) < 1e-12
    assert abs(2.0 * statistics.fmean(values) - expected["mean_splitting_cm-1"]) < 1e-12
    assert abs(2.0 * statistics.stdev(values) - expected["std_splitting_cm-1"]) < 1e-12
    assert abs(statistics.fmean(separations) - expected["mean_separation_A"]) < 1e-12
    assert abs(statistics.stdev(separations) - expected["std_separation_A"]) < 1e-12
    assert abs(statistics.fmean(angles) - expected["mean_dipole_angle_deg"]) < 1e-12
    assert abs(statistics.stdev(angles) - expected["std_dipole_angle_deg"]) < 1e-12
    assert [int(row["frame"]) for row in rows] == list(range(1000))


def validate_qm_inputs() -> None:
    geom = (ROOT / "neo_model/orca_steom/geom_cthrp.xyz").read_text().splitlines()
    field = (ROOT / "neo_model/orca_steom/field.pc").read_text().splitlines()
    assert int(geom[0]) == 44
    assert int(field[0]) == 2350
    assert len(geom[2:]) == 44
    assert len(field[1:]) == 2350


def validate_orca_outputs_if_present() -> None:
    """Cross-check locally retained ORCA outputs against compact references.

    The large licensed-program outputs are intentionally not released, so a
    clean clone skips these checks.  On the production workstation, however,
    their byte hashes and parsed bright-state rows must match the reference
    record exactly enough to catch a stale or partially written rerun.
    """
    ref = json.loads((ROOT / "reference/orca_validation.json").read_text())
    controls: list[dict[str, object]] = []
    embedded = ref.get("orca_tddft_bright_state", {})
    if isinstance(embedded, dict) and embedded.get("source_output"):
        controls.append(embedded)
    gas = ref.get("gas_phase_44_atom_controls", {})
    if isinstance(gas, dict):
        for name in ("tddft_bright_state", "steom_bright_state"):
            record = gas.get(name, {})
            if isinstance(record, dict) and record.get("source_output"):
                controls.append(record)

    tolerances = {
        "energy_eV": 5.1e-7,
        "wavenumber_cm-1": 5.1e-2,
        "wavelength_nm": 5.1e-2,
        "oscillator_strength": 5.1e-10,
    }
    for expected in controls:
        path = ROOT / str(expected["source_output"])
        if not path.exists():
            print(f"SKIP external/generated ORCA output: {path.relative_to(ROOT)}")
            continue
        assert file_sha256(path) == expected["output_sha256"], path
        spectrum = parse_last_electric_spectrum(path)
        bright = max(spectrum, key=lambda row: float(row["oscillator_strength"]))
        assert bright["root"] == expected["root"], (path, bright)
        for key, tolerance in tolerances.items():
            assert abs(float(bright[key]) - float(expected[key])) <= tolerance, (
                path,
                key,
                bright[key],
                expected[key],
            )
        if "transition_dipole_au" in expected:
            observed_vector = np.asarray(bright["transition_dipole_au"], dtype=float)
            expected_vector = np.asarray(expected["transition_dipole_au"], dtype=float)
            assert np.allclose(observed_vector, expected_vector, rtol=0.0, atol=5.1e-6)
            assert abs(
                np.linalg.norm(observed_vector)
                - float(expected["transition_dipole_norm_au"])
            ) < 1.0e-12
        if "final_scf_energy_Eh" in expected:
            assert abs(
                parse_final_scf_energy(path) - float(expected["final_scf_energy_Eh"])
            ) < 5.0e-12
        if "final_cc_total_energy_Eh" in expected:
            assert abs(
                parse_final_cc_total_energy(path)
                - float(expected["final_cc_total_energy_Eh"])
            ) < 5.0e-12
        if "t1_diagnostic" in expected:
            assert abs(
                parse_t1_diagnostic(path) - float(expected["t1_diagnostic"])
            ) < 5.0e-12
        if "terminated_normally" in expected:
            normal = "ORCA TERMINATED NORMALLY" in path.read_text(
                encoding="utf-8", errors="replace"
            )
            assert normal is bool(expected["terminated_normally"]), path


def validate_independent_controls_if_present() -> None:
    """Check the compact ORCA and independent CPU/GPU audit records."""
    ref = json.loads((ROOT / "reference/orca_validation.json").read_text())
    expected = ref.get("independent_crystal_geometry_control")
    if not isinstance(expected, dict):
        return

    orca_path = ROOT / "tmp/validation/cross_code_control/orca_tddft_control.json"
    if orca_path.exists():
        observed = json.loads(orca_path.read_text())
        assert observed["method"] == expected["orca"]["method"]
        assert expected["geometry"] == (
            "two anionic 54-atom fragments at a crystal-derived geometry"
        )
        assert observed["centre_definition"] == expected["centre_definition"]
        for key in (
            "threshold_e_per_bohr3",
            "centre_separation_angstrom",
            "tdc_hartree",
            "tdc_meV",
            "tdc_cm-1",
            "pda_hartree",
            "pda_meV",
            "pda_cm-1",
            "abs_tdc_over_pda",
        ):
            assert float(observed[key]) == float(expected["orca"][key]), key
        assert [fragment["natoms"] for fragment in observed["fragments"]] == [54, 54]
        assert expected["fragments"] == 2
        assert expected["atoms_per_fragment"] == 54
        assert expected["fragment_charge"] == -1
        assert expected["multiplicity"] == 1

    crosscheck_path = ROOT / "tmp/validation/cross_code_control/cpu_gpu_crosscheck.json"
    if crosscheck_path.exists():
        observed = json.loads(crosscheck_path.read_text())
        cpu_gpu = expected["independent_cpu_gpu_crosscheck"]
        for key in (
            "threshold_e_per_bohr3",
            "coarse_bin_width_angstrom",
            "raw_sum_inverse_angstrom",
            "numpy_hartree",
            "opencl_hartree",
            "relative_difference",
            "numpy_meV",
        ):
            assert float(observed[key]) == float(cpu_gpu[key]), key
        assert observed["source_points"] == cpu_gpu["source_points"]
        assert observed["coarse_points"] == cpu_gpu["coarse_points"]


def validate_static_and_multipole() -> None:
    """Keep the static control, multipole table, and density provenance locked."""
    ref = json.loads((ROOT / "reference/orca_validation.json").read_text())

    static_csv = ROOT / "coupling_paper_steom_static/coupling_samples.csv"
    with static_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    row = rows[0]
    expected = ref["static_geometry"]
    field_map = {
        "J_cm": "J_TDC_cm-1",
        "J_pda_cm": "J_PDA_cm-1",
        "separation_A": "separation_A",
        "angle_deg": "dipole_angle_deg",
        "aln_A_rms": "alignment_A_rms_A",
        "aln_B_rms": "alignment_B_rms_A",
    }
    for csv_key, ref_key in field_map.items():
        assert float(row[csv_key]) == expected[ref_key], (csv_key, ref_key)
    assert abs(
        abs(float(row["J_cm"]) / float(row["J_pda_cm"]))
        - expected["TDC_over_PDA"]
    ) < 1.0e-12

    distribution = json.loads(
        (ROOT / "coupling_paper_steom_static/coupling_distribution.json").read_text()
    )
    assert distribution["mean"] == expected["J_TDC_cm-1"]
    assert distribution["samples"] == [expected["J_TDC_cm-1"]]
    assert distribution["tdc_units"] == {
        "status": "corrected",
        "pair_distance_unit": "angstrom",
        "reciprocal_distance_to_atomic_units": 0.529177210903,
    }
    provenance = distribution["density_provenance"]
    density_ref = ref["steom_density"]
    assert provenance["retained_points"] == density_ref["retained_points"] == 259277
    assert provenance["included_nto_pairs"] == len(density_ref["nto_pairs"]) == 7
    assert abs(
        provenance["represented_nto_occupation"]
        - density_ref["represented_nto_occupation"]
    ) < 1.0e-12
    assert provenance["source"].endswith("steom_transdens_capmasked_oldframe.npz")

    dipole_geometry = json.loads(
        (ROOT / "coupling_paper_steom_thermal/dipole_geometry.json").read_text()
    )
    assert dipole_geometry["separation_A"] == expected["separation_A"]
    assert dipole_geometry["angle_deg"] == expected["dipole_angle_deg"]
    assert dipole_geometry["aln_A_rms"] == expected["alignment_A_rms_A"]
    assert dipole_geometry["aln_B_rms"] == expected["alignment_B_rms_A"]
    assert dipole_geometry["source_density"].endswith(
        "steom_transdens_capmasked_oldframe.npz"
    )

    multipole_rows: dict[str, dict[str, str]] = {}
    with (ROOT / "multipole_out_correct/multipole_analysis.csv").open(newline="") as handle:
        for multipole_row in csv.DictReader(handle):
            multipole_rows[multipole_row["term"]] = multipole_row
    multipole_ref = ref["multipole"]
    checks = {
        "dip-dip (PDA)": ("contribution_cm", "dipole_dipole_cm-1"),
        "dip-quad": ("contribution_cm", "dipole_quadrupole_cm-1"),
        "quad-quad": ("contribution_cm", "quadrupole_quadrupole_cm-1"),
        "dip-oct": ("contribution_cm", "dipole_octupole_cm-1"),
        "multipole_sum": ("cumulative_cm", "multipole_sum_cm-1"),
        "full_TDC": ("cumulative_cm", "full_TDC_cm-1"),
        "TDC_over_PDA": ("cumulative_cm", "TDC_over_PDA"),
    }
    for term, (csv_key, ref_key) in checks.items():
        assert abs(float(multipole_rows[term][csv_key]) - multipole_ref[ref_key]) < 5.0e-5

    # These deterministic products are intentionally optional in a clean clone,
    # but when present they must match the released production-density shape.
    for relative in (
        "neo_model/orca_steom/steom_transdens_capmasked.npz",
        "neo_model/orca_steom/steom_transdens_capmasked_oldframe.npz",
    ):
        path = ROOT / relative
        if path.exists():
            with np.load(path, allow_pickle=False) as archive:
                assert archive["pts_ang"].shape == (259277, 3)
                assert archive["q"].shape == (259277,)
                assert archive["mu_au"].shape == (3,)
                points = np.asarray(archive["pts_ang"], dtype=float)
                charges = np.asarray(archive["q"], dtype=float)
                assert abs(charges.sum() - density_ref["net_transition_charge_au"]) < 1.0e-15
                assert abs(
                    np.abs(charges).sum() - density_ref["total_absolute_grid_weight_e"]
                ) < 1.0e-15
                assert abs(
                    abs(charges.sum()) / np.abs(charges).sum()
                    - density_ref["residual_transition_charge_fraction"]
                ) < 1.0e-15
                centered_mu = np.sum(
                    charges[:, None] * (points - points.mean(axis=0)), axis=0
                ) / 0.529177210903
                assert abs(
                    np.linalg.norm(centered_mu)
                    - density_ref["production_centered_dipole_norm_au"]
                ) < 1.0e-12

                # Rigid alignment may rotate the archived metadata vector but
                # must preserve its norm exactly.
                assert abs(
                    np.linalg.norm(archive["mu_au"])
                    - density_ref["archive_metadata_pre_cutoff_dipole_norm_au"]
                ) < 1.0e-12

                if path.name == "steom_transdens_capmasked.npz":
                    retained_source_mu = -np.sum(charges[:, None] * points, axis=0) / 0.529177210903
                    expected_retained_source_mu = np.asarray(
                        density_ref["retained_source_frame_dipole_au"], dtype=float
                    )
                    assert np.allclose(
                        retained_source_mu,
                        expected_retained_source_mu,
                        rtol=0.0,
                        atol=1.0e-12,
                    )
                    assert abs(
                        np.linalg.norm(retained_source_mu)
                        - density_ref["retained_source_frame_dipole_norm_au"]
                    ) < 1.0e-12
                    target_mu = np.asarray(archive["mu_target_au"], dtype=float)
                    norm_error = abs(
                        np.linalg.norm(retained_source_mu) - np.linalg.norm(target_mu)
                    ) / np.linalg.norm(target_mu)
                    assert abs(
                        norm_error
                        - density_ref["retained_source_frame_relative_norm_error"]
                    ) < 1.0e-15
                    cosine = float(
                        np.dot(retained_source_mu, target_mu)
                        / (np.linalg.norm(retained_source_mu) * np.linalg.norm(target_mu))
                    )
                    assert abs(
                        cosine
                        - density_ref["retained_source_frame_dipole_cosine_to_target"]
                    ) < 1.0e-15
                    expected_pairs = np.asarray(density_ref["nto_pairs"], dtype=float)
                    assert np.array_equal(
                        np.asarray(archive["nto_pairs"], dtype=float)[:, :3],
                        expected_pairs,
                    )
                    assert tuple(archive["source_grid_shape"].tolist()) == tuple(
                        density_ref["source_grid_shape"]
                    )
                    spacing = np.linalg.norm(
                        np.asarray(archive["source_grid_axes_bohr"], dtype=float), axis=1
                    ) * 0.529177210903
                    assert np.allclose(
                        spacing,
                        np.asarray(density_ref["source_grid_spacing_A"], dtype=float),
                        rtol=0.0,
                        atol=1.0e-12,
                    )
                    assert abs(
                        np.linalg.norm(archive["mu_unmasked_au"])
                        - density_ref["unmasked_dipole_norm_au"]
                    ) < 1.0e-12
                    assert abs(
                        np.linalg.norm(archive["mu_masked_unscaled_au"])
                        - density_ref["cap_masked_unscaled_dipole_norm_au"]
                    ) < 1.0e-12
                    assert abs(
                        np.linalg.norm(archive["mu_target_au"])
                        - density_ref["target_right_transition_dipole_norm_au"]
                    ) < 1.0e-12
                    assert str(archive["mask_method"]) in density_ref["link_cap_mask"]
                    assert abs(
                        float(archive["retained_grid_fraction"])
                        - density_ref["retained_grid_fraction_before_threshold"]
                    ) < 1.0e-12


if __name__ == "__main__":
    validate_checksums()
    validate_tandem_statistics()
    validate_qm_inputs()
    validate_orca_outputs_if_present()
    validate_independent_controls_if_present()
    validate_static_and_multipole()
    print("reference validation: PASS")

#!/usr/bin/env python3
"""Guard the frozen-geometry gas/QM-MM comparison against method drift."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def substantive_lines(relative: str, *, omit_point_charges: bool = False) -> list[str]:
    lines = []
    for raw in (ROOT / relative).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if omit_point_charges and line.lower().startswith("%pointcharges"):
            continue
        lines.append(line)
    return lines


def test_tddft_gas_control_changes_only_embedding() -> None:
    embedded = "neo_model/orca_dft/tddft_wb97xd3.inp"
    gas = "neo_model/orca_dft/tddft_wb97xd3_gas.inp"
    assert substantive_lines(embedded, omit_point_charges=True) == substantive_lines(gas)
    assert any(
        line.lower().startswith("%pointcharges")
        for line in substantive_lines(embedded)
    )
    assert not any(
        line.lower().startswith("%pointcharges") for line in substantive_lines(gas)
    )


def test_steom_gas_control_changes_only_embedding() -> None:
    embedded = "neo_model/orca_steom/steom_phenol_svpd_robust2.inp"
    gas = "neo_model/orca_steom/steom_phenol_svpd_gas.inp"
    assert substantive_lines(embedded, omit_point_charges=True) == substantive_lines(gas)
    assert any(
        line.lower().startswith("%pointcharges")
        for line in substantive_lines(embedded)
    )
    assert not any(
        line.lower().startswith("%pointcharges") for line in substantive_lines(gas)
    )


def test_gas_controls_use_the_exact_frozen_44_atom_anion() -> None:
    geometry = (ROOT / "neo_model/orca_steom/geom_cthrp.xyz").read_text().splitlines()
    assert int(geometry[0]) == 44
    assert len(geometry[2:]) == 44
    for relative in (
        "neo_model/orca_dft/tddft_wb97xd3_gas.inp",
        "neo_model/orca_steom/steom_phenol_svpd_gas.inp",
    ):
        lines = substantive_lines(relative)
        assert lines[-1].lower() == "* xyzfile -1 1 geom_cthrp.xyz"


def test_qm_region_figure_uses_the_calculation_geometry() -> None:
    def xyz(relative: str) -> tuple[list[str], np.ndarray]:
        lines = (ROOT / relative).read_text().splitlines()
        count = int(lines[0])
        records = [line.split() for line in lines[2:]]
        assert len(records) == count == 44
        return [record[0] for record in records], np.asarray(
            [[float(value) for value in record[1:4]] for record in records]
        )

    calculation_elements, calculation_coordinates = xyz(
        "neo_model/orca_steom/geom_cthrp.xyz"
    )
    figure_elements, figure_coordinates = xyz("neo_model/orca_steom/steom_qm.xyz")
    assert figure_elements == calculation_elements
    assert float(np.max(np.abs(figure_coordinates - calculation_coordinates))) < 1.0e-5
    assert (ROOT / "manuscript/Fig_QM_Region_44.png").is_file()


def test_frozen_geometry_records_its_composite_qmmm_provenance() -> None:
    """The physical CR2 coordinates must be the retained constrained-QM/MM result."""
    opt_lines = (ROOT / "tc_qmmm_opt_constrained/opt.in").read_text().splitlines()
    opt_directives = {line.strip().lower() for line in opt_lines}
    assert {
        "run minimize",
        "basis 6-31g*",
        "method wb97xd3",
        "pointcharges mm_charges.dat",
        "min_coordinates cartesian",
    } <= opt_directives
    constrained = {
        int(line.split()[1])
        for line in opt_lines
        if line.strip().lower().startswith("atom ")
    }
    assert set(range(1, 275)) - constrained == set(range(50, 79))
    freeze_record = {
        int(value)
        for value in (
            ROOT / "tc_qmmm_opt_constrained/cr2only_freeze.txt"
        ).read_text().strip().split(",")
    }
    assert freeze_record == constrained

    qm_lines = (ROOT / "tc_qmmm_opt_constrained/qm_opt.xyz").read_text().splitlines()
    frozen_lines = (ROOT / "neo_model/orca_steom/geom_cthrp.xyz").read_text().splitlines()
    qm_records = [line.split() for line in qm_lines[2:]]
    frozen_records = [line.split() for line in frozen_lines[2:]]

    # opt.in freezes one-based atoms 1--49 and 79--274, leaving the 29
    # physical CR2 atoms 50--78.  The composite builder places those atoms
    # first and appends the three rebuilt link hydrogens as atoms 42--44.
    relaxed_cr2 = qm_records[49:78]
    frozen_cr2 = frozen_records[:29]
    assert [row[0] for row in relaxed_cr2] == [row[0] for row in frozen_cr2]
    assert np.allclose(
        np.asarray([[float(value) for value in row[1:4]] for row in relaxed_cr2]),
        np.asarray([[float(value) for value in row[1:4]] for row in frozen_cr2]),
        rtol=0.0,
        atol=1.0e-8,
    )
    assert [row[0] for row in frozen_records[41:44]] == ["H", "H", "H"]

    validation = json.loads((ROOT / "reference/orca_validation.json").read_text())
    model = validation["qm_model"]
    assert model["complete_44_atom_model_reoptimized"] is False
    assert "constrained TeraChem QM/MM" in model["geometry_provenance"]

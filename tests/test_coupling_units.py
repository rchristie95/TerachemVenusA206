#!/usr/bin/env python3
"""Regression tests for Coulomb inverse-distance unit conversion."""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import coupling_core as cc


def test_inverse_angstrom_conversion() -> None:
    # One unit-charge pair separated by exactly 1 Angstrom has energy
    # 1/(1 Angstrom in bohr) = 0.529177... Hartree.
    converted = cc.inverse_angstrom_to_hartree(1.0)
    assert math.isclose(converted, cc.BOHR_TO_ANGSTROM, rel_tol=0.0, abs_tol=1.0e-15)
    assert math.isclose(converted, 0.529177210903, rel_tol=0.0, abs_tol=1.0e-12)


def test_coordinate_and_reciprocal_factors_are_inverses() -> None:
    assert math.isclose(
        cc.ANGSTROM_TO_BOHR * cc.BOHR_TO_ANGSTROM, 1.0, rel_tol=0.0, abs_tol=1.0e-15
    )
    # The old bug used ANGSTROM_TO_BOHR for a reciprocal distance and inflated
    # every TDC value by this exact factor.
    inflation = cc.ANGSTROM_TO_BOHR / cc.BOHR_TO_ANGSTROM
    assert math.isclose(inflation, 3.571064826093, rel_tol=0.0, abs_tol=1.0e-9)


def test_neutral_dipoles_recover_pda_at_long_range() -> None:
    # Two identical charge-neutral dipoles perpendicular to their separation.
    # At long range, the exact four-charge Coulomb sum must approach the PDA.
    bond_a = 0.02
    separation_a = 100.0
    charges = np.asarray([-1.0, 1.0])
    points_a = np.asarray([[0.0, -bond_a / 2, 0.0], [0.0, bond_a / 2, 0.0]])
    points_b = points_a + np.asarray([separation_a, 0.0, 0.0])
    distance_a = np.linalg.norm(points_a[:, None, :] - points_b[None, :, :], axis=2)
    exact = cc.inverse_angstrom_to_hartree(
        float(np.sum(charges[:, None] * charges[None, :] / distance_a))
    )
    mu = bond_a * cc.ANGSTROM_TO_BOHR
    separation_bohr = separation_a * cc.ANGSTROM_TO_BOHR
    pda = mu**2 / separation_bohr**3
    assert math.isclose(exact, pda, rel_tol=1.0e-7, abs_tol=1.0e-18), (exact, pda)


def test_opencl_one_angstrom_pair_if_available() -> None:
    if not cc._is_opencl_ready():
        return
    value = cc.calculate_coupling_opencl(
        np.asarray([[0.0, 0.0, 0.0]]),
        np.asarray([1.0]),
        np.asarray([[1.0, 0.0, 0.0]]),
        np.asarray([1.0]),
    )
    assert math.isclose(value, cc.BOHR_TO_ANGSTROM, rel_tol=1.0e-12, abs_tol=1.0e-12)


if __name__ == "__main__":
    test_inverse_angstrom_conversion()
    test_coordinate_and_reciprocal_factors_are_inverses()
    test_neutral_dipoles_recover_pda_at_long_range()
    test_opencl_one_angstrom_pair_if_available()
    print("ALL PASSED (4 tests)")

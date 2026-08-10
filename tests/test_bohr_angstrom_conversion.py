"""Guard the Angstrom->bohr conversion that produced a 3.5711x coupling error.

A Coulomb sum ``sum(q_i q_j / r_ij)`` evaluated with r in Angstrom is converted
to Hartree by multiplying by BOHR_TO_ANGSTROM (0.529177), because

    r_bohr = r_angstrom / BOHR_TO_ANGSTROM   =>   1/r_bohr = BOHR_TO_ANGSTROM / r_angstrom

Multiplying instead by ANGSTROM_TO_BOHR (1.8897) - the correct factor for
converting a COORDINATE, not a reciprocal distance - inflates every coupling by
1 / BOHR_TO_ANGSTROM**2 = 3.5711. That error shipped once in this repository
(corrected in commit e1c42f8) and independently in the decoherence repository,
so it is guarded numerically and at source level here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import coupling_core as cc

ROOT = Path(__file__).resolve().parents[1]
INFLATION = 1.0 / cc.BOHR_TO_ANGSTROM**2


def test_constants_are_reciprocal():
    assert cc.BOHR_TO_ANGSTROM == pytest.approx(0.529177210903, rel=0, abs=1e-12)
    assert cc.ANGSTROM_TO_BOHR == pytest.approx(1.0 / cc.BOHR_TO_ANGSTROM, rel=1e-15)
    assert INFLATION == pytest.approx(3.5711, abs=1e-4)


def test_inverse_angstrom_helper_matches_closed_form():
    """Two unit charges 10 A apart: E = 1/(10/0.529177) Hartree."""
    coulomb_sum_per_angstrom = 1.0 * 1.0 / 10.0
    hartree = cc.inverse_angstrom_to_hartree(coulomb_sum_per_angstrom)
    exact = 1.0 / (10.0 / cc.BOHR_TO_ANGSTROM)
    assert hartree == pytest.approx(exact, rel=1e-12)
    # independent route: 14.399645 eV.A / 10 A, expressed in cm^-1
    assert hartree * cc.HARTREE_TO_CM == pytest.approx(14.399645 / 10.0 * 8065.54, rel=1e-4)


def test_using_the_coordinate_factor_would_inflate_by_the_known_amount():
    coulomb_sum_per_angstrom = 1.0 * 1.0 / 10.0
    correct = cc.inverse_angstrom_to_hartree(coulomb_sum_per_angstrom)
    wrong = coulomb_sum_per_angstrom * cc.ANGSTROM_TO_BOHR
    assert wrong / correct == pytest.approx(INFLATION, rel=1e-12)


def test_coupling_kernels_route_through_the_named_helper():
    """The GPU/OpenCL kernels must not open-code the conversion."""
    source = (ROOT / "coupling_core.py").read_text()
    assert "return j_sum * ANGSTROM_TO_BOHR" not in source
    assert source.count("inverse_angstrom_to_hartree(j_sum)") >= 2


def _python_sources():
    skip = {".git", "__pycache__", "orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"}
    for path in ROOT.rglob("*.py"):
        if any(part in skip for part in path.parts):
            continue
        yield path


def test_no_source_multiplies_a_hartree_conversion_by_the_coordinate_factor():
    """ANGSTROM_TO_BOHR and a Hartree conversion on one line is the bug signature.

    Converting a reciprocal distance needs BOHR_TO_ANGSTROM; a line doing both
    ANGSTROM_TO_BOHR and HARTREE_TO_CM is almost certainly the 3.5711x error.
    """
    # Only MULTIPLICATION by the coordinate factor alongside a Hartree
    # conversion is the bug. Dividing by ANGSTROM_TO_BOHR is equivalent to
    # multiplying by BOHR_TO_ANGSTROM and is correct (analyze_solvation_
    # decoherence.py does exactly that).
    pattern = re.compile(
        r"(\*\s*ANGSTROM_TO_BOHR[^/\n]*HARTREE)|(HARTREE[^/\n]*\*\s*ANGSTROM_TO_BOHR)"
    )
    offenders = []
    for path in _python_sources():
        if path.resolve() == Path(__file__).resolve():
            continue
        for number, line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith(("#", "import ", "from ")) or stripped.endswith(","):
                continue
            if pattern.search(line) and "inverse_angstrom_to_hartree" not in line:
                offenders.append(f"{path.relative_to(ROOT)}:{number}: {stripped}")
    assert not offenders, "reciprocal-distance conversion bug signature:\n" + "\n".join(offenders)

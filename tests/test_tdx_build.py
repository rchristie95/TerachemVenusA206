"""Guard the dVenus-TDX construct: one deletion in the sequence, one pair in the structure.

TDX is TD with one barrel's chromophore extinguished.  Addgene #250939 does it
by DELETING the chromophore tyrosine, not by mutating it -- verified against the
deposit on 2026-09-02 and pinned here, because the difference is not cosmetic:
the chain is 510 aa rather than 511, the dark barrel's triad is Gly-Gly rather
than Gly-Xaa-Gly, and every downstream residue number shifts by one.

Three things can silently go wrong and none would be visible in a finished
trajectory:

  * the sequence drifting from the deposit, which would make a TD-vs-TDX
    difference uninterpretable;
  * the dead pair being grafted from a model of the wrong construct -- a TD
    model, or the superseded 511-aa Y->C reading -- which would quietly put a
    residue back where the deposit has none;
  * the graft disturbing the inter-barrel register, which is the one thing the
    crystal structure is there to supply and AlphaFold is not.

The structural tests drive the AlphaFold graft path against a synthetic model
built from the v4 structure itself, so the correct answer is known exactly: the
grafted pair must come back where it started.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "tandem_tdx_build" / "build_tdx_structure.py"
GENERATOR = ROOT / "alphafold_tandem" / "make_tdx_sequence.py"
PARENT = ROOT / "alphafold_tandem" / "pdVenus-TD_Addgene250938_511aa.fasta"
DEPOSIT = ROOT / "alphafold_tandem" / "pdVenus-TDX_Addgene250939_deposit_510aa.fasta"
V4 = ROOT / "tandem_v4_correct_seq" / "tandem_dimer_v4_af3tails.pdb"


def _module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build = _module(BUILDER, "build_tdx_structure")


def _fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


# --------------------------------------------------------------------------- #
# sequence
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("dark, position", [(1, 67), (2, 339)])
def test_tdx_deletes_exactly_one_chromophore_tyrosine(tmp_path, dark, position):
    subprocess.run([sys.executable, str(GENERATOR), "--dead-barrel", str(dark),
                    "--outdir", str(tmp_path)], check=True, cwd=ROOT,
                   capture_output=True, text=True)
    td = _fasta(PARENT)
    tdx = _fasta(tmp_path / f"pdVenus-TDX_derived_510aa_dark{dark}.fasta")

    assert len(td) == 511 and len(tdx) == 510
    # deleting residue `position` and nothing else
    assert tdx == td[:position - 1] + td[position:]
    assert td[position - 1] == "Y"
    assert td[position - 2:position + 1] == "GYG"
    # the dark barrel is left with Gly-Gly, and the other triad is untouched
    assert tdx[position - 2:position] == "GG"
    survivor = 67 if dark == 2 else 338
    assert [m.start() + 2 for m in re.finditer("GYG", tdx)] == [survivor]


def test_dark_barrel_2_reproduces_the_addgene_deposit(tmp_path):
    """The deposited construct darkens the C-terminal barrel; this must BE it."""
    subprocess.run([sys.executable, str(GENERATOR), "--outdir", str(tmp_path)],
                   check=True, cwd=ROOT, capture_output=True, text=True)
    assert _fasta(tmp_path / "pdVenus-TDX_derived_510aa_dark2.fasta") == _fasta(DEPOSIT)


def test_substitute_mode_is_available_but_is_not_the_default(tmp_path):
    """The Y->C reading is a different protein (Amber); keep it reproducible and named."""
    subprocess.run([sys.executable, str(GENERATOR), "--mode", "substitute",
                    "--outdir", str(tmp_path)], check=True, cwd=ROOT,
                   capture_output=True, text=True)
    amber = _fasta(tmp_path / "pdVenus-TDX_derived_511aa_dark2_amberY67C.fasta")
    td = _fasta(PARENT)
    assert len(amber) == 511
    assert amber == td[:338] + "C" + td[339:]
    assert amber != _fasta(DEPOSIT)


def test_alphafold_job_request_is_submittable(tmp_path):
    subprocess.run([sys.executable, str(GENERATOR), "--outdir", str(tmp_path)],
                   check=True, cwd=ROOT, capture_output=True, text=True)
    job = json.loads((tmp_path / "af3_tdx_job_request_dark2.json").read_text())
    assert isinstance(job, list) and len(job) == 1
    assert job[0]["dialect"] == "alphafoldserver" and job[0]["version"] == 3
    chain = job[0]["sequences"][0]["proteinChain"]
    assert chain["sequence"] == _fasta(tmp_path / "pdVenus-TDX_derived_510aa_dark2.fasta")
    assert chain["count"] == 1
    # same seed as the TD run in af3_server_run/, so the two are comparable
    td_job = json.loads((ROOT / "alphafold_tandem" / "af3_server_run" /
                         "fold_dvenus_td_tandem_job_request.json").read_text())
    assert job[0]["modelSeeds"] == td_job[0]["modelSeeds"]


# --------------------------------------------------------------------------- #
# structure
# --------------------------------------------------------------------------- #

# Mature CR2 -> the Gly-Tyr-Gly it formed from, for writing a surviving
# chromophore back into a synthetic AlphaFold model of the unreacted chain.
LIVE_TRIAD_ATOMS = (
    ("GLY", (("N", "N1"), ("CA", "CA1"), ("C", "C1"))),
    ("TYR", (("N", "N2"), ("CA", "CA2"), ("CB", "CB2"), ("CG", "CG2"),
             ("CD1", "CD1"), ("CD2", "CD2"), ("CE1", "CE1"), ("CE2", "CE2"),
             ("CZ", "CZ"), ("OH", "OH"), ("C", "C2"), ("O", "O2"))),
    ("GLY", (("N", "N3"), ("CA", "CA3"), ("C", "C3"), ("O", "O3"))),
)


def _synthetic_af_model(path: Path, dark: int = 2, middle: str | None = None) -> None:
    """An AF3-format model of the TDX sequence, made from the v4 structure itself.

    Non-chromophore residues are renumbered into the 510-aa numbering AF uses;
    the dead barrel's pair is written as the Gly-Gly the in-place collapse
    produces.  Grafting from this must therefore be the identity.

    `middle` inserts an extra residue of that name between the two glycines,
    which is what a model of TD (TYR) or of the superseded Y->C reading (CYS)
    looks like, and must be rejected.
    """
    residues = build.read_pdb(V4)
    cr2 = [i for i, residue in enumerate(residues) if residue["name"] == "CR2"]
    dead_index = cr2[dark - 1]
    to_af = build.sequence_numbering(residues, dead_index)

    rows = []
    for index, residue in enumerate(residues):
        if index == dead_index:
            start = to_af[index]
            pair = build.collapse_in_place(residue)
            if middle:
                pair = [pair[0],
                        (middle, [("N", build.atom_map(residue)["N2"]),
                                  ("CA", build.atom_map(residue)["CA2"]),
                                  ("C", build.atom_map(residue)["C2"])]),
                        pair[1]]
            for offset, (name, atoms) in enumerate(pair):
                for atom, xyz in atoms:
                    rows.append((atom, name, start + offset, xyz))
            continue
        if residue["name"] == "CR2":
            # the surviving chromophore, written back out as the three residues
            # it formed from so the synthetic model has the same chain length
            # AlphaFold would produce for this construct
            xyz = build.atom_map(residue)
            start = to_af[index]
            for offset, (name, pairs) in enumerate(LIVE_TRIAD_ATOMS):
                for target, source in pairs:
                    rows.append((target, name, start + offset, xyz[source]))
            continue
        number = to_af[index]
        for atom in residue["atoms"]:
            rows.append((atom["name"], residue["name"], number, atom["xyz"]))

    header = ["data_synthetic", "loop_"] + [
        f"_atom_site.{field}" for field in
        ("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
         "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
         "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
         "B_iso_or_equiv", "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num")]
    lines = list(header)
    for serial, (atom, name, number, xyz) in enumerate(rows, start=1):
        lines.append(f"ATOM {serial} {atom[0]} {atom} . {name} A 1 {number} ? "
                     f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} 1.00 90.00 {number} A 1")
    path.write_text("\n".join(lines) + "\n")


def _run(args, **kwargs):
    return subprocess.run([sys.executable, str(BUILDER), *args],
                          cwd=ROOT, capture_output=True, text=True, **kwargs)


def test_in_place_collapse_leaves_a_real_peptide_bond(tmp_path):
    """Deleting the tyrosine turns the imidazolinone C1-N3 ring bond into the
    Gly-Gly peptide bond, and frees the lobe the carbonyl oxygen needs.  That is
    why this route is merely unrelaxed here, where the Y->C reading of the
    construct made it geometrically impossible."""
    output = tmp_path / "inplace.pdb"
    result = _run(["--output", str(output)])
    assert result.returncode == 0, result.stderr
    assert "CLASH" not in result.stdout
    bond = float(re.search(r"Gly\(1\) C - Gly\(2\) N: ([\d.]+) A", result.stdout).group(1))
    assert 1.20 <= bond <= 1.50
    # nothing moves, so the junctions stay exactly as the crystal had them
    for label in ("N-terminal", "C-terminal"):
        distance = float(re.search(rf"{label} junction C-N: ([\d.]+) A",
                                   result.stdout).group(1))
        assert distance < 1.45


def test_alphafold_graft_reproduces_a_known_pair(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    result = _run(["--from-af3", str(model), "--output", str(output)])
    assert result.returncode == 0, result.stderr

    rmsd = float(re.search(r"superposition RMSD ([\d.]+) A", result.stdout).group(1))
    assert rmsd < 0.01

    grafted = {(residue["number"], atom["name"]): atom["xyz"]
               for residue in build.read_pdb(output)
               for atom in residue["atoms"] if 336 <= residue["number"] <= 337}
    expected = build.collapse_in_place(
        next(r for r in build.read_pdb(V4) if r["name"] == "CR2" and r["number"] == 336))
    for offset, (_, atoms) in enumerate(expected):
        for atom, xyz in atoms:
            assert grafted[(336 + offset, atom)] == pytest.approx(xyz, abs=2e-3)


def test_the_surviving_chromophore_and_the_register_are_untouched(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    assert _run(["--from-af3", str(model), "--output", str(output)]).returncode == 0

    before = {r["number"]: r for r in build.read_pdb(V4)}
    after = {r["number"]: r for r in build.read_pdb(output)}

    live = [r for r in after.values() if r["name"] == "CR2"]
    assert len(live) == 1 and live[0]["number"] == 66
    for atom, reference in zip(live[0]["atoms"], before[66]["atoms"]):
        assert atom["name"] == reference["name"]
        assert atom["xyz"] == pytest.approx(reference["xyz"], abs=1e-3)

    # 507 residues become 508: one CR2 stands for two, not three.
    assert len(after) == len(before) + 1
    # Nothing outside the replaced pair moves at all -- this is the register.
    for number in (1, 65, 100, 335):
        for atom, reference in zip(after[number]["atoms"], before[number]["atoms"]):
            assert atom["xyz"] == pytest.approx(reference["xyz"], abs=1e-3)


def test_dead_pair_is_two_glycines_with_no_tyrosine_left(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    assert _run(["--from-af3", str(model), "--output", str(output)]).returncode == 0
    residues = {r["number"]: r for r in build.read_pdb(output)}
    assert [residues[n]["name"] for n in (336, 337)] == ["GLY", "GLY"]
    atoms = {atom["name"] for n in (336, 337) for atom in residues[n]["atoms"]}
    # no side chain at all: no cysteine sulfur, no phenol, not even a CB
    assert not atoms & {"SG", "CB", "CG", "CZ", "OH", "CD1", "CE1"}


@pytest.mark.parametrize("middle, label", [("TYR", "TD"), ("CYS", "the Y->C reading")])
def test_graft_rejects_a_model_of_the_wrong_construct(tmp_path, middle, label):
    """A residue between the glycines means the model is not of the deposit."""
    model = tmp_path / f"wrong_{middle}.cif"
    _synthetic_af_model(model, middle=middle)
    result = _run(["--from-af3", str(model), "--output", str(tmp_path / "x.pdb")])
    assert result.returncode != 0
    assert "expected ['GLY', 'GLY']" in result.stderr


def test_graft_rejects_the_real_td_model():
    """The TD model has a tyrosine there and a 511-aa chain; both must be caught."""
    td_model = (ROOT / "alphafold_tandem" / "af3_server_run" /
                "fold_dvenus_td_tandem_model_0.cif")
    result = _run(["--from-af3", str(td_model)])
    assert result.returncode != 0
    assert "expected ['GLY', 'GLY']" in result.stderr or "residues, expected" in result.stderr


def test_the_real_build_keeps_the_register_and_makes_no_clashes(tmp_path):
    """The production build, from the AlphaFold Server run on this branch."""
    model = (ROOT / "alphafold_tandem" / "af3_tdx_server_run" /
             "fold_dvenus_tdx_dark2_model_0.cif")
    if not model.exists():
        pytest.skip("af3_tdx_server_run/ not populated")
    output = tmp_path / "tdx.pdb"
    result = _run(["--from-af3", str(model), "--output", str(output)])
    assert result.returncode == 0, result.stderr
    assert "CLASH" not in result.stdout

    rmsd = float(re.search(r"superposition RMSD ([\d.]+) A", result.stdout).group(1))
    assert rmsd < 0.6
    assert "surviving CR2: 1 residue, max coordinate change 0.000 A" in result.stdout
    assert "all other residues: max coordinate change 0.000 A" in result.stdout
    assert len(build.read_pdb(output)) == 508


# --------------------------------------------------------------------------- #
# the MD driver
# --------------------------------------------------------------------------- #

def test_run_nvt_accepts_a_single_chromophore():
    """TDX has one CR2; run_nvt.py used to demand exactly two."""
    source = (ROOT / "run_nvt.py").read_text()
    assert "needs exactly two CR2 residues" not in source
    assert source.count('if not target_cr2:\n        raise RuntimeError("Target topology has no CR2 residue")') == 2


def test_tdx_md_numbers_close():
    """1 fs x 2000 steps x 500 frames = 1 ns, and the script must say so."""
    script = (ROOT / "tandem_tdx_build" / "run_md_tdx.sh").read_text()
    dt = float(re.search(r"^DT=\$\{DT:-([\d.]+)\}", script, re.M).group(1))
    per_frame = int(re.search(r"^STEPS_PER_FRAME=\$\{STEPS_PER_FRAME:-(\d+)\}", script, re.M).group(1))
    frames = int(re.search(r"^FRAMES=\$\{FRAMES:-(\d+)\}", script, re.M).group(1))
    assert (dt, per_frame, frames) == (1.0, 2000, 500)
    assert dt * per_frame * frames / 1e6 == pytest.approx(1.0)   # nanoseconds

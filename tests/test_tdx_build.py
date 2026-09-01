"""Guard the dVenus-TDX construct: one residue in the sequence, one triad in the structure.

TDX is TD with one barrel's chromophore extinguished (Amber, barrel-local
Y67C).  Two things can silently go wrong and both would be invisible in a
finished trajectory:

  * the sequence drifting from the Addgene #250938 parent by more than the one
    intended residue, which would make a TD-vs-TDX difference uninterpretable;
  * the dead barrel's triad being built by un-maturing the chromophore in
    place, which cannot work -- the mature C1 already carries three
    substituents (CA1, N2 and the ring N3), so the lobe the restored carbonyl
    oxygen needs is the one N3 occupies, and the rebuilt O lands 0.21 A from
    it.  The builder must refuse that structure rather than hand a 0.21 A
    overlap to a minimiser.

The structural tests drive the AlphaFold graft path against a synthetic model
built from the v4 structure itself, so the correct answer is known exactly: the
grafted triad must come back where it started.
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
def test_tdx_differs_from_the_deposit_by_one_chromophore_tyrosine(tmp_path, dark, position):
    subprocess.run([sys.executable, str(GENERATOR), "--dead-barrel", str(dark),
                    "--outdir", str(tmp_path)], check=True, cwd=ROOT,
                   capture_output=True, text=True)
    td = _fasta(PARENT)
    tdx = _fasta(tmp_path / f"pdVenus-TDX_derived_511aa_dark{dark}.fasta")

    assert len(tdx) == len(td) == 511
    differences = [i + 1 for i, (a, b) in enumerate(zip(td, tdx)) if a != b]
    assert differences == [position]
    assert td[position - 1] == "Y" and tdx[position - 1] == "C"
    # the mutated residue is the middle of a chromophore-forming triad ...
    assert td[position - 2:position + 1] == "GYG"
    # ... and the other triad is left alone.
    assert [m.start() + 2 for m in re.finditer("GYG", tdx)] == [339 if dark == 1 else 67]


def test_alphafold_job_request_is_submittable(tmp_path):
    subprocess.run([sys.executable, str(GENERATOR), "--outdir", str(tmp_path)],
                   check=True, cwd=ROOT, capture_output=True, text=True)
    job = json.loads((tmp_path / "af3_tdx_job_request_dark2.json").read_text())
    assert isinstance(job, list) and len(job) == 1
    assert job[0]["dialect"] == "alphafoldserver" and job[0]["version"] == 3
    chain = job[0]["sequences"][0]["proteinChain"]
    assert chain["sequence"] == _fasta(tmp_path / "pdVenus-TDX_derived_511aa_dark2.fasta")
    assert chain["count"] == 1
    # same seed as the TD run in af3_server_run/, so the two are comparable
    td_job = json.loads((ROOT / "alphafold_tandem" / "af3_server_run" /
                         "fold_dvenus_td_tandem_job_request.json").read_text())
    assert job[0]["modelSeeds"] == td_job[0]["modelSeeds"]


# --------------------------------------------------------------------------- #
# structure
# --------------------------------------------------------------------------- #

def _synthetic_af_model(path: Path, dark: int = 2) -> None:
    """An AF3-format model of the TDX sequence, made from the v4 structure itself.

    Non-chromophore residues are renumbered into the 511-aa numbering AF uses;
    the dead barrel's triad is written as the Gly-Cys-Gly the in-place mapping
    produces.  Grafting from this must therefore be the identity.
    """
    residues = build.read_pdb(V4)
    cr2 = [i for i, residue in enumerate(residues) if residue["name"] == "CR2"]
    dead_index = cr2[dark - 1]

    def to_af(number: int) -> int | None:
        if number <= 65:
            return number
        if number == 66:
            return None
        if number <= 335:
            return number + 2
        if number == 336:
            return None
        return number + 4

    rows = []
    for index, residue in enumerate(residues):
        if index == dead_index:
            start = to_af(residue["number"] - 1) + 1
            for offset, (name, atoms) in enumerate(build.unmature_in_place(residue)):
                for atom, xyz in atoms:
                    rows.append((atom, name, start + offset, xyz))
            continue
        if residue["name"] == "CR2":
            continue
        number = to_af(residue["number"])
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


def test_in_place_unmaturation_is_refused():
    """The carbonyl oxygen has nowhere to go, and the builder must say so."""
    result = _run([])
    assert result.returncode != 0
    assert "refusing to write" in result.stderr
    assert "CLASH 0.21 A: GLY336 O ... GLY338 N" in result.stdout


def test_alphafold_graft_reproduces_a_known_triad(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    result = _run(["--from-af3", str(model), "--output", str(output), "--force-strained"])
    assert result.returncode == 0, result.stderr

    rmsd = float(re.search(r"superposition RMSD ([\d.]+) A", result.stdout).group(1))
    assert rmsd < 0.01

    grafted = {(residue["number"], atom["name"]): atom["xyz"]
               for residue in build.read_pdb(output)
               for atom in residue["atoms"] if 336 <= residue["number"] <= 338}
    expected = build.unmature_in_place(
        next(r for r in build.read_pdb(V4) if r["name"] == "CR2" and r["number"] == 336))
    for offset, (_, atoms) in enumerate(expected):
        for atom, xyz in atoms:
            assert grafted[(336 + offset, atom)] == pytest.approx(xyz, abs=2e-3)


def test_the_surviving_chromophore_and_the_register_are_untouched(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    assert _run(["--from-af3", str(model), "--output", str(output),
                 "--force-strained"]).returncode == 0

    before = {r["number"]: r for r in build.read_pdb(V4)}
    after = {r["number"]: r for r in build.read_pdb(output)}

    live = [r for r in after.values() if r["name"] == "CR2"]
    assert len(live) == 1 and live[0]["number"] == 66
    for atom, reference in zip(live[0]["atoms"], before[66]["atoms"]):
        assert atom["name"] == reference["name"]
        assert atom["xyz"] == pytest.approx(reference["xyz"], abs=1e-3)

    # 507 residues become 509: one CR2 stands for three.
    assert len(after) == len(before) + 2
    # Nothing upstream of the dead triad moves at all.
    for number in (1, 65, 100, 335):
        for atom, reference in zip(after[number]["atoms"], before[number]["atoms"]):
            assert atom["xyz"] == pytest.approx(reference["xyz"], abs=1e-3)


def test_dead_triad_is_glycine_cysteine_glycine(tmp_path):
    model = tmp_path / "synthetic_tdx_model.cif"
    _synthetic_af_model(model)
    output = tmp_path / "tdx.pdb"
    assert _run(["--from-af3", str(model), "--output", str(output),
                 "--force-strained"]).returncode == 0
    residues = {r["number"]: r for r in build.read_pdb(output)}
    assert [residues[n]["name"] for n in (336, 337, 338)] == ["GLY", "CYS", "GLY"]
    assert "SG" in {atom["name"] for atom in residues[337]["atoms"]}
    # the phenol is gone
    assert not {"CZ", "OH"} & {atom["name"] for atom in residues[337]["atoms"]}


def test_graft_rejects_a_model_of_the_wrong_construct(tmp_path):
    """The TD model has a tyrosine there; grafting it would rebuild the chromophore."""
    td_model = ROOT / "alphafold_tandem" / "af3_server_run" / "fold_dvenus_td_tandem_model_0.cif"
    result = _run(["--from-af3", str(td_model), "--output", str(tmp_path / "x.pdb")])
    assert result.returncode != 0
    assert "expected GLY CYS GLY" in result.stderr


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

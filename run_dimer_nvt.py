#!/usr/bin/env python3
"""
Run OpenMM NVT relaxation for the dimer structure and render a video.
Standalone: does not call any other python scripts.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import shlex
import shutil
import subprocess
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def discover_openmm_plugin_dirs() -> List[Path]:
    candidates: List[Path] = []
    seen = set()

    def add(path: Optional[Path | str]) -> None:
        if path is None:
            return
        p = Path(path).expanduser()
        if not p.exists() or not p.is_dir():
            return
        resolved = p.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        candidates.append(resolved)

    env_plugin_dir = os.environ.get("OPENMM_PLUGIN_DIR")
    if env_plugin_dir:
        add(env_plugin_dir)

    exe_prefix = Path(sys.executable).resolve().parents[1]
    add(exe_prefix / "lib" / "plugins")

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        add(Path(conda_prefix) / "lib" / "plugins")

    openmm_spec = importlib.util.find_spec("openmm")
    if openmm_spec and openmm_spec.origin:
        origin = Path(openmm_spec.origin).resolve()
        for parent in origin.parents:
            add(parent / "lib" / "plugins")

    return candidates


def score_openmm_plugin_dir(plugin_dir: Path) -> int:
    if not plugin_dir.exists() or not plugin_dir.is_dir():
        return -1
    names = {entry.name for entry in plugin_dir.iterdir() if entry.is_file()}
    score = 0
    if "libOpenMMCUDA.so" in names:
        score += 100
    if "libOpenMMCPU.so" in names:
        score += 20
    if any("Reference" in name for name in names):
        score += 5
    score += min(len(names), 50)
    return score


def choose_best_openmm_plugin_dir() -> Optional[Path]:
    candidates = discover_openmm_plugin_dirs()
    if not candidates:
        return None
    ranked = sorted(candidates, key=score_openmm_plugin_dir, reverse=True)
    return ranked[0]


def configure_openmm_env() -> None:
    plugin_dir = choose_best_openmm_plugin_dir()
    if plugin_dir is None:
        return
    lib_dir = plugin_dir.parent
    os.environ["OPENMM_PLUGIN_DIR"] = str(plugin_dir)
    os.environ["OPENMM_LIB_PATH"] = str(lib_dir)


configure_openmm_env()

import pdbfixer
import openmm
from openmm import Platform, unit
from openmm.app import (
    CutoffNonPeriodic,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    PDBReporter,
    PME,
    Simulation,
)

WATER_RESIDUE_NAMES = {"HOH", "WAT", "SOL"}
STANDARD_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "LYN",
    "CYM",
    "CYX",
    "ACE",
    "NME",
    "NHE",
}
COMMON_ION_RESIDUES = {
    "NA",
    "K",
    "CL",
    "CA",
    "MG",
    "ZN",
    "MN",
    "FE",
    "CU",
    "NI",
    "LI",
    "CS",
    "RB",
    "I",
}
LJ_BY_ELEMENT = {
    "H": (0.250, 0.0157),
    "C": (0.340, 0.2761),
    "N": (0.325, 0.1700),
    "O": (0.296, 0.2100),
    "S": (0.356, 1.0460),
    "P": (0.374, 0.8368),
}
ATOMIC_MASS_DALTON = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
}


def list_platform_names() -> List[str]:
    return [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]


def load_additional_openmm_plugins() -> None:
    for plugin_dir in discover_openmm_plugin_dirs():
        try:
            Platform.loadPluginsFromDirectory(str(plugin_dir))
        except Exception:
            continue


def pick_platform(preferred_name: str = "CUDA", strict: bool = False) -> Platform:
    if preferred_name:
        try:
            return Platform.getPlatformByName(preferred_name)
        except Exception as exc:
            load_additional_openmm_plugins()
            try:
                return Platform.getPlatformByName(preferred_name)
            except Exception:
                if strict:
                    available = ", ".join(list_platform_names()) or "none"
                    raise RuntimeError(
                        f"Requested OpenMM platform '{preferred_name}' is not available. Available: {available}"
                    ) from exc
            if preferred_name.upper() == "CUDA":
                for fallback_name in ("OpenCL", "CPU", "Reference"):
                    try:
                        fallback = Platform.getPlatformByName(fallback_name)
                        print(
                            f"    - Warning: requested OpenMM platform '{preferred_name}' unavailable; "
                            f"falling back to '{fallback_name}'"
                        )
                        return fallback
                    except Exception:
                        continue
            available = ", ".join(list_platform_names()) or "none"
            raise RuntimeError(
                f"Requested OpenMM platform '{preferred_name}' is not available. Available: {available}"
            ) from exc
    for name in ("CUDA", "OpenCL", "CPU", "Reference"):
        try:
            return Platform.getPlatformByName(name)
        except Exception:
            continue
    if Platform.getNumPlatforms() > 0:
        return Platform.getPlatform(0)
    raise RuntimeError("No OpenMM platform is available")


def safe_remove_directory(path: Path, workspace_root: Path) -> None:
    path = Path(path)
    workspace_root = Path(workspace_root).resolve()
    if path.is_symlink():
        raise RuntimeError(f"Refusing to delete symlinked workdir: {path}")
    resolved = path.resolve()
    if not resolved.is_dir():
        raise RuntimeError(f"Refusing to delete non-directory workdir path: {resolved}")
    forbidden = {Path("/"), Path.home().resolve(), workspace_root}
    if resolved in forbidden:
        raise RuntimeError(f"Refusing to delete unsafe workdir path: {resolved}")
    if workspace_root not in resolved.parents:
        raise RuntimeError(f"Refusing to delete workdir outside workspace: {resolved}")
    shutil.rmtree(resolved)


def get_periodic_box_lengths_ang(topology) -> Optional[np.ndarray]:
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return None
    lengths = np.array([np.linalg.norm(v.value_in_unit(unit.angstrom)) for v in vectors], dtype=float)
    if np.any(lengths <= 1.0e-8):
        return None
    return lengths


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def element_class(symbol: str) -> str:
    mapping = {"H": "H", "C": "C", "N": "N", "O": "O", "S": "S", "P": "P"}
    return mapping.get(symbol, safe_name(symbol))


def atom_mass_dalton(atom, symbol: str) -> float:
    if atom.element is not None and atom.element.mass is not None:
        return float(atom.element.mass.value_in_unit(unit.dalton))
    return float(ATOMIC_MASS_DALTON.get(symbol, ATOMIC_MASS_DALTON["C"]))


def find_nonstandard_residues(topology) -> List[str]:
    names = set()
    for residue in topology.residues():
        rname = residue.name
        if rname in STANDARD_PROTEIN_RESIDUES:
            continue
        if rname in WATER_RESIDUE_NAMES:
            continue
        if rname in COMMON_ION_RESIDUES:
            continue
        names.add(rname)
    return sorted(names)


def residue_names_in_forcefield_xml(xml_path: Path) -> set:
    text = Path(xml_path).read_text(errors="replace")
    return set(re.findall(r"<Residue\s+name=\"([^\"]+)\"", text))


def write_generic_forcefield_xml(topology, positions, residue_names: List[str], xml_path: Path) -> None:
    residue_names = sorted(set(residue_names))
    residues_by_name: Dict[str, object] = {}
    for residue in topology.residues():
        residues_by_name.setdefault(residue.name, residue)

    coords_nm = np.array([p.value_in_unit(unit.nanometer) for p in positions])
    atom_type_data = []
    residue_blocks = []
    bond_params = []
    angle_params = []
    nonbonded_params = []

    for residue_name in residue_names:
        residue = residues_by_name.get(residue_name)
        if residue is None:
            continue

        atoms = list(residue.atoms())
        atom_name_to_atom = {a.name: a for a in atoms}
        type_by_atom_name = {}

        for atom in atoms:
            symbol = atom.element.symbol if atom.element is not None else "C"
            type_name = safe_name(f"{residue_name}_{atom.name}")
            class_name = element_class(symbol)
            type_by_atom_name[atom.name] = type_name
            atom_type_data.append((type_name, class_name, symbol, atom_mass_dalton(atom, symbol)))
            sigma, epsilon = LJ_BY_ELEMENT.get(symbol, (0.340, 0.2000))
            nonbonded_params.append((type_name, 0.0, sigma, epsilon))

        bonds_local = []
        neighbors = {atom.name: set() for atom in atoms}
        external_bond_atoms = set()

        for bond in topology.bonds():
            a1 = bond.atom1
            a2 = bond.atom2
            in1 = a1.residue == residue
            in2 = a2.residue == residue
            if in1 and in2:
                bonds_local.append((a1, a2))
                neighbors[a1.name].add(a2.name)
                neighbors[a2.name].add(a1.name)
                distance_nm = np.linalg.norm(coords_nm[a1.index] - coords_nm[a2.index])
                bond_params.append((type_by_atom_name[a1.name], type_by_atom_name[a2.name], distance_nm, 300000.0))
            elif in1 ^ in2:
                external_bond_atoms.add(a1.name if in1 else a2.name)

        for center_name, neighbor_names in neighbors.items():
            if len(neighbor_names) < 2:
                continue
            center_atom = atom_name_to_atom[center_name]
            for name_i, name_k in combinations(sorted(neighbor_names), 2):
                atom_i = atom_name_to_atom[name_i]
                atom_k = atom_name_to_atom[name_k]
                vec_i = coords_nm[atom_i.index] - coords_nm[center_atom.index]
                vec_k = coords_nm[atom_k.index] - coords_nm[center_atom.index]
                norm_i = np.linalg.norm(vec_i)
                norm_k = np.linalg.norm(vec_k)
                if norm_i < 1e-8 or norm_k < 1e-8:
                    continue
                cosang = np.dot(vec_i, vec_k) / (norm_i * norm_k)
                cosang = float(np.clip(cosang, -1.0, 1.0))
                theta = float(np.arccos(cosang))
                angle_params.append(
                    (
                        type_by_atom_name[name_i],
                        type_by_atom_name[center_name],
                        type_by_atom_name[name_k],
                        theta,
                        300.0,
                    )
                )

        residue_blocks.append((residue_name, atoms, bonds_local, external_bond_atoms, type_by_atom_name))

    with open(xml_path, "w", encoding="utf-8") as handle:
        handle.write("<ForceField>\n")
        handle.write("  <AtomTypes>\n")
        for type_name, class_name, symbol, mass in atom_type_data:
            handle.write(
                f"    <Type name=\"{type_name}\" class=\"{class_name}\" element=\"{symbol}\" mass=\"{mass:.6f}\"/>\n"
            )
        handle.write("  </AtomTypes>\n")

        handle.write("  <Residues>\n")
        for residue_name, atoms, bonds_local, external_bond_atoms, type_by_atom_name in residue_blocks:
            handle.write(f"    <Residue name=\"{residue_name}\">\n")
            for atom in atoms:
                handle.write(f"      <Atom name=\"{atom.name}\" type=\"{type_by_atom_name[atom.name]}\"/>\n")
            for atom1, atom2 in bonds_local:
                handle.write(f"      <Bond atomName1=\"{atom1.name}\" atomName2=\"{atom2.name}\"/>\n")
            for atom_name in sorted(external_bond_atoms):
                handle.write(f"      <ExternalBond atomName=\"{atom_name}\"/>\n")
            handle.write("    </Residue>\n")
        handle.write("  </Residues>\n")

        handle.write("  <HarmonicBondForce>\n")
        for type1, type2, length, k_value in bond_params:
            handle.write(
                f"    <Bond type1=\"{type1}\" type2=\"{type2}\" length=\"{length:.6f}\" k=\"{k_value:.1f}\"/>\n"
            )
        handle.write("  </HarmonicBondForce>\n")

        handle.write("  <HarmonicAngleForce>\n")
        for type1, type2, type3, angle, k_value in angle_params:
            handle.write(
                f"    <Angle type1=\"{type1}\" type2=\"{type2}\" type3=\"{type3}\" angle=\"{angle:.6f}\" k=\"{k_value:.1f}\"/>\n"
            )
        handle.write("  </HarmonicAngleForce>\n")

        handle.write("  <NonbondedForce coulomb14scale=\"0.833333\" lj14scale=\"0.5\">\n")
        for type_name, charge, sigma, epsilon in nonbonded_params:
            handle.write(
                f"    <Atom type=\"{type_name}\" charge=\"{charge:.6f}\" sigma=\"{sigma:.6f}\" epsilon=\"{epsilon:.6f}\"/>\n"
            )
        handle.write("  </NonbondedForce>\n")
        handle.write("</ForceField>\n")


def choose_cutoff_from_box(topology, default_nm: float = 1.0):
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return default_nm * unit.nanometer
    lengths_nm = [np.linalg.norm(v.value_in_unit(unit.nanometer)) for v in vectors]
    max_allowed = 0.49 * min(lengths_nm)
    cutoff_nm = min(default_nm, max_allowed)
    cutoff_nm = max(cutoff_nm, 0.1)
    return cutoff_nm * unit.nanometer


def resolve_path(path: Path) -> Path:
    path = path.expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    return path


def preprocess_pdb_for_solvation(pdb_path: Path, workdir: Path) -> Path:
    """Insert TER records at residue gaps or protein/nonstandard boundaries to help solvation."""
    output_path = workdir / "pdb_for_solvation.pdb"
    inserted = False

    prev_chain = None
    prev_resi = None
    prev_resn = None
    prev_is_protein = None

    def is_protein_residue(resn: str) -> bool:
        return resn in STANDARD_PROTEIN_RESIDUES

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle, output_path.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("TER"):
                out_handle.write(line + "\n")
                prev_chain = prev_resi = prev_resn = prev_is_protein = None
                continue

            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip()
                resi = line[22:26].strip()
                resn = line[17:20].strip().upper()
                is_protein = is_protein_residue(resn)

                if prev_chain is not None:
                    new_residue = (chain, resi, resn) != (prev_chain, prev_resi, prev_resn)
                    if new_residue and chain == prev_chain:
                        gap = False
                        if prev_resi and resi and prev_resi.isdigit() and resi.isdigit():
                            gap = int(resi) != int(prev_resi) + 1
                        boundary = (prev_is_protein != is_protein)
                        if gap or boundary:
                            out_handle.write("TER\n")
                            inserted = True

                if (chain, resi, resn) != (prev_chain, prev_resi, prev_resn):
                    prev_chain, prev_resi, prev_resn, prev_is_protein = chain, resi, resn, is_protein

            out_handle.write(line + "\n")

    if not inserted:
        try:
            output_path.unlink()
        except Exception:
            pass
        return pdb_path
    return output_path


def count_pdb_models(pdb_path: Path) -> int:
    model_count = 0
    atom_seen = False
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("MODEL"):
                model_count += 1
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                atom_seen = True
    if model_count > 0:
        return model_count
    return 1 if atom_seen else 0


def write_pymol_script(
    script_path: Path,
    traj_path: Path,
    frame_prefix: Path,
    width: int,
    height: int,
    align_states: bool,
    zoom_target: str,
    zoom_buffer_a: float,
    ray: bool,
    qm_residue: str,
    qm_protein_cutoff_a: float,
    qm_nearest_waters: int,
) -> None:
    ray_flag = "1" if ray else "0"
    ray_int = 1 if ray else 0
    align_states_flag = "True" if align_states else "False"
    script_text = f"""reinitialize
set bg_rgb, [1, 1, 1]
set grid_mode, 0
set transparency, 0.5
set cartoon_transparency, 0.6
set cartoon_fancy_helices, 1
set line_width, 1.0
set depth_cue, 0
set auto_zoom, 0
set ray_opaque_background, off
set stick_radius, 0.15
set sphere_scale, 0.20
set line_smooth, 1
set antialias, 4
set ray_trace_mode, 1
set ray_trace_frames, {ray_flag}
viewport {width}, {height}
load {traj_path.as_posix()}, qm_min
hide everything, qm_min
python
from pymol import cmd, util
import math

WATER_NAMES = {{"HOH", "WAT", "SOL"}}
qm_residue_name = "{qm_residue.upper()}"
protein_cutoff = float({qm_protein_cutoff_a})
nearest_waters = max(int({qm_nearest_waters}), 0)
zoom_target_mode = "{zoom_target}"
zoom_buffer = float({zoom_buffer_a})
align_states = {align_states_flag}

def min_distance(points_a, points_b):
    best = 1.0e9
    for ax, ay, az in points_a:
        for bx, by, bz in points_b:
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < best:
                best = dist
    return best

model = cmd.get_model("qm_min", state=1)
if not model.atom:
    raise RuntimeError("Trajectory loaded but no atoms found.")

# Group atom coordinates by residue key.
residue_atoms = {{}}
core_points = []
for atom in model.atom:
    chain = (atom.chain or "").strip()
    resi = str(atom.resi).strip()
    resn = (atom.resn or "").strip().upper()
    key = (chain, resi, resn)
    coord = (float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]))
    residue_atoms.setdefault(key, []).append((atom.name.strip().upper(), atom.symbol.strip().upper(), coord))
    if resn == qm_residue_name:
        core_points.append(coord)

qm_core_keys = {{key for key in residue_atoms if key[2] == qm_residue_name}}
qm_cutoff_keys = set()
water_candidates = []

if core_points:
    for key, atoms in residue_atoms.items():
        if key in qm_core_keys:
            continue
        resn = key[2]
        atom_points = [entry[2] for entry in atoms]
        if resn in WATER_NAMES:
            oxygen_points = [entry[2] for entry in atoms if entry[1] == "O" or entry[0].startswith("O")]
            probe_points = oxygen_points if oxygen_points else atom_points
            water_candidates.append((min_distance(probe_points, core_points), key))
        else:
            if min_distance(atom_points, core_points) <= protein_cutoff:
                qm_cutoff_keys.add(key)

water_candidates.sort(key=lambda item: item[0])
qm_water_keys = []
for dist, key in water_candidates:
    if key in qm_core_keys or key in qm_cutoff_keys:
        continue
    qm_water_keys.append((dist, key))
    if len(qm_water_keys) >= nearest_waters:
        break

def select_from_keys(name, keys):
    keys = list(keys)
    if not keys:
        cmd.select(name, "none")
        return
    clauses = []
    for key in keys:
        chain, resi, resn = key
        resi = str(resi).strip()
        clause = f"(resn {{resn}} and resi {{resi}})"
        if chain:
            clause = f"({{clause}} and chain {{chain}})"
        clauses.append(clause)

    chunk_size = 50
    for i in range(0, len(clauses), chunk_size):
        chunk = clauses[i:i+chunk_size]
        sel_str = " or ".join(chunk)
        if i == 0:
            cmd.select(name, f"qm_min and ({{sel_str}})")
        else:
            cmd.select(name, f"{{name}} or ({{sel_str}})")

cmd.select("qm_core", "none")
select_from_keys("qm_core", sorted(qm_core_keys))
select_from_keys("qm_cutoff", sorted(qm_cutoff_keys))
select_from_keys("qm_waters", [key for _, key in qm_water_keys])
# Do not include solvent in the displayed QM region.
cmd.select("qm_region", "qm_core or qm_cutoff")
cmd.select("solvent_all", "qm_min and resn HOH+WAT+SOL")

# Match visualise_dimer.pml logic: non-QM protein/ligand as lines+spectrum
cmd.select("non_qm_context", "qm_min and not qm_region and not resn HOH+WAT+SOL")
cmd.select("context_structure", "qm_region or non_qm_context")
if zoom_target_mode == "qm":
    cmd.select("zoom_target", "qm_region")
elif zoom_target_mode == "context":
    cmd.select("zoom_target", "context_structure")
else:
    cmd.select("zoom_target", "qm_min and ((polymer.protein and not resn HOH+WAT+SOL) or qm_region)")
if cmd.count_atoms("zoom_target") < 1:
    cmd.select("zoom_target", "context_structure")
    zoom_target_mode = "context-fallback"

n_states = cmd.count_states("qm_min")
if n_states < 1:
    raise RuntimeError("No states found in trajectory.")

if align_states and n_states > 1:
    try:
        cmd.intra_fit("zoom_target", 1)
    except Exception as exc:
        print(f"Warning: intra_fit alignment failed: {{exc}}")

cmd.create("zoom_ref", "zoom_target", 1, 1)
if cmd.count_atoms("zoom_ref") < 1:
    cmd.create("zoom_ref", "context_structure", 1, 1)
    zoom_target_mode = "context-fallback"

mn, mx = cmd.get_extent("zoom_ref")
max_span = max(mx[0]-mn[0], mx[1]-mn[1], mx[2]-mn[2])
# Standard zoom logic without aggressive negative buffers
effective_zoom_buffer = zoom_buffer

cmd.mset(f"1 x{{n_states}}")
cmd.frame(1)
cmd.orient("zoom_ref")
cmd.zoom("zoom_ref", buffer=effective_zoom_buffer)
cmd.delete("zoom_ref")

farthest_water = max([dist for dist, _ in qm_water_keys], default=0.0)
print(
    f"Tagged QM residues using cutoff {{protein_cutoff:.2f}} A around {{qm_residue_name}}: "
    f"core={{len(qm_core_keys)}}, protein={{len(qm_cutoff_keys)}}, waters={{len(qm_water_keys)}} "
    f"(farthest water O distance={{farthest_water:.3f}} A)"
)
print(
    f"Camera zoom target: {{zoom_target_mode}} (atoms={{cmd.count_atoms('zoom_target')}}), "
    f"max-span={{max_span:.2f}} A, requested-buffer={{zoom_buffer:.2f}} A, aligned={{align_states}}"
)
cmd.hide("everything", "qm_min")
cmd.show("cartoon", "context_structure")
cmd.hide("cartoon", "qm_region")
cmd.color("gray85", "context_structure")

# MM lines + rainbow spectrum (matching visualise_dimer.pml)
cmd.select("mm_lines", "non_qm_context and not resn HOH+WAT+SOL")
cmd.hide("cartoon", "mm_lines")
cmd.show("lines", "mm_lines")
cmd.spectrum("count", "rainbow", "mm_lines")

# QM region sticks + spheres
cmd.show("sticks", "qm_region")
cmd.show("spheres", "qm_region")
cmd.set("sphere_scale", 0.25, "qm_region")
cmd.set("stick_radius", 0.16, "qm_region")
util.cbaw("qm_region")
cmd.color("yellow", "qm_region and elem H")
cmd.hide("everything", "solvent_all")
frame_prefix = "{frame_prefix.as_posix()}"
for state in range(1, n_states + 1):
    cmd.frame(state)
    cmd.set("state", state)
    cmd.refresh()
    cmd.png(f"{{frame_prefix}}{{state:04d}}.png", width={width}, height={height}, ray={ray_int})
    print(f"Rendered frame {{state}}/{{n_states}}", flush=True)
python end
quit
"""
    script_path.write_text(script_text, encoding="utf-8")


def collect_frame_pngs(frames_dir: Path, frame_stem: str = "frame_") -> List[Path]:
    frames = sorted(frames_dir.glob(f"{frame_stem}*.png"))
    return [frame.resolve() for frame in frames if frame.is_file()]


def encode_mp4_with_ffmpeg(
    frame_pattern: str,
    fps: int,
    output_path: Path,
    rotate_90: str = "cw",
    crf: int = 16,
    preset: str = "slow",
) -> None:
    vf_chain = []
    if rotate_90 == "cw":
        vf_chain.append("transpose=1")
    elif rotate_90 == "ccw":
        vf_chain.append("transpose=2")
    vf_chain.append("format=yuv420p")

    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-i",
        frame_pattern,
        "-vf",
        ",".join(vf_chain),
        "-c:v",
        "libx264",
        "-preset",
        str(preset),
        "-crf",
        str(crf),
        "-movflags",
        "+faststart",
        str(output_path),
    ]
    subprocess.run(cmd, check=True)


def encode_gif_with_pillow(frame_paths: List[Path], fps: int, output_path: Path) -> None:
    from PIL import Image

    duration_ms = max(int(round(1000 / max(fps, 1))), 1)
    images = [Image.open(path).convert("P", palette=Image.ADAPTIVE) for path in frame_paths]
    first, rest = images[0], images[1:]
    first.save(
        output_path,
        save_all=True,
        append_images=rest,
        duration=duration_ms,
        loop=0,
        optimize=False,
    )


def cleanup_intermediates(paths: List[Path]) -> None:
    for path in paths:
        try:
            if path.is_file():
                path.unlink()
        except Exception:
            pass


def render_frames_with_pymol(pymol_cmd: List[str], pymol_log: Path, append: bool = False) -> int:
    mode = "a" if append else "w"
    with pymol_log.open(mode, encoding="utf-8") as log_handle:
        if append:
            log_handle.write("\n\n=== PyMOL fallback render ===\n")
        result = subprocess.run(
            pymol_cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return int(result.returncode)


def frames_look_blank(frame_paths: List[Path], sample_count: int = 3) -> bool:
    if not frame_paths:
        return False
    from PIL import Image

    sample = frame_paths[: max(1, int(sample_count))]
    for frame_path in sample:
        try:
            with Image.open(frame_path) as img:
                extrema = img.convert("RGB").getextrema()
            if not all(lo == hi for lo, hi in extrema):
                return False
        except Exception:
            return False
    return True


def clear_rendered_frames(frames_dir: Path) -> None:
    for stale in frames_dir.glob("frame_*.png"):
        try:
            stale.unlink()
        except Exception:
            pass


def render_video(
    traj_path: Path,
    output_path: Path,
    frames_dir: Path,
    fps: int,
    width: int,
    height: int,
    align_states: bool,
    zoom_target: str,
    zoom_buffer_a: float,
    ray: bool,
    qm_residue: str,
    qm_protein_cutoff_a: float,
    qm_nearest_waters: int,
    rotate_90: str,
    crf: int,
    preset: str,
    keep_frames: bool,
    allow_gif_fallback: bool,
    pymol_launch: str,
) -> int:
    if not traj_path.exists() or traj_path.stat().st_size == 0:
        print(f"[ERROR] Trajectory missing/empty: {traj_path}")
        return 1

    n_frames = count_pdb_models(traj_path)
    if n_frames < 1:
        print(f"[ERROR] Trajectory has no complete frames: {traj_path}")
        return 1

    pymol_exe = shutil.which("pymol")
    if pymol_exe is None:
        print("[ERROR] PyMOL executable not found in PATH.")
        return 1

    frames_dir = frames_dir.resolve()
    frames_dir.mkdir(parents=True, exist_ok=True)
    frame_prefix = frames_dir / "frame_"
    pymol_script = frames_dir / "_render_minimization_movie.pml"
    pymol_log = frames_dir / "pymol_render.log"

    clear_rendered_frames(frames_dir)

    write_pymol_script(
        script_path=pymol_script,
        traj_path=traj_path,
        frame_prefix=frame_prefix,
        width=width,
        height=height,
        align_states=align_states,
        zoom_target=zoom_target,
        zoom_buffer_a=zoom_buffer_a,
        ray=ray,
        qm_residue=qm_residue,
        qm_protein_cutoff_a=qm_protein_cutoff_a,
        qm_nearest_waters=qm_nearest_waters,
    )

    print(f"[*] Trajectory: {traj_path}")
    print(f"[*] Frames detected: {n_frames}")
    display_env = os.environ.get("DISPLAY", "").strip()
    launch_mode = pymol_launch
    if launch_mode == "auto":
        if (not ray) and display_env:
            launch_mode = "gui"
        else:
            launch_mode = "headless"

    if launch_mode == "gui":
        if not display_env:
            print("[ERROR] --pymol-launch=gui requested, but DISPLAY is not set.")
            return 1
        pymol_cmd = [pymol_exe, "-q", str(pymol_script)]
    else:
        pymol_cmd = [pymol_exe, "-cq", str(pymol_script)]

    mode_desc = "gui/OpenGL (GPU-capable)" if launch_mode == "gui" else "headless/offscreen (typically CPU)"
    if ray:
        print("[*] PyMOL ray tracing is enabled; frame rendering is CPU-based.")
    print(f"[*] Rendering PNG frames with PyMOL ({mode_desc})...")
    render_code = render_frames_with_pymol(pymol_cmd, pymol_log, append=False)
    frame_paths = collect_frame_pngs(frames_dir)

    if render_code != 0 and launch_mode == "gui":
        print(
            f"[!] GUI/OpenGL PyMOL render failed ({render_code}); retrying in headless mode. "
            f"See {pymol_log}"
        )
        clear_rendered_frames(frames_dir)
        fallback_cmd = [pymol_exe, "-cq", str(pymol_script)]
        render_code = render_frames_with_pymol(fallback_cmd, pymol_log, append=True)
        frame_paths = collect_frame_pngs(frames_dir)

    if render_code == 0 and launch_mode == "gui" and (not ray) and frames_look_blank(frame_paths):
        print(
            "[!] PyMOL GPU/OpenGL frames appear blank; retrying with headless ray tracing "
            "(CPU fallback for correctness)."
        )
        clear_rendered_frames(frames_dir)
        write_pymol_script(
            script_path=pymol_script,
            traj_path=traj_path,
            frame_prefix=frame_prefix,
            width=width,
            height=height,
            align_states=align_states,
            zoom_target=zoom_target,
            zoom_buffer_a=zoom_buffer_a,
            ray=True,
            qm_residue=qm_residue,
            qm_protein_cutoff_a=qm_protein_cutoff_a,
            qm_nearest_waters=qm_nearest_waters,
        )
        fallback_cmd = [pymol_exe, "-cq", str(pymol_script)]
        render_code = render_frames_with_pymol(fallback_cmd, pymol_log, append=True)
        frame_paths = collect_frame_pngs(frames_dir)

    if render_code != 0:
        print(f"[ERROR] PyMOL render failed ({render_code}). See {pymol_log}")
        return render_code

    if not frame_paths:
        print(f"[ERROR] No frames were rendered in {frames_dir}. See {pymol_log}")
        return 1

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ffmpeg_exe = shutil.which("ffmpeg")

    suffix = output_path.suffix.lower()
    if suffix in ("", ".mp4"):
        if suffix == "":
            output_path = output_path.with_suffix(".mp4")
        if ffmpeg_exe is None:
            if not allow_gif_fallback:
                print(
                    "[ERROR] ffmpeg not found for MP4 output. Install ffmpeg or re-run with "
                    "--allow-gif-fallback."
                )
                return 1
            gif_path = output_path.with_suffix(".gif")
            print("[!] ffmpeg not found; writing GIF fallback instead of MP4.")
            encode_gif_with_pillow(frame_paths, fps, gif_path)
            print(f"[Done] Wrote animation: {gif_path}")
        else:
            print("[*] Encoding MP4 with ffmpeg...")
            pattern = str((frames_dir / "frame_%04d.png").resolve())
            try:
                encode_mp4_with_ffmpeg(
                    pattern,
                    fps,
                    output_path,
                    rotate_90=rotate_90,
                    crf=crf,
                    preset=preset,
                )
                print(f"[Done] Wrote video: {output_path}")
            except subprocess.CalledProcessError as exc:
                print(f"[ERROR] ffmpeg failed ({exc.returncode}). See terminal output.")
                return exc.returncode
    elif suffix == ".gif":
        print("[*] Writing GIF animation...")
        encode_gif_with_pillow(frame_paths, fps, output_path)
        print(f"[Done] Wrote animation: {output_path}")
    else:
        print(f"[ERROR] Unsupported output extension: {suffix}. Use .mp4 or .gif")
        return 1

    if not keep_frames:
        cleanup_targets = frame_paths + [pymol_script]
        cleanup_intermediates(cleanup_targets)
        print(f"[*] Cleaned up intermediate frames in {frames_dir}")
    else:
        print(f"[*] Kept frames and script in {frames_dir}")

    print(f"[*] PyMOL log: {pymol_log}")
    return 0


def parse_video_args(raw: str) -> Dict[str, object]:
    if not raw:
        return {}
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--fps", type=int)
    parser.add_argument("--width", type=int)
    parser.add_argument("--height", type=int)
    parser.add_argument("--align-states", dest="align_states", action="store_true")
    parser.add_argument("--no-align-states", dest="align_states", action="store_false")
    parser.add_argument("--zoom-target", choices=("protein", "context", "qm"))
    parser.add_argument("--zoom-buffer-a", type=float)
    parser.add_argument("--ray", dest="ray", action="store_true")
    parser.add_argument("--no-ray", dest="ray", action="store_false")
    parser.add_argument("--pymol-launch", choices=("auto", "headless", "gui"))
    parser.add_argument("--rotate-90", choices=("cw", "ccw", "none"))
    parser.add_argument("--crf", type=int)
    parser.add_argument("--preset")
    parser.add_argument("--keep-frames", dest="keep_frames", action="store_true")
    parser.add_argument("--cleanup-frames", dest="keep_frames", action="store_false")
    parser.add_argument("--allow-gif-fallback", dest="allow_gif_fallback", action="store_true")
    parser.add_argument("--qm-residue")
    parser.add_argument("--qm-protein-cutoff-a", type=float)
    parser.add_argument("--qm-nearest-waters", type=int)
    parser.set_defaults(
        align_states=None,
        ray=None,
        keep_frames=None,
        allow_gif_fallback=None,
    )
    overrides = vars(parser.parse_args(shlex.split(raw)))
    return {key: value for key, value in overrides.items() if value is not None}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run dimer OpenMM NVT relaxation and render a video (standalone)."
    )
    parser.add_argument("--pdb", type=Path, default=Path("venus_dimer.pdb"), help="Input dimer PDB")
    parser.add_argument("--workdir", type=Path, default=Path("tc_dimer_nvt"), help="Output workdir")
    parser.add_argument(
        "--overwrite-workdir",
        action="store_true",
        default=True,
        help="Overwrite existing workdir (default: enabled)",
    )
    parser.add_argument(
        "--no-overwrite-workdir",
        dest="overwrite_workdir",
        action="store_false",
        help="Do not overwrite existing workdir",
    )

    parser.add_argument("--ph", type=float, default=7.0, help="pH for protonation")
    parser.add_argument("--padding-a", type=float, default=10.0, help="Solvent padding (Angstrom)")
    parser.add_argument(
        "--skip-solvation",
        action="store_true",
        default=False,
        help="Skip addSolvent and use dimer composition as provided (default: disabled)",
    )
    parser.add_argument(
        "--add-solvent",
        dest="skip_solvation",
        action="store_false",
        help="Enable addSolvent step (periodic box + ions) (default: enabled)",
    )
    parser.add_argument("--ionic-strength-m", type=float, default=0.15, help="Ionic strength (M)")
    parser.add_argument("--temperature-k", type=float, default=300, help="NVT temperature (K)")
    parser.add_argument("--friction-ps", type=float, default=1.0, help="Langevin friction (1/ps)")
    parser.add_argument("--timestep-fs", type=float, default=2.0, help="NVT timestep (fs)")
    parser.add_argument("--nvt-steps", type=int, default=100000, help="NVT MD steps")
    parser.add_argument("--minimize-iters", type=int, default=10000, help="Max minimization iterations")
    parser.add_argument("--platform", default="CUDA", help="OpenMM platform (default: CUDA)")
    parser.add_argument("--cuda-device", default="0", help="CUDA device index")
    parser.add_argument(
        "--strict-gpu",
        action="store_true",
        default=True,
        help="Require requested platform exactly (default: enabled)",
    )
    parser.add_argument(
        "--allow-platform-fallback",
        dest="strict_gpu",
        action="store_false",
        help="Allow OpenMM fallback platform selection",
    )
    parser.add_argument(
        "--openmm-trajectory-file",
        default="dimer_nvt_trajectory.pdb",
        help="Trajectory filename saved in workdir",
    )
    parser.add_argument(
        "--openmm-trajectory-interval",
        type=int,
        default=500,
        help="Trajectory report interval in steps",
    )
    parser.add_argument(
        "--ignore-external-bonds",
        action="store_true",
        default=True,
        help="Ignore external-bond checks in forcefield matching (default: enabled for dimer fragments)",
    )
    parser.add_argument(
        "--strict-external-bonds",
        dest="ignore_external_bonds",
        action="store_false",
        help="Require strict external-bond matching in forcefield templates",
    )
    parser.add_argument(
        "--extra-simple-args",
        default="",
        help="Deprecated (ignored) for standalone mode",
    )

    parser.add_argument(
        "--make-video",
        action="store_true",
        default=True,
        help="Render a video from the generated dimer NVT trajectory (default: enabled)",
    )
    parser.add_argument(
        "--no-video",
        dest="make_video",
        action="store_false",
        help="Skip video rendering step",
    )
    parser.add_argument(
        "--video-output",
        type=Path,
        default=None,
        help="Output video path (default: <workdir>/dimer_nvt.mp4)",
    )
    parser.add_argument(
        "--video-args",
        default="",
        help="Additional raw args for video rendering (standalone parser)",
    )

    parser.add_argument("--frames-dir", type=Path, default=Path("minimization_frames"))
    parser.add_argument("--fps", type=int, default=12, help="Frames per second (default: 12)")
    parser.add_argument("--width", type=int, default=1920, help="Render width in pixels")
    parser.add_argument("--height", type=int, default=1440, help="Render height in pixels")
    parser.add_argument(
        "--align-states",
        action="store_true",
        default=True,
        help="Align all trajectory states to state 1 using zoom-target atoms (default: enabled)",
    )
    parser.add_argument(
        "--no-align-states",
        dest="align_states",
        action="store_false",
        help="Disable alignment of states to state 1 before rendering",
    )
    parser.add_argument(
        "--zoom-target",
        choices=("protein", "context", "qm"),
        default="context",
        help="Selection used for orient/zoom camera fit (default: context)",
    )
    parser.add_argument(
        "--zoom-buffer-a",
        type=float,
        default=-8.0,
        help="PyMOL zoom buffer in Angstrom (more negative = closer zoom)",
    )
    parser.add_argument(
        "--ray",
        action="store_true",
        default=True,
        help="Use ray tracing for each frame (higher quality, slower) (default: enabled)",
    )
    parser.add_argument(
        "--no-ray",
        dest="ray",
        action="store_false",
        help="Disable ray tracing for faster rendering",
    )
    parser.add_argument(
        "--pymol-launch",
        choices=("auto", "headless", "gui"),
        default="auto",
        help="PyMOL launch mode (default: auto)",
    )
    parser.add_argument(
        "--rotate-90",
        choices=("cw", "ccw", "none"),
        default="none",
        help="Rotate output video by 90 degrees (default: none)",
    )
    parser.add_argument(
        "--crf",
        type=int,
        default=16,
        help="x264 quality factor for MP4 (lower is higher quality; default: 16)",
    )
    parser.add_argument(
        "--preset",
        default="slow",
        help="x264 preset for MP4 encoding (default: slow)",
    )
    parser.add_argument(
        "--keep-frames",
        action="store_true",
        default=True,
        help="Keep intermediate PNG frames and temporary PyMOL script (default: enabled)",
    )
    parser.add_argument(
        "--allow-gif-fallback",
        action="store_true",
        help="If ffmpeg is unavailable for MP4 output, write a GIF fallback instead of failing.",
    )
    parser.add_argument("--qm-residue", default="CR2", help="QM core residue name for visualization tagging")
    parser.add_argument(
        "--qm-protein-cutoff-a",
        type=float,
        default=2.65,
        help="Protein cutoff from QM core used to tag QM residues (Angstrom)",
    )
    parser.add_argument(
        "--qm-nearest-waters",
        type=int,
        default=5,
        help="Nearest waters from core to tag as QM waters",
    )

    return parser.parse_args(argv)


def run_openmm_nvt(args) -> tuple[Path | None, Path, Path]:
    pdb_path = resolve_path(args.pdb)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb_path}")

    workdir = resolve_path(args.workdir)
    if workdir.exists():
        if not args.overwrite_workdir:
            raise RuntimeError(f"Workdir already exists: {workdir}. Use --overwrite-workdir to replace it.")
        safe_remove_directory(workdir, Path.cwd())
    workdir.mkdir(parents=True)

    print("[*] Stage 1/2: Protonate and solvate full PDB")
    prepared_pdb = preprocess_pdb_for_solvation(pdb_path, workdir) if not args.skip_solvation else pdb_path
    fixer = pdbfixer.PDBFixer(filename=str(prepared_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(args.ph)

    base_topology = fixer.topology
    base_positions = fixer.positions

    nonstandard_names = find_nonstandard_residues(base_topology)
    print(f"    - Nonstandard residues: {nonstandard_names if nonstandard_names else 'none'}")

    fallback_nonstandard_xml = None
    if nonstandard_names:
        fallback_nonstandard_xml = workdir / "nonstandard_residues_generic.xml"
        write_generic_forcefield_xml(base_topology, base_positions, nonstandard_names, fallback_nonstandard_xml)
        print(
            f"    - Warning: using generic fallback FF for {nonstandard_names}; "
            "this is approximate and may reduce physical accuracy"
        )
        print(f"    - Wrote fallback nonstandard FF: {fallback_nonstandard_xml}")

    ff_inputs = ["amber14-all.xml", "amber14/tip3pfb.xml"]
    if fallback_nonstandard_xml is not None:
        ff_inputs.append(str(fallback_nonstandard_xml))
    forcefield = ForceField(*ff_inputs)

    modeller = Modeller(base_topology, base_positions)
    if args.skip_solvation:
        print("    - Solvation: skipped (--skip-solvation)")
    else:
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=args.padding_a * unit.angstroms,
            neutralize=True,
            ionicStrength=args.ionic_strength_m * unit.molar,
        )

    solvated_pdb = workdir / "solvated_protonated.pdb"
    with open(solvated_pdb, "w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    print(f"    - Saved: {solvated_pdb}")

    print("[*] Stage 2/2: Classical GPU relaxation (minimize + NVT + minimize)")
    nonbonded_cutoff = choose_cutoff_from_box(modeller.topology, default_nm=1.0)
    has_periodic_box = get_periodic_box_lengths_ang(modeller.topology) is not None
    nonbonded_method = PME if has_periodic_box else CutoffNonPeriodic
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=nonbonded_cutoff,
        constraints=HBonds,
        ignoreExternalBonds=args.ignore_external_bonds,
    )

    platform = pick_platform(args.platform, strict=args.strict_gpu)
    properties: Dict[str, str] = {}
    if platform.getName() == "CUDA":
        properties["DeviceIndex"] = str(args.cuda_device)
        properties["Precision"] = "mixed"
    elif platform.getName() == "OpenCL":
        properties["Precision"] = "mixed"

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        args.friction_ps / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    try:
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception as exc:
        if "CUDA_ERROR_UNSUPPORTED_PTX_VERSION" in str(exc):
            raise RuntimeError(
                "CUDA PTX mismatch detected. Update the NVIDIA driver, or pin this conda env to a "
                "matching CUDA runtime."
            ) from exc
        raise
    simulation.context.setPositions(modeller.positions)

    class _NvtEnergyReporter:
        def __init__(self, interval: int):
            self._interval = max(1, int(interval))
            self._frame = 0

        def describeNextReport(self, simulation):
            steps = self._interval - (simulation.currentStep % self._interval)
            return (steps, False, False, False, True, None)

        def report(self, simulation, state):
            self._frame += 1
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            print(
                f"      [OpenMM NVT] frame={self._frame} step={simulation.currentStep} "
                f"energy={energy:.6f} kJ/mol",
                flush=True,
            )

    minimize_report_interval = 500
    can_stream_minimization = (
        minimize_report_interval > 0
        and hasattr(openmm, "LocalEnergyMinimizer")
        and hasattr(openmm, "MinimizationReporter")
    )

    if hasattr(openmm, "MinimizationReporter"):
        class _MinimizationProgressReporter(openmm.MinimizationReporter):
            def __init__(self, label: str, interval: int):
                super().__init__()
                self._label = str(label)
                self._interval = max(1, int(interval))
                self._frame = 0
                self._last_iteration = None
                self._cycle = 1

            def report(self, iteration, x, grad, stats=None):
                self._frame += 1
                if self._last_iteration is not None and iteration <= self._last_iteration:
                    self._cycle += 1
                self._last_iteration = iteration
                if self._frame == 1 or (self._frame % self._interval == 0):
                    energy = float("nan")
                    if stats is not None:
                        try:
                            energy = float(stats["system energy"])
                        except Exception:
                            pass
                    print(
                        f"      [OpenMM min {self._label}] frame={self._frame} step={iteration} "
                        f"cycle={self._cycle} energy={energy:.6f} kJ/mol",
                        flush=True,
                    )
                return False
    else:
        _MinimizationProgressReporter = None

    openmm_traj_path = None
    if args.nvt_steps > 0 and args.openmm_trajectory_interval and int(args.openmm_trajectory_interval) > 0:
        openmm_traj_path = workdir / args.openmm_trajectory_file
        simulation.reporters.append(PDBReporter(str(openmm_traj_path), int(args.openmm_trajectory_interval)))
        simulation.reporters.append(_NvtEnergyReporter(int(args.openmm_trajectory_interval)))

    def potential_kj_per_mol() -> float:
        state = simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    e_before_min1 = potential_kj_per_mol()
    if can_stream_minimization and _MinimizationProgressReporter is not None:
        reporter = _MinimizationProgressReporter("#1", minimize_report_interval)
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            maxIterations=args.minimize_iters,
            reporter=reporter,
        )
    else:
        simulation.minimizeEnergy(maxIterations=args.minimize_iters)
    e_after_min1 = potential_kj_per_mol()

    if args.nvt_steps > 0:
        simulation.step(int(args.nvt_steps))
        e_before_min2 = potential_kj_per_mol()
    else:
        print("    - NVT skipped (--nvt-steps=0)")
        e_before_min2 = e_after_min1

    if can_stream_minimization and _MinimizationProgressReporter is not None:
        reporter = _MinimizationProgressReporter("#2", minimize_report_interval)
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            maxIterations=args.minimize_iters,
            reporter=reporter,
        )
    else:
        simulation.minimizeEnergy(maxIterations=args.minimize_iters)
    e_after_min2 = potential_kj_per_mol()

    relaxed_state = simulation.context.getState(getPositions=True)
    relaxed_positions = relaxed_state.getPositions()

    relaxed_pdb = workdir / "classical_relaxed.pdb"
    with open(relaxed_pdb, "w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, relaxed_positions, handle, keepIds=True)

    print(f"    - Platform: {platform.getName()}")
    print(f"    - Nonbonded method: {'PME' if has_periodic_box else 'CutoffNonPeriodic'}")
    print(f"    - Nonbonded cutoff: {nonbonded_cutoff.value_in_unit(unit.nanometer):.3f} nm")
    if args.ignore_external_bonds:
        print("    - ForceField matching: ignoreExternalBonds=True")
    print(f"    - OpenMM potential energy before minimization: {e_before_min1:.6f} kJ/mol")
    print(f"    - OpenMM potential energy after minimization #1: {e_after_min1:.6f} kJ/mol")
    if args.nvt_steps > 0:
        print(f"    - OpenMM potential energy before minimization #2 (after NVT): {e_before_min2:.6f} kJ/mol")
    print(f"    - OpenMM potential energy after minimization #2: {e_after_min2:.6f} kJ/mol")
    print(f"    - Saved: {relaxed_pdb}")
    if openmm_traj_path is not None:
        print(
            f"    - OpenMM NVT trajectory: {openmm_traj_path} "
            f"(interval={int(args.openmm_trajectory_interval)} steps)"
        )

    return openmm_traj_path, workdir, relaxed_pdb


def main(argv=None) -> int:
    args = parse_args(argv)

    if args.extra_simple_args:
        print("[!] Warning: --extra-simple-args is ignored in standalone mode.")

    if args.video_args:
        overrides = parse_video_args(args.video_args)
        for key, value in overrides.items():
            setattr(args, key, value)

    try:
        traj_path, workdir, _ = run_openmm_nvt(args)
    except Exception as exc:
        print(f"[ERROR] {exc}")
        return 1

    if args.make_video:
        if traj_path is None:
            print("[ERROR] NVT trajectory was not written; cannot render video.")
            return 1
        if not traj_path.exists() or traj_path.stat().st_size == 0:
            print(f"[ERROR] NVT trajectory missing/empty: {traj_path}")
            return 1

        video_output = args.video_output
        if video_output is None:
            video_output = workdir / "dimer_nvt.mp4"
        video_output = resolve_path(video_output)

        frames_dir = resolve_path(args.frames_dir)

        print("[*] Rendering dimer NVT video...")
        code = render_video(
            traj_path=traj_path,
            output_path=video_output,
            frames_dir=frames_dir,
            fps=int(args.fps),
            width=int(args.width),
            height=int(args.height),
            align_states=bool(args.align_states),
            zoom_target=str(args.zoom_target),
            zoom_buffer_a=float(args.zoom_buffer_a),
            ray=bool(args.ray),
            qm_residue=str(args.qm_residue),
            qm_protein_cutoff_a=float(args.qm_protein_cutoff_a),
            qm_nearest_waters=int(args.qm_nearest_waters),
            rotate_90=str(args.rotate_90),
            crf=int(args.crf),
            preset=str(args.preset),
            keep_frames=bool(args.keep_frames),
            allow_gif_fallback=bool(args.allow_gif_fallback),
            pymol_launch=str(args.pymol_launch),
        )
        if code != 0:
            return code
        print(f"    - Video: {video_output}")

    print("[Done] Dimer NVT completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

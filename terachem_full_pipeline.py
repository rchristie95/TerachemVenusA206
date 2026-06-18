#!/usr/bin/env python3
"""
Standalone full pipeline for old-current workflow.

This file embeds and runs stage logic directly (no dependency on local stage script files):
1) terachem_simple_old.py logic
2) terachem_tddft_old_current.py logic
3) terachem_davydov_coupling_old_current.py logic
"""

# ===== Embedded Stage 1 Code (from terachem_simple_old.py) =====

"""
Simple OpenMM + TeraChem QM/MM setup workflow with true electrostatic embedding.

Pipeline:
1) Protonate and solvate the full PDB with OpenMM/PDBFixer.
2) Relax the full solvated structure on GPU with classical minimization.
3) Build a chemically safer QM boundary (whole residues + user-specified key residues).
4) Prepare QM/MM inputs on deprotonated CR2 + waters:
   - MM point-charge electrostatic embedding
   - Link-atom capping on QM/MM covalent cuts
   - No TeraChem QM optimization run
"""

import argparse
import gc
import importlib.util
import os
import random
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from itertools import combinations
from pathlib import Path

import numpy as np

# ===== Method/Basis Profile Toggle =====
# Production profile
TC_METHOD = "wb97xd3"
TC_BASIS = "6-311g**"
# Fast test profile
CHEAP_METHOD = "hf"
CHEAP_BASIS = "3-21g"
# Toggle here: False -> production profile, True -> cheap test profile
USE_CHEAP_METHOD = False

if USE_CHEAP_METHOD:
    ACTIVE_TC_METHOD = CHEAP_METHOD
    ACTIVE_TC_BASIS = CHEAP_BASIS
else:
    ACTIVE_TC_METHOD = TC_METHOD
    ACTIVE_TC_BASIS = TC_BASIS

# ===== Optional Run Seed Toggle =====
# Paste your seed here and set USE_FIXED_RUN_SEED=True to force reuse.
# CLI --seed still takes priority when provided.
USE_FIXED_RUN_SEED = True
FIXED_RUN_SEED = 1342234088

RUN_SEED = None


def _generate_seed():
    seed = int.from_bytes(os.urandom(8), "big") & 0x7FFFFFFF
    return seed if seed != 0 else 1


def resolve_seed_arg(seed_arg):
    if seed_arg is not None:
        return int(seed_arg)
    if USE_FIXED_RUN_SEED:
        return int(FIXED_RUN_SEED)
    return None


def init_run_seed(seed_arg, workdir=None, announce=True, always_write=False):
    global RUN_SEED
    created = False
    source = "existing"
    if RUN_SEED is None:
        if seed_arg is None:
            seed = _generate_seed()
            source = "auto"
        else:
            seed = int(seed_arg)
            source = "user"
        RUN_SEED = seed
        created = True
        random.seed(seed)
        np.random.seed(seed)
    else:
        seed = RUN_SEED

    if created and announce:
        print(f"[*] Random seed ({source}): {seed}")

    if workdir is not None and (created or always_write):
        seed_path = Path(workdir) / "random_seed.txt"
        with open(seed_path, "w") as handle:
            handle.write(f"{seed}\n")
        if created and announce:
            print(f"    - Seed saved: {seed_path}")

    return seed

def discover_openmm_plugin_dirs():
    candidates = []
    seen = set()

    def add(path):
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


def score_openmm_plugin_dir(plugin_dir):
    plugin_dir = Path(plugin_dir)
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


def choose_best_openmm_plugin_dir():
    candidates = discover_openmm_plugin_dirs()
    if not candidates:
        return None
    ranked = sorted(candidates, key=score_openmm_plugin_dir, reverse=True)
    return ranked[0]


def preferred_openmm_plugin_dir():
    env_dir = os.environ.get("OPENMM_PLUGIN_DIR")
    if env_dir:
        p = Path(env_dir).expanduser()
        if p.exists() and p.is_dir():
            return p.resolve()
    return None


def configure_openmm_env():
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
from openmm.app import ForceField, Modeller, PDBFile, PME, CutoffNonPeriodic, HBonds, Simulation

WATER_RESIDUE_NAMES = {"HOH", "WAT", "SOL"}
STANDARD_PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH",
    "LYN", "CYM", "CYX", "ACE", "NME", "NHE",
}
COMMON_ION_RESIDUES = {
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "NI", "LI", "CS", "RB", "I",
}
FORMAL_CHARGES = {
    "ASP": -1,
    "GLU": -1,
    "ARG": 1,
    "LYS": 1,
    "HIP": 1,
    "HIS": 0,
    "HID": 0,
    "HIE": 0,
}
CCD_EMBEDDED_FALLBACK = {
    "CR2": {
        "formal_charge": 0,
        "oxygen_to_hydrogens": {
            "OH": ["HOH"],
            "OXT": ["HXT"],
        },
    }
}
LJ_BY_ELEMENT = {
    "H": (0.250, 0.0157),
    "C": (0.340, 0.2761),
    "N": (0.325, 0.1700),
    "O": (0.296, 0.2100),
    "S": (0.356, 1.0460),
    "P": (0.374, 0.8368),
}
LINK_BOND_DISTANCE_A = {
    "C": 1.09,
    "N": 1.01,
    "O": 0.96,
    "S": 1.34,
    "P": 1.42,
}
ATOMIC_MASS_DALTON = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
}


def stage1_parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Simple OpenMM + TeraChem QM/MM setup workflow with electrostatic embedding")
    parser.add_argument("--pdb", default="1MYW.pdb", help="Input PDB file")
    parser.add_argument("--workdir", default="tc_simple_old", help="Output working directory")
    parser.add_argument(
        "--overwrite-workdir",
        action="store_true",
        help="Allow deleting and recreating an existing workdir",
    )
    parser.add_argument("--qm-residue", default="CR2", help="QM core residue name")
    parser.add_argument("--ph", type=float, default=7.0, help="pH for protonation")
    parser.add_argument("--padding-a", type=float, default=10.0, help="Solvent padding in Angstrom")
    parser.add_argument(
        "--skip-solvation",
        action="store_true",
        help="Skip adding solvent/ions and keep the input composition as-is",
    )
    parser.add_argument("--ionic-strength-m", type=float, default=0.0, help="Ionic strength in molar")
    parser.add_argument("--minimize-iters", type=int, default=1000000, help="Max minimization iterations")
    parser.add_argument(
        "--ignore-external-bonds",
        action="store_true",
        help="Allow forcefield template matching to ignore external bonds (useful for truncated chains/fragments)",
    )
    parser.add_argument("--qm-protein-cutoff-a", type=float, default=2.65, help="Protein cutoff from CR2 in Angstrom")
    parser.add_argument(
        "--qm-nearest-waters",
        type=int,
        default=5,
        help="Number of nearest water residues to include in QM region (default: 5; others treated by PCM)",
    )
    parser.add_argument(
        "--qm-include-resids",
        default="",
        help="Comma-separated residue IDs to force-include in QM region (e.g. 65,66,150)",
    )
    parser.add_argument(
        "--qm-include-resnames",
        default="",
        help="Comma-separated residue names to force-include in QM region (e.g. HIS,ARG,HOH)",
    )
    parser.add_argument("--tc-method", default=ACTIVE_TC_METHOD, help="TeraChem method for downstream TD-DFT")
    parser.add_argument("--tc-basis", default=ACTIVE_TC_BASIS, help="TeraChem basis for downstream TD-DFT")
    parser.add_argument(
        "--tc-pcm",
        choices=("none", "cosmo", "xppcm"),
        default="cosmo",
        help="TeraChem PCM model (User Guide: `pcm` keyword); use 'none' to disable PCM",
    )
    parser.add_argument(
        "--tc-epsilon",
        type=float,
        default=78.39,
        help="TeraChem PCM solvent dielectric (User Guide: `epsilon` keyword)",
    )
    parser.add_argument(
        "--tc-pcm-grid",
        choices=("polyhedron", "lebedev", "iswig", "swig", "sphere"),
        default="iswig",
        help="TeraChem PCM cavity/grid type (User Guide: `pcm_grid` keyword)",
    )
    parser.add_argument(
        "--tc-solvent-radius",
        type=float,
        default=1.40,
        help="TeraChem PCM solvent probe radius in Angstrom (water default: 1.40 A)",
    )
    parser.add_argument("--tc-charge", type=int, default=None, help="Override QM total charge")
    parser.add_argument(
        "--strict-qm-charge",
        action="store_true",
        help="Require explicit --tc-charge (disable automatic QM charge estimation)",
    )
    parser.add_argument("--tc-spinmult", type=int, default=1, help="QM spin multiplicity")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (default: auto-generate each run unless USE_FIXED_RUN_SEED is enabled)",
    )
    parser.add_argument("--platform", default="CUDA", help="OpenMM platform to use (default: CUDA)")
    parser.add_argument(
        "--strict-platform",
        action="store_true",
        help="Require the requested OpenMM platform exactly (disable fallback to other platforms)",
    )
    parser.add_argument("--cuda-device", default="0", help="CUDA device index (used when platform=CUDA)")
    parser.add_argument(
        "--embedding-cutoff-a",
        type=float,
        default=0.0,
        help="MM point-charge cutoff from QM atoms in Angstrom (default 0 = no cutoff/all MM atoms)",
    )
    parser.add_argument(
        "--embedding-min-distance-a",
        type=float,
        default=1.2,
        help="Exclude MM point charges closer than this distance to any QM atom (Angstrom)",
    )
    parser.add_argument(
        "--embedding-repulsion-inner-a",
        type=float,
        default=1.8,
        help="Short-range QM/MM safeguard: MM point charges at or inside this distance are zeroed (Angstrom)",
    )
    parser.add_argument(
        "--embedding-repulsion-outer-a",
        type=float,
        default=2.8,
        help="Short-range QM/MM safeguard: MM point charges between inner/outer distances are smoothly restored to full charge (Angstrom)",
    )
    parser.add_argument(
        "--embedding-preserve-residue-charge",
        dest="embedding_preserve_residue_charge",
        action="store_true",
        default=True,
        help="Preserve per-residue MM net charge after short-range damping (recommended)",
    )
    parser.add_argument(
        "--embedding-no-preserve-residue-charge",
        dest="embedding_preserve_residue_charge",
        action="store_false",
        help="Disable per-residue MM charge preservation after short-range damping",
    )
    parser.add_argument(
        "--embedding-preserve-total-charge",
        dest="embedding_preserve_total_charge",
        action="store_true",
        default=True,
        help="Preserve total MM point-charge net charge after short-range damping (recommended)",
    )
    parser.add_argument(
        "--embedding-no-preserve-total-charge",
        dest="embedding_preserve_total_charge",
        action="store_false",
        help="Disable total MM net-charge preservation after short-range damping",
    )
    parser.add_argument(
        "--embedding-exclusion-hops",
        type=int,
        default=2,
        help="Exclude MM atoms this many bond hops out from the QM/MM boundary",
    )
    parser.add_argument(
        "--embedding-max-point-charges",
        type=int,
        default=0,
        help="Maximum MM point charges for TeraChem (default 0 = no cap/no truncation)",
    )
    parser.add_argument(
        "--embedding-include-nonqm-water",
        action="store_true",
        help="Include non-QM waters in MM point charges even when PCM is enabled (not recommended)",
    )
    parser.add_argument(
        "--distance-boundary",
        choices=("auto", "nonperiodic", "periodic"),
        default="auto",
        help="Distance convention for QM/MM selection (auto=periodic when box is present)",
    )
    parser.add_argument(
        "--deprotonate-atom-name",
        default="OH",
        help="Preferred CR2 oxygen atom name to deprotonate from (default: OH)",
    )
    parser.add_argument(
        "--deprotonation-max-oh-a",
        type=float,
        default=1.25,
        help="Maximum O-H distance in Angstrom for deprotonation matching",
    )
    parser.add_argument(
        "--allow-heuristic-deprotonation",
        action="store_true",
        help="If preferred deprotonation atom fails, allow heuristic fallback over CR2 oxygens",
    )
    parser.add_argument(
        "--strict-deprotonation",
        action="store_true",
        help="Fail when CR2 deprotonation cannot be assigned",
    )
    parser.add_argument(
        "--nonstandard-ff-xml",
        default="",
        help="Comma-separated FF XML files with validated parameters for nonstandard residues",
    )
    parser.add_argument(
        "--strict-nonstandard-ff",
        action="store_true",
        help="Fail when nonstandard residue templates are missing instead of auto-generating a generic fallback FF",
    )
    parser.add_argument(
        "--disable-ccd",
        action="store_true",
        help="Disable automatic RCSB CCD lookup for the QM core residue",
    )
    parser.add_argument(
        "--ccd-cache-dir",
        default=".ccd_cache",
        help="Directory for cached CCD CIF files",
    )
    parser.add_argument(
        "--stop-after-openmm",
        action="store_true",
        help="Stop after OpenMM protonation/solvation/classical relaxation (skip QM/MM setup)",
    )
    return parser.parse_args(argv)


def list_platform_names():
    return [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]


def load_additional_openmm_plugins():
    for plugin_dir in discover_openmm_plugin_dirs():
        try:
            Platform.loadPluginsFromDirectory(str(plugin_dir))
        except Exception:
            continue


def pick_platform(preferred_name="CUDA", strict=False):
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


def to_numpy_angstrom(positions):
    return np.array(positions.value_in_unit(unit.angstroms))


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_]", "_", text)


def element_class(symbol):
    mapping = {"H": "H", "C": "C", "N": "N", "O": "O", "S": "S", "P": "P"}
    return mapping.get(symbol, safe_name(symbol))


def safe_remove_directory(path, workspace_root):
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


def get_periodic_box_lengths_ang(topology):
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return None
    lengths = np.array([np.linalg.norm(v.value_in_unit(unit.angstrom)) for v in vectors], dtype=float)
    if np.any(lengths <= 1.0e-8):
        return None
    return lengths


def min_distances_to_reference(points, reference_points, box_lengths_a=None, use_periodic=False):
    points = np.asarray(points, dtype=float)
    reference_points = np.asarray(reference_points, dtype=float)
    if points.size == 0 or reference_points.size == 0:
        return np.array([], dtype=float)

    points_count = points.shape[0]
    refs_count = reference_points.shape[0]
    if points_count * refs_count <= 2_000_000:
        deltas = points[:, np.newaxis, :] - reference_points[np.newaxis, :, :]
        if use_periodic and box_lengths_a is not None:
            deltas = deltas - box_lengths_a * np.round(deltas / box_lengths_a)
        distances = np.linalg.norm(deltas, axis=2)
        return np.min(distances, axis=1)

    min_distances = np.empty(points_count, dtype=float)
    chunk_size = 4096
    for start in range(0, points_count, chunk_size):
        end = min(start + chunk_size, points_count)
        chunk = points[start:end]
        deltas = chunk[:, np.newaxis, :] - reference_points[np.newaxis, :, :]
        if use_periodic and box_lengths_a is not None:
            deltas = deltas - box_lengths_a * np.round(deltas / box_lengths_a)
        distances = np.linalg.norm(deltas, axis=2)
        min_distances[start:end] = np.min(distances, axis=1)
    return min_distances


def atom_mass_dalton(atom, symbol):
    if atom.element is not None and atom.element.mass is not None:
        return float(atom.element.mass.value_in_unit(unit.dalton))
    return float(ATOMIC_MASS_DALTON.get(symbol, ATOMIC_MASS_DALTON["C"]))


def find_nonstandard_residues(topology):
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


def write_generic_forcefield_xml(topology, positions, residue_names, xml_path):
    residue_names = sorted(set(residue_names))
    residues_by_name = {}
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
                angle_params.append((
                    type_by_atom_name[name_i],
                    type_by_atom_name[center_name],
                    type_by_atom_name[name_k],
                    theta,
                    300.0,
                ))

        residue_blocks.append((residue_name, atoms, bonds_local, external_bond_atoms, type_by_atom_name))

    with open(xml_path, "w") as handle:
        handle.write("<ForceField>\n")
        handle.write("  <AtomTypes>\n")
        for type_name, class_name, symbol, mass in atom_type_data:
            handle.write(f"    <Type name=\"{type_name}\" class=\"{class_name}\" element=\"{symbol}\" mass=\"{mass:.6f}\"/>\n")
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


def parse_csv_tokens(raw_text):
    if not raw_text:
        return set()
    tokens = [token.strip() for token in raw_text.split(",")]
    return {token for token in tokens if token}


def parse_csv_list(raw_text):
    if not raw_text:
        return []
    tokens = [token.strip() for token in raw_text.split(",")]
    return [token for token in tokens if token]


def residue_names_in_forcefield_xml(xml_path):
    text = Path(xml_path).read_text(errors="replace")
    return set(re.findall(r"<Residue\s+name=\"([^\"]+)\"", text))


def extract_mmcif_loop(lines, required_columns):
    i = 0
    total = len(lines)
    required_columns = list(required_columns)
    while i < total:
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        headers = []
        while i < total and lines[i].strip().startswith("_"):
            headers.append(lines[i].strip())
            i += 1
        if not headers:
            continue
        if not all(col in headers for col in required_columns):
            while i < total:
                token = lines[i].strip()
                if token == "loop_" or token.startswith("_"):
                    break
                i += 1
            continue
        rows = []
        while i < total:
            token = lines[i].strip()
            if not token or token == "#":
                i += 1
                if token == "#":
                    break
                continue
            if token == "loop_" or token.startswith("_"):
                break
            if token.startswith(";"):
                i += 1
                while i < total and not lines[i].startswith(";"):
                    i += 1
                if i < total:
                    i += 1
                continue
            rows.append(token.split())
            i += 1
        return headers, rows
    return [], []


def parse_ccd_cif(cif_path, comp_id):
    comp_id = str(comp_id).strip().upper()
    text = Path(cif_path).read_text(errors="replace")
    lines = text.splitlines()

    formal_charge = None
    match = re.search(r"_chem_comp\.pdbx_formal_charge\s+([^\s]+)", text)
    if match:
        raw = match.group(1).strip()
        if raw not in {"?", "."}:
            try:
                formal_charge = int(float(raw))
            except Exception:
                formal_charge = None

    atom_headers, atom_rows = extract_mmcif_loop(
        lines,
        (
            "_chem_comp_atom.comp_id",
            "_chem_comp_atom.atom_id",
            "_chem_comp_atom.type_symbol",
        ),
    )
    atoms = {}
    if atom_headers:
        comp_idx = atom_headers.index("_chem_comp_atom.comp_id")
        atom_idx = atom_headers.index("_chem_comp_atom.atom_id")
        type_idx = atom_headers.index("_chem_comp_atom.type_symbol")
        for row in atom_rows:
            if len(row) <= max(comp_idx, atom_idx, type_idx):
                continue
            if row[comp_idx].strip().upper() != comp_id:
                continue
            atoms[row[atom_idx].strip()] = row[type_idx].strip().upper()

    bond_headers, bond_rows = extract_mmcif_loop(
        lines,
        (
            "_chem_comp_bond.comp_id",
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
        ),
    )
    neighbors = {}
    if bond_headers:
        comp_idx = bond_headers.index("_chem_comp_bond.comp_id")
        atom1_idx = bond_headers.index("_chem_comp_bond.atom_id_1")
        atom2_idx = bond_headers.index("_chem_comp_bond.atom_id_2")
        for row in bond_rows:
            if len(row) <= max(comp_idx, atom1_idx, atom2_idx):
                continue
            if row[comp_idx].strip().upper() != comp_id:
                continue
            atom1 = row[atom1_idx].strip()
            atom2 = row[atom2_idx].strip()
            neighbors.setdefault(atom1, set()).add(atom2)
            neighbors.setdefault(atom2, set()).add(atom1)

    oxygen_to_hydrogens = {}
    for atom_name, atom_type in atoms.items():
        if atom_type != "O":
            continue
        attached_hydrogens = sorted(
            neighbor for neighbor in neighbors.get(atom_name, set())
            if atoms.get(neighbor) == "H"
        )
        if attached_hydrogens:
            oxygen_to_hydrogens[atom_name] = attached_hydrogens

    return {
        "comp_id": comp_id,
        "formal_charge": formal_charge,
        "atoms": atoms,
        "oxygen_to_hydrogens": oxygen_to_hydrogens,
        "cif_path": str(cif_path),
    }


def load_ccd_component(comp_id, cache_dir):
    comp_id = str(comp_id).strip().upper()
    cache_dir = Path(cache_dir).expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cif_path = cache_dir / f"{comp_id}.cif"

    embedded = CCD_EMBEDDED_FALLBACK.get(comp_id)

    if not cif_path.exists():
        url = f"https://files.rcsb.org/ligands/view/{comp_id}.cif"
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                cif_path.write_bytes(response.read())
        except urllib.error.URLError as exc:
            if embedded is not None:
                return {
                    "comp_id": comp_id,
                    "formal_charge": int(embedded.get("formal_charge", 0)),
                    "atoms": {},
                    "oxygen_to_hydrogens": dict(embedded.get("oxygen_to_hydrogens", {})),
                    "cif_path": None,
                    "source": "embedded",
                }, None
            return None, f"download failed: {exc}"
        except Exception as exc:
            if embedded is not None:
                return {
                    "comp_id": comp_id,
                    "formal_charge": int(embedded.get("formal_charge", 0)),
                    "atoms": {},
                    "oxygen_to_hydrogens": dict(embedded.get("oxygen_to_hydrogens", {})),
                    "cif_path": None,
                    "source": "embedded",
                }, None
            return None, f"download failed: {exc}"

    try:
        info = parse_ccd_cif(cif_path, comp_id)
    except Exception as exc:
        if embedded is not None:
            return {
                "comp_id": comp_id,
                "formal_charge": int(embedded.get("formal_charge", 0)),
                "atoms": {},
                "oxygen_to_hydrogens": dict(embedded.get("oxygen_to_hydrogens", {})),
                "cif_path": None,
                "source": "embedded",
            }, None
        return None, f"parse failed: {exc}"

    if not info.get("atoms"):
        if embedded is not None:
            return {
                "comp_id": comp_id,
                "formal_charge": int(embedded.get("formal_charge", 0)),
                "atoms": {},
                "oxygen_to_hydrogens": dict(embedded.get("oxygen_to_hydrogens", {})),
                "cif_path": None,
                "source": "embedded",
            }, None
        return None, "atom table missing in CCD entry"
    info["source"] = "downloaded"
    return info, None


def select_qm_residues(
    topology,
    positions_ang,
    qm_residue_name,
    protein_cutoff_a,
    nearest_waters=5,
    include_resids=None,
    include_resnames=None,
    box_lengths_a=None,
    use_periodic=False,
):
    residues = list(topology.residues())
    include_resids = include_resids or set()
    include_resnames = include_resnames or set()
    cr2_residues = [r for r in residues if r.name == qm_residue_name]
    if not cr2_residues:
        raise RuntimeError(f"No residue named {qm_residue_name} found for QM core")

    cr2_indices = [a.index for residue in cr2_residues for a in residue.atoms()]
    cr2_coords = positions_ang[cr2_indices]

    selected = set(cr2_residues)
    water_candidates = []
    for residue in residues:
        if residue in selected:
            continue

        if residue.name in WATER_RESIDUE_NAMES:
            residue_indices = [a.index for a in residue.atoms()]
            residue_coords = positions_ang[residue_indices]
            oxygen_indices = [a.index for a in residue.atoms() if a.element is not None and a.element.symbol == "O"]
            probe_coords = positions_ang[oxygen_indices] if oxygen_indices else residue_coords
            distances = min_distances_to_reference(
                probe_coords,
                cr2_coords,
                box_lengths_a=box_lengths_a,
                use_periodic=use_periodic,
            )
            if distances.size:
                water_candidates.append((float(np.min(distances)), residue))
            continue

        residue_indices = [a.index for a in residue.atoms()]
        residue_coords = positions_ang[residue_indices]
        distances = min_distances_to_reference(
            residue_coords,
            cr2_coords,
            box_lengths_a=box_lengths_a,
            use_periodic=use_periodic,
        )
        if np.min(distances) <= protein_cutoff_a:
            selected.add(residue)

    for residue in residues:
        if residue.id in include_resids or residue.name in include_resnames:
            selected.add(residue)

    nearest_waters = max(int(nearest_waters), 0)
    available_water_candidates = sorted(
        ((distance, residue) for distance, residue in water_candidates if residue not in selected),
        key=lambda item: item[0],
    )
    selected_water_candidates = available_water_candidates[:nearest_waters]
    for _, residue in selected_water_candidates:
        selected.add(residue)

    farthest_selected_water_a = None
    if selected_water_candidates:
        farthest_selected_water_a = float(max(distance for distance, _ in selected_water_candidates))

    return (
        selected,
        cr2_residues,
        len(selected_water_candidates),
        len(available_water_candidates),
        farthest_selected_water_a,
    )


def build_qm_atom_records(topology, positions_ang, qm_residues):
    qm_residue_ids = {id(residue) for residue in qm_residues}
    records = []
    for atom in topology.atoms():
        if id(atom.residue) not in qm_residue_ids:
            continue
        symbol = atom.element.symbol if atom.element is not None else "C"
        records.append(
            {
                "global_index": atom.index,
                "atom_name": atom.name,
                "residue_name": atom.residue.name,
                "residue_id": atom.residue.id,
                "symbol": symbol,
                "coord": positions_ang[atom.index].copy(),
                "is_link": False,
            }
        )
    return records


def build_link_atom_records(topology, positions_ang, qm_atom_indices, excluded_mm_indices=None):
    qm_set = set(qm_atom_indices)
    excluded_mm_indices = set(excluded_mm_indices or [])
    if not qm_set:
        return [], []

    links = []
    cut_bonds = []
    link_counter = 1
    for bond in topology.bonds():
        atom1 = bond.atom1
        atom2 = bond.atom2
        in1 = atom1.index in qm_set
        in2 = atom2.index in qm_set
        if in1 == in2:
            continue

        qm_atom = atom1 if in1 else atom2
        mm_atom = atom2 if in1 else atom1
        if mm_atom.index in excluded_mm_indices:
            continue
        qm_coord = positions_ang[qm_atom.index]
        mm_coord = positions_ang[mm_atom.index]
        direction = mm_coord - qm_coord
        norm = float(np.linalg.norm(direction))
        if norm < 1e-8:
            continue
        direction_unit = direction / norm
        qm_symbol = qm_atom.element.symbol if qm_atom.element is not None else "C"
        link_distance = LINK_BOND_DISTANCE_A.get(qm_symbol, 1.09)
        link_coord = qm_coord + link_distance * direction_unit

        links.append(
            {
                "global_index": -link_counter,
                "atom_name": f"L{link_counter}",
                "residue_name": "LNK",
                "residue_id": "0",
                "symbol": "H",
                "coord": link_coord,
                "is_link": True,
                "capped_qm_global_index": qm_atom.index,
                "mm_global_index": mm_atom.index,
            }
        )
        cut_bonds.append((qm_atom, mm_atom, link_distance))
        link_counter += 1

    return links, cut_bonds


def deprotonate_cr2(
    records,
    qm_residue_name,
    preferred_oxygen_name="OH",
    max_oh_distance=1.25,
    allow_heuristic=False,
    ccd_info=None,
):
    cr2_atoms = [r for r in records if r["residue_name"] == qm_residue_name]
    fallback_oxygens = [r for r in cr2_atoms if r["symbol"] == "O"]
    hydrogens = [r for r in cr2_atoms if r["symbol"] == "H"]
    oxygen_by_name = {atom["atom_name"]: atom for atom in fallback_oxygens}
    hydrogen_by_name = {atom["atom_name"]: atom for atom in hydrogens}

    if not fallback_oxygens:
        return records, None, "CR2 oxygen atoms not found in QM region records"
    if not hydrogens:
        return records, None, "already deprotonated (no CR2 hydrogen atoms found)"

    def find_best_pair(oxygen_candidates, hydrogen_candidates):
        best_pair = None
        for oxygen in oxygen_candidates:
            for hydrogen in hydrogen_candidates:
                distance = np.linalg.norm(oxygen["coord"] - hydrogen["coord"])
                if distance > max_oh_distance:
                    continue
                if best_pair is None or distance < best_pair[0]:
                    best_pair = (distance, oxygen, hydrogen)
        return best_pair

    preferred_oxygens = [r for r in fallback_oxygens if r["atom_name"] == preferred_oxygen_name]
    if preferred_oxygen_name and not preferred_oxygens:
        if not allow_heuristic:
            return records, None, f"preferred oxygen '{preferred_oxygen_name}' not found in {qm_residue_name}"

    # Use CCD oxygen-hydrogen connectivity when available.
    ccd_oxygen_to_hydrogens = {}
    if ccd_info:
        ccd_oxygen_to_hydrogens = dict(ccd_info.get("oxygen_to_hydrogens", {}))

    best = None
    mode = "preferred"
    if preferred_oxygens:
        preferred_oxygen = preferred_oxygens[0]
        ccd_h_names = ccd_oxygen_to_hydrogens.get(preferred_oxygen["atom_name"], [])
        if ccd_h_names:
            ccd_hydrogens = [hydrogen_by_name[name] for name in ccd_h_names if name in hydrogen_by_name]
            best = find_best_pair([preferred_oxygen], ccd_hydrogens)
            if best is not None:
                mode = "ccd_preferred"
        if best is None:
            best = find_best_pair([preferred_oxygen], hydrogens)
            if best is not None:
                mode = "preferred"

    if best is None and preferred_oxygens and not allow_heuristic:
        return (
            records,
            None,
            f"no hydrogen within {max_oh_distance:.2f} A of {qm_residue_name}:{preferred_oxygen_name}",
        )

    named_oxygens = [r for r in fallback_oxygens if r["atom_name"] in {"OH", "O", "O1", "O2"}]
    if best is None and named_oxygens:
        best = find_best_pair(named_oxygens, hydrogens)
        mode = "named_oxygen"

    if best is None and ccd_oxygen_to_hydrogens:
        ccd_candidates = []
        for oxygen_name, hydrogen_names in ccd_oxygen_to_hydrogens.items():
            oxygen = oxygen_by_name.get(oxygen_name)
            if oxygen is None:
                continue
            hydrogen_candidates = [hydrogen_by_name[name] for name in hydrogen_names if name in hydrogen_by_name]
            pair = find_best_pair([oxygen], hydrogen_candidates or hydrogens)
            if pair is not None:
                ccd_candidates.append(pair)
        if ccd_candidates:
            best = min(ccd_candidates, key=lambda item: item[0])
            mode = "ccd_oxygen"

    if best is None:
        best = find_best_pair(fallback_oxygens, hydrogens)
        mode = "fallback_oxygen"

    if best is None:
        return records, None, f"no bonded O-H pair found within {max_oh_distance:.2f} A"

    distance, oxygen, removed = best
    removed["deprotonated_from_atom"] = oxygen["atom_name"]
    removed["deprotonation_distance_a"] = float(distance)
    removed["deprotonation_mode"] = mode
    filtered = [r for r in records if r is not removed]
    return filtered, removed, None


def estimate_qm_charge(qm_residues, qm_residue_name, core_is_deprotonated, core_base_charge=0):
    total_charge = 0
    for residue in qm_residues:
        if residue.name == qm_residue_name:
            total_charge += int(core_base_charge)
            if core_is_deprotonated:
                total_charge += -1
        else:
            total_charge += FORMAL_CHARGES.get(residue.name, 0)
    return total_charge


def write_xyz(records, xyz_path, title):
    with open(xyz_path, "w") as handle:
        handle.write(f"{len(records)}\n")
        handle.write(f"{title}\n")
        for rec in records:
            x, y, z = rec["coord"]
            handle.write(f"{rec['symbol']} {x:.8f} {y:.8f} {z:.8f}\n")


def read_xyz_coords(xyz_path):
    # Read the last complete XYZ frame and ignore TeraChem point-charge pseudoatoms ("pnt").
    last_atom_lines = None
    with open(xyz_path, "r") as handle:
        while True:
            natoms_line = handle.readline()
            if not natoms_line:
                break
            natoms_line = natoms_line.strip()
            if not natoms_line:
                continue
            try:
                natoms = int(natoms_line.split()[0])
            except ValueError:
                continue

            comment = handle.readline()
            if not comment:
                break

            atom_lines = []
            for _ in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    atom_lines = []
                    break
                atom_lines.append(atom_line.strip())
            if len(atom_lines) != natoms:
                break
            last_atom_lines = atom_lines

    if last_atom_lines is None:
        raise ValueError(f"No complete XYZ frame found in {xyz_path}")

    coords = []
    symbols = []
    for line in last_atom_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbol = parts[0]
        if symbol.lower() == "pnt":
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        symbols.append(symbol)
        coords.append([x, y, z])

    if not coords:
        raise ValueError(f"No atomic coordinates parsed from {xyz_path}")
    return symbols, np.array(coords, dtype=float)


def choose_cutoff_from_box(topology, default_nm=1.0):
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return default_nm * unit.nanometer
    lengths_nm = [np.linalg.norm(v.value_in_unit(unit.nanometer)) for v in vectors]
    max_allowed = 0.49 * min(lengths_nm)
    cutoff_nm = min(default_nm, max_allowed)
    cutoff_nm = max(cutoff_nm, 0.1)
    return cutoff_nm * unit.nanometer


def get_atomic_charges_from_system(system):
    nonbonded_force = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nonbonded_force = force
            break
    if nonbonded_force is not None:
        charges = np.zeros(system.getNumParticles(), dtype=float)
        for atom_index in range(system.getNumParticles()):
            charge, _, _ = nonbonded_force.getParticleParameters(atom_index)
            charges[atom_index] = charge.value_in_unit(unit.elementary_charge)
        return charges, "NonbondedForce"

    amoeba_cls = getattr(openmm, "AmoebaMultipoleForce", None)
    if amoeba_cls is not None:
        for force in system.getForces():
            if isinstance(force, amoeba_cls):
                charges = np.zeros(system.getNumParticles(), dtype=float)
                for atom_index in range(system.getNumParticles()):
                    params = force.getMultipoleParameters(atom_index)
                    charges[atom_index] = params[0].value_in_unit(unit.elementary_charge)
                return charges, "AmoebaMultipoleForce"

    raise RuntimeError("No supported electrostatic force found (expected NonbondedForce or AmoebaMultipoleForce)")


def select_mm_embedding_indices(
    topology,
    positions_ang,
    qm_indices,
    removed_atom_indices,
    cutoff_a,
    min_distance_a,
    exclusion_hops,
    max_point_charges,
    box_lengths_a=None,
    use_periodic=False,
):
    qm_set = set(qm_indices)
    if not qm_set:
        return [], set(), 0, None

    exclusion_hops = max(int(exclusion_hops), 0)
    cutoff_a = float(cutoff_a)
    min_distance_a = max(float(min_distance_a), 0.0)
    max_point_charges = int(max_point_charges)

    bonded_neighbors = {}
    for atom in topology.atoms():
        bonded_neighbors[atom.index] = set()
    for bond in topology.bonds():
        atom_i = bond.atom1.index
        atom_j = bond.atom2.index
        bonded_neighbors[atom_i].add(atom_j)
        bonded_neighbors[atom_j].add(atom_i)

    excluded = set(qm_set)
    excluded.update(removed_atom_indices)
    boundary_excluded = set()

    if exclusion_hops > 0:
        frontier = set()
        for qm_idx in qm_set:
            for neighbor in bonded_neighbors.get(qm_idx, set()):
                if neighbor not in qm_set:
                    frontier.add(neighbor)
        boundary_excluded.update(frontier)
        current = frontier
        for _ in range(1, exclusion_hops):
            next_frontier = set()
            for atom_idx in current:
                for neighbor in bonded_neighbors.get(atom_idx, set()):
                    if neighbor in qm_set or neighbor in boundary_excluded:
                        continue
                    next_frontier.add(neighbor)
            if not next_frontier:
                break
            boundary_excluded.update(next_frontier)
            current = next_frontier

    excluded.update(boundary_excluded)
    qm_coords = positions_ang[sorted(qm_set)]

    total_atoms = len(positions_ang)
    excluded_mask = np.zeros(total_atoms, dtype=bool)
    if excluded:
        excluded_mask[list(excluded)] = True
    candidate_indices = np.where(~excluded_mask)[0]
    if candidate_indices.size == 0:
        return [], boundary_excluded, 0, None

    candidate_coords = positions_ang[candidate_indices]
    candidate_min_distances = min_distances_to_reference(
        candidate_coords,
        qm_coords,
        box_lengths_a=box_lengths_a,
        use_periodic=use_periodic,
    )

    keep_mask = np.ones(candidate_indices.size, dtype=bool)
    if cutoff_a > 0.0:
        keep_mask &= candidate_min_distances <= cutoff_a
    if min_distance_a > 0.0:
        keep_mask &= candidate_min_distances >= min_distance_a

    selected_indices = candidate_indices[keep_mask]
    selected_min_distances = candidate_min_distances[keep_mask]
    selected_before_cap = int(selected_indices.size)
    cap_distance = None
    if max_point_charges > 0 and selected_before_cap > max_point_charges:
        order = np.argsort(selected_min_distances)
        keep_order = order[:max_point_charges]
        selected_indices = selected_indices[keep_order]
        selected_min_distances = selected_min_distances[keep_order]
        cap_distance = float(np.max(selected_min_distances))

    if selected_indices.size > 1:
        order = np.argsort(selected_min_distances)
        selected_indices = selected_indices[order]

    return selected_indices.tolist(), boundary_excluded, selected_before_cap, cap_distance


def apply_mm_short_range_repulsion(
    atom_indices,
    positions_ang,
    charges,
    qm_reference_coords,
    inner_a,
    outer_a,
    box_lengths_a=None,
    use_periodic=False,
    topology=None,
    preserve_residue_charge=True,
    preserve_total_charge=True,
):
    atom_indices = list(atom_indices)
    empty_stats = {
        "zeroed_count": 0,
        "scaled_count": 0,
        "full_count": 0,
        "closest_distance_a": None,
        "closest_retained_distance_a": None,
        "dropped_nearly_zero_count": 0,
        "residue_rebalance_count": 0,
        "residue_unresolved_count": 0,
        "residue_redistributed_abs_e": 0.0,
        "total_redistributed_abs_e": 0.0,
        "unresolved_total_delta_e": 0.0,
        "net_charge_before_e": 0.0,
        "net_charge_after_e": 0.0,
        "max_abs_charge_delta_e": 0.0,
    }
    if not atom_indices:
        return [], np.array([], dtype=float), empty_stats

    inner_a = max(float(inner_a), 0.0)
    outer_a = max(float(outer_a), 0.0)
    if 0.0 < outer_a < inner_a:
        outer_a = inner_a

    index_array = np.asarray(atom_indices, dtype=int)
    mm_coords = positions_ang[index_array]
    distances = min_distances_to_reference(
        mm_coords,
        qm_reference_coords,
        box_lengths_a=box_lengths_a,
        use_periodic=use_periodic,
    )

    factors = np.ones_like(distances)
    if inner_a > 0.0:
        factors[distances <= inner_a] = 0.0
    if outer_a > inner_a:
        ramp_mask = (distances > inner_a) & (distances < outer_a)
        if np.any(ramp_mask):
            t = (distances[ramp_mask] - inner_a) / (outer_a - inner_a)
            factors[ramp_mask] = t * t * (3.0 - 2.0 * t)

    original_charges = np.asarray(charges[index_array], dtype=float)
    effective_charges = original_charges * factors

    residue_rebalance_count = 0
    residue_unresolved_count = 0
    residue_redistributed_abs = 0.0
    if preserve_residue_charge and topology is not None and index_array.size:
        residue_id_by_atom = np.full(len(charges), -1, dtype=int)
        for atom in topology.atoms():
            residue_id_by_atom[atom.index] = atom.residue.index
        selected_residue_ids = residue_id_by_atom[index_array]
        for residue_id in np.unique(selected_residue_ids):
            if residue_id < 0:
                continue
            residue_mask = selected_residue_ids == residue_id
            delta = float(np.sum(original_charges[residue_mask]) - np.sum(effective_charges[residue_mask]))
            if abs(delta) <= 1.0e-12:
                continue
            candidate_mask = residue_mask & (factors >= 1.0 - 1.0e-12)
            if not np.any(candidate_mask):
                candidate_mask = residue_mask & (factors > 1.0e-12)
            if not np.any(candidate_mask):
                residue_unresolved_count += 1
                continue
            weights = np.abs(original_charges[candidate_mask])
            weight_sum = float(np.sum(weights))
            if weight_sum <= 1.0e-12:
                weights = np.ones(int(np.count_nonzero(candidate_mask)), dtype=float)
                weight_sum = float(weights.size)
            effective_charges[candidate_mask] += delta * (weights / weight_sum)
            residue_rebalance_count += 1
            residue_redistributed_abs += abs(delta)

    total_redistributed_abs = 0.0
    unresolved_total_delta = 0.0
    total_delta = float(np.sum(original_charges) - np.sum(effective_charges))
    if preserve_total_charge and abs(total_delta) > 1.0e-12:
        candidate_mask = factors >= 1.0 - 1.0e-12
        if not np.any(candidate_mask):
            candidate_mask = factors > 1.0e-12
        if np.any(candidate_mask):
            weights = np.abs(original_charges[candidate_mask])
            weight_sum = float(np.sum(weights))
            if weight_sum <= 1.0e-12:
                weights = np.ones(int(np.count_nonzero(candidate_mask)), dtype=float)
                weight_sum = float(weights.size)
            effective_charges[candidate_mask] += total_delta * (weights / weight_sum)
            total_redistributed_abs = abs(total_delta)
        else:
            unresolved_total_delta = total_delta

    keep_mask = np.abs(effective_charges) > 1.0e-12
    kept_indices = index_array[keep_mask].tolist()
    kept_charges = effective_charges[keep_mask]

    retained_distances = distances[keep_mask]
    stats = {
        "zeroed_count": int(np.count_nonzero(factors <= 1.0e-12)),
        "scaled_count": int(np.count_nonzero((factors > 1.0e-12) & (factors < 1.0 - 1.0e-12))),
        "full_count": int(np.count_nonzero(factors >= 1.0 - 1.0e-12)),
        "closest_distance_a": float(np.min(distances)) if distances.size else None,
        "closest_retained_distance_a": float(np.min(retained_distances)) if retained_distances.size else None,
        "dropped_nearly_zero_count": int(np.count_nonzero(~keep_mask)),
        "residue_rebalance_count": residue_rebalance_count,
        "residue_unresolved_count": residue_unresolved_count,
        "residue_redistributed_abs_e": residue_redistributed_abs,
        "total_redistributed_abs_e": total_redistributed_abs,
        "unresolved_total_delta_e": unresolved_total_delta,
        "net_charge_before_e": float(np.sum(original_charges)),
        "net_charge_after_e": float(np.sum(effective_charges)),
        "max_abs_charge_delta_e": float(np.max(np.abs(effective_charges - original_charges))),
    }
    return kept_indices, kept_charges, stats


def write_mm_pointcharges(pc_path, positions_ang, charges, atom_indices, charge_values=None):
    if charge_values is not None and len(charge_values) != len(atom_indices):
        raise ValueError("charge_values length must match atom_indices length")
    with open(pc_path, "w") as handle:
        handle.write(f"{len(atom_indices)}\n")
        handle.write("MM Point Charges\n")
        for i, atom_index in enumerate(atom_indices):
            x, y, z = positions_ang[atom_index]
            charge = charge_values[i] if charge_values is not None else charges[atom_index]
            handle.write(f"{charge:.8f} {x:.8f} {y:.8f} {z:.8f}\n")


def write_qm_setup_settings(
    out_path,
    method,
    basis,
    charge,
    spinmult,
    pcm_model="cosmo",
    pcm_epsilon=78.39,
    pcm_grid="iswig",
    pcm_solvent_radius=1.40,
):
    with open(out_path, "w") as handle:
        handle.write("# Stage 1 QM setup settings for downstream TD-DFT\n")
        handle.write(f"method {method}\n")
        handle.write(f"basis {basis}\n")
        handle.write(f"charge {int(charge)}\n")
        handle.write(f"spinmult {int(spinmult)}\n")
        handle.write(f"pcm {pcm_model}\n")
        handle.write(f"epsilon {float(pcm_epsilon):.6f}\n")
        handle.write(f"pcm_grid {pcm_grid}\n")
        handle.write(f"solvent_radius {float(pcm_solvent_radius):.6f}\n")


def stage1_main(argv=None):
    args = stage1_parse_args(argv)
    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb_path}")

    workdir = Path(args.workdir)
    if workdir.exists():
        if not args.overwrite_workdir:
            raise RuntimeError(
                f"Workdir already exists: {workdir}. "
                "Use --overwrite-workdir to replace it."
            )
        safe_remove_directory(workdir, Path.cwd())
    workdir.mkdir(parents=True)

    init_run_seed(resolve_seed_arg(args.seed), workdir=workdir, announce=True, always_write=True)

    print("[*] Stage 1/4: Protonate and solvate full PDB")
    fixer = pdbfixer.PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(args.ph)

    base_topology = fixer.topology
    base_positions = fixer.positions

    nonstandard_names = find_nonstandard_residues(base_topology)
    print(f"    - Nonstandard residues: {nonstandard_names if nonstandard_names else 'none'}")

    ccd_info = None
    qm_core_base_charge = 0
    if not args.disable_ccd:
        ccd_info, ccd_error = load_ccd_component(args.qm_residue, args.ccd_cache_dir)
        if ccd_info is not None:
            ccd_charge = ccd_info.get("formal_charge")
            if ccd_charge is not None:
                qm_core_base_charge = int(ccd_charge)
            ccd_source = ccd_info.get("source", "unknown")
            oh_sites = ccd_info.get("oxygen_to_hydrogens", {})
            if oh_sites:
                preview = ", ".join(f"{oxygen}->{'/'.join(hs)}" for oxygen, hs in sorted(oh_sites.items())[:3])
                print(
                    f"    - CCD {args.qm_residue} ({ccd_source}): formal_charge={qm_core_base_charge}, "
                    f"oxygen-H sites: {preview}"
                )
            else:
                print(f"    - CCD {args.qm_residue} ({ccd_source}): formal_charge={qm_core_base_charge}")
        else:
            print(
                f"    - Warning: CCD lookup for {args.qm_residue} unavailable ({ccd_error}); "
                "using local residue heuristics"
            )

    model_topology = base_topology
    model_positions = base_positions

    provided_nonstandard_xmls = [Path(token).expanduser() for token in parse_csv_list(args.nonstandard_ff_xml)]
    for xml_path in provided_nonstandard_xmls:
        if not xml_path.exists():
            raise FileNotFoundError(f"Nonstandard FF XML not found: {xml_path}")

    provided_residue_templates = set()
    for xml_path in provided_nonstandard_xmls:
        provided_residue_templates.update(residue_names_in_forcefield_xml(xml_path))

    fallback_nonstandard_xml = None
    if nonstandard_names:
        missing_templates = sorted(set(nonstandard_names) - provided_residue_templates)
        if missing_templates:
            if args.strict_nonstandard_ff:
                raise RuntimeError(
                    "Missing nonstandard forcefield templates for "
                    + ", ".join(missing_templates)
                    + ". Supply them via --nonstandard-ff-xml, or rerun without --strict-nonstandard-ff "
                    "to allow approximate fallback parameters."
                )
            fallback_nonstandard_xml = workdir / "nonstandard_residues_generic.xml"
            write_generic_forcefield_xml(model_topology, model_positions, missing_templates, fallback_nonstandard_xml)
            print(
                f"    - Warning: using generic fallback FF for {missing_templates}; "
                "this is approximate and may reduce physical accuracy"
            )
            print(f"    - Wrote fallback nonstandard FF: {fallback_nonstandard_xml}")

    ff_inputs = ["amber14-all.xml", "amber14/tip3pfb.xml"]
    ff_inputs.extend(str(path) for path in provided_nonstandard_xmls)
    if fallback_nonstandard_xml is not None:
        ff_inputs.append(str(fallback_nonstandard_xml))
    if provided_nonstandard_xmls:
        print(f"    - Using user nonstandard FF XMLs: {[str(p) for p in provided_nonstandard_xmls]}")
    forcefield = ForceField(*ff_inputs)

    modeller = Modeller(model_topology, model_positions)
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
    with open(solvated_pdb, "w") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    print(f"    - Saved: {solvated_pdb}")

    print("[*] Stage 2/4: Classical GPU relaxation (minimization)")
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

    platform = pick_platform(args.platform, strict=args.strict_platform)
    properties = {}
    if platform.getName() == "CUDA":
        properties["DeviceIndex"] = str(args.cuda_device)
        properties["Precision"] = "mixed"
    elif platform.getName() == "OpenCL":
        properties["Precision"] = "mixed"

    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    try:
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception as exc:
        if "CUDA_ERROR_UNSUPPORTED_PTX_VERSION" in str(exc):
            raise RuntimeError(
                "CUDA PTX mismatch detected. Your OpenMM runtime is generating newer PTX than the installed "
                "NVIDIA driver supports. Update the NVIDIA driver, or pin this conda env to a matching CUDA "
                "runtime (for example cuda-version/cuda-nvrtc compatible with your driver)."
            ) from exc
        raise
    simulation.context.setPositions(modeller.positions)

    def potential_kj_per_mol(sim_context):
        state = sim_context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    minimize_report_interval = 500
    can_stream_minimization = (
        minimize_report_interval > 0
        and hasattr(openmm, "LocalEnergyMinimizer")
        and hasattr(openmm, "MinimizationReporter")
    )

    if hasattr(openmm, "MinimizationReporter"):
        class _MinimizationProgressReporter(openmm.MinimizationReporter):
            def __init__(self, interval: int):
                super().__init__()
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
                        f"      [OpenMM min] frame={self._frame} step={iteration} "
                        f"cycle={self._cycle} energy={energy:.6f} kJ/mol",
                        flush=True,
                    )
                return False
    else:
        _MinimizationProgressReporter = None

    e_before = potential_kj_per_mol(simulation.context)
    if can_stream_minimization and _MinimizationProgressReporter is not None:
        reporter = _MinimizationProgressReporter(minimize_report_interval)
        openmm.LocalEnergyMinimizer.minimize(
            simulation.context,
            maxIterations=args.minimize_iters,
            reporter=reporter,
        )
    else:
        simulation.minimizeEnergy(maxIterations=args.minimize_iters)
    e_after = potential_kj_per_mol(simulation.context)

    relaxed_state = simulation.context.getState(getPositions=True)
    relaxed_positions = relaxed_state.getPositions()
    relaxed_positions_ang = to_numpy_angstrom(relaxed_positions)

    relaxed_pdb = workdir / "classical_relaxed.pdb"
    with open(relaxed_pdb, "w") as handle:
        PDBFile.writeFile(modeller.topology, relaxed_positions, handle, keepIds=True)
    print(f"    - Platform: {platform.getName()}")
    print(f"    - Nonbonded method: {'PME' if has_periodic_box else 'CutoffNonPeriodic'}")
    print(f"    - Nonbonded cutoff: {nonbonded_cutoff.value_in_unit(unit.nanometer):.3f} nm")
    if args.ignore_external_bonds:
        print("    - ForceField matching: ignoreExternalBonds=True")
    print(f"    - OpenMM potential energy before minimization: {e_before:.6f} kJ/mol")
    print(f"    - OpenMM potential energy after minimization: {e_after:.6f} kJ/mol")
    print(f"    - Saved: {relaxed_pdb}")

    # Free the OpenMM runtime context before launching TeraChem on the same GPU.
    del simulation
    del integrator
    gc.collect()

    if args.stop_after_openmm:
        print("    - Stopping after OpenMM relaxation (--stop-after-openmm)")
        return

    print("[*] Stage 3/4: Generate QM-region forcefield")
    box_lengths_a = get_periodic_box_lengths_ang(modeller.topology)
    if args.distance_boundary == "auto":
        use_periodic_distances = box_lengths_a is not None
    elif args.distance_boundary == "periodic":
        use_periodic_distances = box_lengths_a is not None
        if box_lengths_a is None:
            print("    - Warning: periodic distance mode requested but no periodic box found; using nonperiodic distances")
    else:
        use_periodic_distances = False
    distance_mode_label = "periodic" if use_periodic_distances else "nonperiodic"

    include_resids = parse_csv_tokens(args.qm_include_resids)
    include_resnames = {name.upper() for name in parse_csv_tokens(args.qm_include_resnames)}

    (
        qm_residues,
        _,
        selected_water_count,
        available_water_count,
        farthest_selected_water_a,
    ) = select_qm_residues(
        modeller.topology,
        relaxed_positions_ang,
        args.qm_residue,
        args.qm_protein_cutoff_a,
        nearest_waters=args.qm_nearest_waters,
        include_resids=include_resids,
        include_resnames=include_resnames,
        box_lengths_a=box_lengths_a,
        use_periodic=use_periodic_distances,
    )
    qm_residue_names = sorted({residue.name for residue in qm_residues})
    qm_ff_xml = workdir / "qm_region_forcefield.xml"
    qm_custom_residue_names = sorted(
        {
            name for name in qm_residue_names
            if name not in STANDARD_PROTEIN_RESIDUES
            and name not in WATER_RESIDUE_NAMES
            and name not in COMMON_ION_RESIDUES
        }
    )
    if qm_custom_residue_names:
        write_generic_forcefield_xml(modeller.topology, relaxed_positions, qm_custom_residue_names, qm_ff_xml)
    else:
        qm_ff_xml.write_text("<ForceField/>\n")
    print(f"    - QM residues: {len(qm_residues)} ({', '.join(qm_residue_names)})")
    print(
        f"    - QM nearest waters from {args.qm_residue}: "
        f"{selected_water_count}/{available_water_count} selected (target={args.qm_nearest_waters})"
    )
    if farthest_selected_water_a is not None:
        print(f"    - Farthest selected QM water O distance: {farthest_selected_water_a:.3f} A")
    print(f"    - QM custom residue templates: {qm_custom_residue_names if qm_custom_residue_names else 'none'}")
    print(f"    - Saved: {qm_ff_xml}")

    qm_records = build_qm_atom_records(modeller.topology, relaxed_positions_ang, qm_residues)
    qm_records_deprot, removed_h, deprotonation_note = deprotonate_cr2(
        qm_records,
        args.qm_residue,
        preferred_oxygen_name=args.deprotonate_atom_name,
        max_oh_distance=args.deprotonation_max_oh_a,
        allow_heuristic=args.allow_heuristic_deprotonation,
        ccd_info=ccd_info,
    )
    if removed_h is None:
        note = deprotonation_note or "no bound hydrogen found"
        already_deprotonated = "already deprotonated" in note.lower()
        if args.strict_deprotonation:
            raise RuntimeError(
                f"{args.qm_residue} deprotonation failed: {note}. "
                "Adjust deprotonation options or rerun without --strict-deprotonation."
            )
        if already_deprotonated:
            print(f"    - {args.qm_residue} deprotonation: {note}")
        else:
            print(f"    - Warning: {args.qm_residue} deprotonation not applied ({note})")
    else:
        from_atom = removed_h.get("deprotonated_from_atom", "?")
        bond_len = removed_h.get("deprotonation_distance_a", float("nan"))
        mode = removed_h.get("deprotonation_mode", "unknown")
        print(
            f"    - {args.qm_residue} deprotonation: removed {removed_h['atom_name']} "
            f"(global atom {removed_h['global_index']}) from {from_atom} "
            f"(O-H {bond_len:.3f} A, mode={mode})"
        )

    print("[*] Stage 4/4: QM/MM setup with electrostatic embedding (no TeraChem QM optimization)")
    qm_xyz = workdir / "qm_deprotonated.xyz"
    removed_atom_indices = {removed_h["global_index"]} if removed_h is not None else set()
    qm_global_indices = {record["global_index"] for record in qm_records_deprot}
    link_records, cut_bonds = build_link_atom_records(
        modeller.topology,
        relaxed_positions_ang,
        qm_global_indices,
        excluded_mm_indices=removed_atom_indices,
    )
    qm_records_with_links = qm_records_deprot + link_records
    write_xyz(qm_records_with_links, qm_xyz, "QM region (with link atoms) for downstream TD-DFT")

    qm_summary = workdir / "qm_region_atoms.txt"
    with open(qm_summary, "w") as handle:
        handle.write("local_index\tglobal_index\tresidue\tresid\tatom\tsymbol\tstatus\n")
        for i, rec in enumerate(qm_records_with_links, start=1):
            status = "link_atom" if rec.get("is_link") else "qm_atom"
            handle.write(
                f"{i}\t{rec['global_index']}\t{rec['residue_name']}\t{rec['residue_id']}\t{rec['atom_name']}\t{rec['symbol']}\t{status}\n"
            )

    cuts_report = workdir / "qm_boundary_cuts.txt"
    with open(cuts_report, "w") as handle:
        handle.write("qm_global\tqm_res\tqm_atom\tmm_global\tmm_res\tmm_atom\tlink_distance_a\n")
        for qm_atom, mm_atom, link_distance in cut_bonds:
            handle.write(
                f"{qm_atom.index}\t{qm_atom.residue.name}:{qm_atom.residue.id}\t{qm_atom.name}\t"
                f"{mm_atom.index}\t{mm_atom.residue.name}:{mm_atom.residue.id}\t{mm_atom.name}\t{link_distance:.3f}\n"
            )

    core_is_deprotonated = removed_h is not None
    if not core_is_deprotonated and deprotonation_note:
        if "already deprotonated" in deprotonation_note.lower():
            core_is_deprotonated = True
    atom_charges, charge_source = get_atomic_charges_from_system(system)
    formal_charge_estimate = estimate_qm_charge(
        qm_residues,
        args.qm_residue,
        core_is_deprotonated,
        core_base_charge=qm_core_base_charge,
    )
    qm_partial_charge = float(np.sum(atom_charges[sorted(qm_global_indices)]))
    qm_partial_charge_rounded = int(np.rint(qm_partial_charge))
    if args.tc_charge is None:
        if args.strict_qm_charge:
            raise RuntimeError(
                "QM charge is not explicitly set. Provide --tc-charge for a physically controlled setup. "
                f"Estimates: formal={formal_charge_estimate}, MM-rounded={qm_partial_charge_rounded} "
                f"(raw MM sum={qm_partial_charge:+.4f})."
            )
        tc_charge = formal_charge_estimate
        if formal_charge_estimate != qm_partial_charge_rounded:
            print(
                "    - Warning: QM charge estimates disagree; using formal estimate "
                f"{formal_charge_estimate} (MM-rounded={qm_partial_charge_rounded}, "
                f"raw MM sum={qm_partial_charge:+.4f}). Use --tc-charge to override."
            )
        else:
            print(
                "    - QM charge estimated consistently from topology/MM charges: "
                f"{tc_charge} (raw MM sum={qm_partial_charge:+.4f})"
            )
    else:
        tc_charge = args.tc_charge
        print(
            f"    - QM charge set explicitly: {tc_charge} "
            f"(formal-estimate={formal_charge_estimate}, MM-rounded={qm_partial_charge_rounded})"
        )

    pcm_active = bool(args.tc_pcm and str(args.tc_pcm).lower() != "none")
    non_qm_water_atom_indices = set()
    embedding_removed_atom_indices = set(removed_atom_indices)
    if pcm_active and not args.embedding_include_nonqm_water:
        non_qm_water_atom_indices = {
            atom.index
            for atom in modeller.topology.atoms()
            if atom.residue.name in WATER_RESIDUE_NAMES and atom.index not in qm_global_indices
        }
        embedding_removed_atom_indices.update(non_qm_water_atom_indices)

    mm_pc_indices, boundary_excluded, mm_pc_uncapped_count, cap_distance = select_mm_embedding_indices(
        modeller.topology,
        relaxed_positions_ang,
        qm_global_indices,
        embedding_removed_atom_indices,
        args.embedding_cutoff_a,
        args.embedding_min_distance_a,
        args.embedding_exclusion_hops,
        args.embedding_max_point_charges,
        box_lengths_a=box_lengths_a,
        use_periodic=use_periodic_distances,
    )
    mm_pc_file = workdir / "mm_charges.dat"
    selected_mm_pc_indices = list(mm_pc_indices)
    selected_mm_pc_total_charge = float(np.sum(atom_charges[selected_mm_pc_indices])) if selected_mm_pc_indices else 0.0
    qm_reference_coords = np.array([record["coord"] for record in qm_records_with_links], dtype=float)
    active_mm_pc_indices, active_mm_pc_charges, mm_repulsion_stats = apply_mm_short_range_repulsion(
        selected_mm_pc_indices,
        relaxed_positions_ang,
        atom_charges,
        qm_reference_coords,
        args.embedding_repulsion_inner_a,
        args.embedding_repulsion_outer_a,
        box_lengths_a=box_lengths_a,
        use_periodic=use_periodic_distances,
        topology=modeller.topology,
        preserve_residue_charge=args.embedding_preserve_residue_charge,
        preserve_total_charge=args.embedding_preserve_total_charge,
    )
    write_mm_pointcharges(
        mm_pc_file,
        relaxed_positions_ang,
        atom_charges,
        active_mm_pc_indices,
        charge_values=active_mm_pc_charges,
    )
    mm_pc_total_charge = float(np.sum(active_mm_pc_charges)) if active_mm_pc_charges.size else 0.0

    qm_settings_file = workdir / "qm_setup_settings.in"
    write_qm_setup_settings(
        qm_settings_file,
        args.tc_method,
        args.tc_basis,
        tc_charge,
        args.tc_spinmult,
        args.tc_pcm,
        args.tc_epsilon,
        args.tc_pcm_grid,
        args.tc_solvent_radius,
    )

    print(f"    - QM atoms (real): {len(qm_records_deprot)}")
    print(f"    - Link atoms added: {len(link_records)}")
    print(f"    - Covalent boundary cuts: {len(cut_bonds)}")
    if cut_bonds:
        cut_residues = sorted({f"{mm_atom.residue.name}:{mm_atom.residue.id}" for _, mm_atom, _ in cut_bonds})
        print(f"    - MM-side cut residues: {', '.join(cut_residues)}")
        print("    - Tip: add chemically active residues via --qm-include-resids/--qm-include-resnames")
    print(f"    - Saved: {cuts_report}")
    print(f"    - Embedding source: {charge_source}")
    print(
        f"    - MM point charges: {len(active_mm_pc_indices)} "
        f"(raw {len(selected_mm_pc_indices)}; net {mm_pc_total_charge:+.4f} e, raw net {selected_mm_pc_total_charge:+.4f} e)"
    )
    if non_qm_water_atom_indices:
        print(f"    - Non-QM water atoms excluded from MM embedding (PCM active): {len(non_qm_water_atom_indices)}")
    if mm_pc_uncapped_count > len(selected_mm_pc_indices):
        print(
            f"    - MM point-charge cap applied: {len(selected_mm_pc_indices)}/{mm_pc_uncapped_count} "
            f"(nearest to QM, cutoff ~{cap_distance:.2f} A)"
        )
    if len(active_mm_pc_indices) != len(selected_mm_pc_indices):
        print(
            f"    - MM short-range repulsion filter: kept {len(active_mm_pc_indices)}/{len(selected_mm_pc_indices)} "
            f"(zeroed={mm_repulsion_stats['zeroed_count']}, tapered={mm_repulsion_stats['scaled_count']})"
        )
    print(
        f"    - MM repulsion inner/outer: {args.embedding_repulsion_inner_a:.2f} / "
        f"{args.embedding_repulsion_outer_a:.2f} A"
    )
    print(
        f"    - MM charge preservation: residue={args.embedding_preserve_residue_charge} "
        f"total={args.embedding_preserve_total_charge}"
    )
    if mm_repulsion_stats["closest_distance_a"] is not None:
        print(f"    - Closest MM-QM distance before repulsion filter: {mm_repulsion_stats['closest_distance_a']:.4f} A")
    if mm_repulsion_stats["closest_retained_distance_a"] is not None:
        print(f"    - Closest MM-QM distance after repulsion filter: {mm_repulsion_stats['closest_retained_distance_a']:.4f} A")
    if mm_repulsion_stats["residue_rebalance_count"] or mm_repulsion_stats["total_redistributed_abs_e"]:
        print(
            f"    - MM charge redistribution: residue-groups={mm_repulsion_stats['residue_rebalance_count']} "
            f"(unresolved={mm_repulsion_stats['residue_unresolved_count']}), "
            f"|dQ|_res={mm_repulsion_stats['residue_redistributed_abs_e']:.6f} e, "
            f"|dQ|_total={mm_repulsion_stats['total_redistributed_abs_e']:.6f} e"
        )
    if abs(mm_repulsion_stats["unresolved_total_delta_e"]) > 1.0e-12:
        print(f"    - Warning: unresolved MM total-charge delta {mm_repulsion_stats['unresolved_total_delta_e']:+.6e} e")
    print(
        f"    - MM net charge before/after repulsion: {mm_repulsion_stats['net_charge_before_e']:+.6f} / "
        f"{mm_repulsion_stats['net_charge_after_e']:+.6f} e"
    )
    print(f"    - Max |delta q_i| from repulsion model: {mm_repulsion_stats['max_abs_charge_delta_e']:.6f} e")
    if not active_mm_pc_indices:
        print("    - Warning: no MM point charges selected (check embedding cutoff/exclusions)")
    print(
        f"    - Embedding cutoff/min-distance: {args.embedding_cutoff_a:.2f} / "
        f"{args.embedding_min_distance_a:.2f} A"
    )
    print(f"    - Distance boundary mode: {distance_mode_label} (requested={args.distance_boundary})")
    print(f"    - Boundary-excluded MM atoms: {len(boundary_excluded)}")
    print(f"    - Saved: {mm_pc_file}")
    print(f"    - Saved: {qm_settings_file}")
    print(
        f"    - QM waters kept: {selected_water_count}/{available_water_count} nearest "
        f"to {args.qm_residue} (target={args.qm_nearest_waters})"
    )
    if farthest_selected_water_a is not None:
        print(f"    - Farthest QM-water O distance from {args.qm_residue}: {farthest_selected_water_a:.3f} A")
    if pcm_active:
        print(
            f"    - PCM: model={args.tc_pcm} epsilon={args.tc_epsilon:.3f} "
            f"radius={args.tc_solvent_radius:.3f} A grid={args.tc_pcm_grid}"
        )
    else:
        print("    - PCM: disabled")
    print(f"    - Charge/Spin: {tc_charge}/{args.tc_spinmult}")

    print("    - TeraChem QM optimization step removed; keeping OpenMM-relaxed QM geometry")

    final_topology = modeller.topology
    final_positions_quantity = relaxed_positions_ang * unit.angstroms
    if removed_h is not None:
        atom_to_delete = list(final_topology.atoms())[removed_h["global_index"]]
        deprot_modeller = Modeller(final_topology, final_positions_quantity)
        deprot_modeller.delete([atom_to_delete])
        final_topology = deprot_modeller.topology
        final_positions_quantity = deprot_modeller.positions

    final_pdb = workdir / "final_qmmm_setup_relaxed.pdb"
    with open(final_pdb, "w") as handle:
        PDBFile.writeFile(final_topology, final_positions_quantity, handle, keepIds=True)

    print(f"    - Saved final structure (no QM optimization): {final_pdb}")


# ===== Embedded Stage 2 Code (from terachem_tddft_old_current.py) =====

"""
terachem_tddft_old_current.py

Purpose:
  Run TD-DFT analysis on the most recent geometry from Stage 1 QM/MM setup
  in tc_simple_old, while importing the generated forcefield XML.

Workflow:
  1) Read the newest QM geometry snapshot (prefers qm_deprotonated.xyz).
  2) Import/validate QM forcefield XML via OpenMM ForceField.
  3) Run TD-DFT energy and select brightest state.
  4) Generate transition/difference densities for that state.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from openmm.app import ForceField

# --- Configuration ---
INPUT_DIR = Path("tc_simple_old")
WORKDIR_PREFIX = "tc_tddft_old_current"
NUM_STATES = 10
TC_PATH = os.environ.get("TC_PATH", "terachem")


def parse_qm_setup_settings(input_dir):
    settings = {
        "method": ACTIVE_TC_METHOD,
        "basis": ACTIVE_TC_BASIS,
        "charge": 0,
        "spinmult": 1,
    }
    settings_sources = [input_dir / "qm_setup_settings.in"]
    settings_path = next((p for p in settings_sources if p.exists()), None)
    if settings_path is None:
        return settings

    for raw_line in settings_path.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        key = parts[0].lower()
        value = parts[1]

        if key == "method":
            settings["method"] = value
        elif key == "basis":
            settings["basis"] = value
        elif key == "charge":
            try:
                settings["charge"] = int(float(value))
            except ValueError:
                pass
        elif key == "spinmult":
            try:
                settings["spinmult"] = int(float(value))
            except ValueError:
                pass

    return settings


def read_last_xyz_frame(xyz_path):
    last_frame = None
    frame_index = -1

    with open(xyz_path, "r") as handle:
        while True:
            natoms_line = handle.readline()
            if not natoms_line:
                break
            natoms_line = natoms_line.strip()
            if not natoms_line:
                continue
            try:
                natoms = int(natoms_line.split()[0])
            except ValueError:
                continue

            comment = handle.readline()
            if not comment:
                break

            atom_lines = []
            for _ in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    atom_lines = []
                    break
                atom_lines.append(atom_line.rstrip("\n"))
            if len(atom_lines) != natoms:
                break

            frame_index += 1
            last_frame = (frame_index, comment.strip(), atom_lines)

    return last_frame


def sanitize_symbol(raw_symbol):
    if not raw_symbol:
        return None
    if raw_symbol.lower() == "pnt":
        return None
    if not re.match(r"^[A-Za-z]{1,3}$", raw_symbol):
        return None
    return raw_symbol[0].upper() + raw_symbol[1:].lower()


def records_from_atom_lines(atom_lines):
    records = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbol = sanitize_symbol(parts[0])
        if symbol is None:
            continue
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        records.append((symbol, x, y, z))
    return records


def get_latest_qm_records(input_dir):
    candidates = [
        input_dir / "qm_deprotonated.xyz",
        input_dir / "qm_step.xyz",
    ]
    for candidate in candidates:
        if not candidate.exists() or candidate.stat().st_size == 0:
            continue
        last = read_last_xyz_frame(candidate)
        if last is None:
            continue
        _, _, atom_lines = last
        records = records_from_atom_lines(atom_lines)
        if records:
            return records, "current", candidate

    raise FileNotFoundError(
        f"No usable QM geometry found in {input_dir}. "
        "Checked qm_deprotonated.xyz and qm_step.xyz."
    )


def write_xyz_records(records, xyz_path, title):
    with open(xyz_path, "w") as handle:
        handle.write(f"{len(records)}\n")
        handle.write(f"{title}\n")
        for symbol, x, y, z in records:
            handle.write(f"{symbol:<2} {x: .10f} {y: .10f} {z: .10f}\n")


def import_forcefield(input_dir, workdir):
    qm_xml = input_dir / "qm_region_forcefield.xml"
    nonstandard_xml = input_dir / "nonstandard_residues.xml"

    existing_xmls = [p for p in (nonstandard_xml, qm_xml) if p.exists()]
    if not existing_xmls:
        raise FileNotFoundError(
            f"No local forcefield XML found in {input_dir}. "
            "Expected qm_region_forcefield.xml (and optional nonstandard_residues.xml)."
        )

    validated_xmls = []
    notes = []

    if qm_xml.exists():
        try:
            ForceField("amber14-all.xml", "amber14/tip3pfb.xml", str(qm_xml))
            validated_xmls = [qm_xml]
        except Exception as exc:
            if "same override level" in str(exc):
                ForceField(str(qm_xml))
                validated_xmls = [qm_xml]
                notes.append("Skipped Amber XML merge due to duplicate residue templates in qm_region_forcefield.xml")
            else:
                raise
    elif nonstandard_xml.exists():
        ForceField("amber14-all.xml", "amber14/tip3pfb.xml", str(nonstandard_xml))
        validated_xmls = [nonstandard_xml]

    for xml_path in existing_xmls:
        shutil.copy(xml_path, workdir / xml_path.name)

    return validated_xmls, notes


def write_energy_input(inp_path, settings, use_pointcharges, seed=None):
    with open(inp_path, "w") as handle:
        handle.write("coordinates geometry.xyz\n")
        handle.write("run energy\n")
        handle.write("cis yes\n")
        handle.write(f"basis {settings['basis']}\n")
        handle.write(f"method {settings['method']}\n")
        handle.write(f"charge {settings['charge']}\n")
        handle.write(f"spinmult {settings['spinmult']}\n")
        if seed is not None:
            handle.write(f"seed {int(seed)}\n")
        if use_pointcharges:
            handle.write("pointcharges mm_charges.dat\n")
        handle.write("scrdir scr_energy\n")
        handle.write(f"cisnumstates {NUM_STATES}\n")
        handle.write("cismaxiter 200\n")
        handle.write("cismax 500\n")
        handle.write("scf diis+a\n")
        handle.write("threall 1.0e-13\n")
        handle.write("end\n")


def run_terachem(workdir, inp_name, out_name):
    out_path = workdir / out_name
    with open(out_path, "w") as log:
        result = subprocess.run([TC_PATH, inp_name], cwd=workdir, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"TeraChem failed ({result.returncode}). See {out_path}")
    return out_path


def run_tddft(geom_records, step_label, geom_source, settings):
    workdir = Path(f"{WORKDIR_PREFIX}_{step_label}")
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    seed = init_run_seed(resolve_seed_arg(None), workdir=workdir, announce=True, always_write=True)

    geom_path = workdir / "geometry.xyz"
    write_xyz_records(
        geom_records,
        geom_path,
        f"Latest QM geometry from {geom_source}",
    )
    print(f"[*] Workspace created: {workdir}")
    print(f"    - Geometry source: {geom_source}")
    print(f"    - QM atoms in snapshot: {len(geom_records)}")

    validated_xmls, ff_notes = import_forcefield(INPUT_DIR, workdir)
    print(f"    - Imported forcefield XML: {', '.join(p.name for p in validated_xmls)}")
    for note in ff_notes:
        print(f"    - Forcefield note: {note}")

    src_charges = INPUT_DIR / "mm_charges.dat"
    use_pointcharges = src_charges.exists()
    if use_pointcharges:
        shutil.copy(src_charges, workdir / "mm_charges.dat")
        print("    - Using MM electrostatic embedding from mm_charges.dat")
    else:
        print("    - Warning: mm_charges.dat not found; running gas-phase TD-DFT")

    print(
        f"    - TDDFT settings: {settings['method']}/{settings['basis']} "
        f"charge={settings['charge']} spin={settings['spinmult']}"
    )

    energy_in = workdir / "energy.in"
    write_energy_input(energy_in, settings, use_pointcharges, seed=seed)
    print(f"[*] Running excited-state calculation for {step_label}...")
    out_file = run_terachem(workdir, energy_in.name, "energy.out")

    import terachem_tddft_analysis_big as ta

    states = ta.parse_excitation_energies(out_file)
    print_excited_state_table(states)
    target_root = ta.select_brightest_state(states)

    if target_root is None:
        print("[!] No suitable excited state found.")
        return

    ta.WORKDIR = workdir
    ta.TC_METHOD = settings["method"]
    ta.TC_BASIS = settings["basis"]
    ta.TC_CHARGE = settings["charge"]
    ta.TC_SPIN = settings["spinmult"]
    ta.NUM_STATES = max(NUM_STATES, target_root + 1)
    ta.TC_PATH = TC_PATH

    print(f"[*] Generating transition density for Root {target_root}...")
    ta.generate_densities(Path("geometry.xyz"), Path("mm_charges.dat"), target_root)
    print(f"\n[Done] Snapshot analysis complete in {workdir}")


def stage2_main(input_dir=None, workdir_prefix=None, output_dir=None):
    global INPUT_DIR, WORKDIR_PREFIX
    if input_dir is not None:
        INPUT_DIR = Path(input_dir)
    if workdir_prefix is not None:
        WORKDIR_PREFIX = str(workdir_prefix)
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)

    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Input directory not found: {INPUT_DIR}")

    settings = parse_qm_setup_settings(INPUT_DIR)
    qm_records, step_label, geom_source = get_latest_qm_records(INPUT_DIR)
    prev_cwd = Path.cwd()
    try:
        os.chdir(output_dir)
        run_tddft(qm_records, step_label, geom_source, settings)
    finally:
        os.chdir(prev_cwd)


# ===== Embedded Stage 3 Code (from terachem_davydov_coupling_old_current.py) =====

"""
terachem_davydov_coupling_old_current.py

Calculates the Coulombic coupling (J) and Davydov splitting for a dimer
using the Transition Density Coupling (TDC) method, targeting outputs from
the old-current TDDFT workflow.

Usage:
  python3 terachem_davydov_coupling_old_current.py --workdir tc_tddft_old_current_frame0249
"""

# ===== Pure TDC / coupling building blocks (moved to coupling_core.py) =====
# These functions were extracted verbatim into coupling_core.py so that
# lightweight analysis scripts (sample_coupling_md.py, multipole_decomposition.py,
# lineshape_cd.py) can reuse them without importing OpenMM. Behaviour is unchanged.
import time  # retained: used by re-exported GPU coupling routines' timing
from coupling_core import (
    BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR, HARTREE_TO_EV, HARTREE_TO_CM,
    oscillator_to_dipole_au, print_excited_state_table, parse_pdb_residue_atoms,
    get_super_matrices_with_pymol, apply_pymol_matrix, transition_dipole_au,
    print_transform_matrix, read_dx, calculate_coupling_gpu,
    calculate_coupling_opencl, calculate_coupling, parse_excited_state_candidates,
    autodetect_workdir_and_candidates, find_transition_density_file,
    select_target_state_and_density,
)


def stage3_main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", type=Path, default=Path("tc_tddft_old_current"))
    parser.add_argument(
        "--root",
        type=int,
        help="Excited-state root. If omitted, auto-select using oscillator strength and available density files.",
    )
    parser.add_argument("--monomer", type=Path, default=Path("tc_simple_old/classical_relaxed.pdb"))
    parser.add_argument("--dimer", type=Path, default=Path("venus_dimer.pdb"))
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--thresh", type=float, default=1e-7)
    parser.add_argument("--epsilon", type=float, default=1.77)
    parser.add_argument(
        "--backend",
        choices=["auto", "gpu", "opencl"],
        default="auto",
        help="Integral backend: auto (CUDA->OpenCL), gpu (Numba CUDA), or opencl (PyOpenCL).",
    )
    parser.add_argument(
        "--gpu-chunk",
        type=int,
        default=10000,
        help="Number of source grid points per CUDA chunk.",
    )
    parser.add_argument(
        "--opencl-chunk",
        type=int,
        default=20000,
        help="Number of source grid points per OpenCL chunk.",
    )
    parser.add_argument(
        "--opencl-platform",
        type=int,
        default=None,
        help="OpenCL platform index.",
    )
    parser.add_argument(
        "--opencl-device",
        type=int,
        default=None,
        help="OpenCL device index.",
    )
    parser.add_argument(
        "--density-mode",
        choices=["signed", "auto", "abs"],
        default="auto",
        help="Transition-density source.",
    )
    args = parser.parse_args(argv)

    if not args.monomer.exists():
        monomer_fallbacks = [
            Path("tc_simple_old/final_qmmm_setup_relaxed.pdb"),
            Path("tc_simple_old/classical_relaxed.pdb"),
            Path("tc_simple_new/final_optimized.pdb"),
        ]
        for fallback in monomer_fallbacks:
            if fallback.exists():
                print(f"    - Monomer not found at {args.monomer}; using {fallback}")
                args.monomer = fallback
                break

    if not args.monomer.exists():
        print(f"[!] Monomer PDB not found: {args.monomer}")
        sys.exit(1)
    if not args.dimer.exists():
        print(f"[!] Dimer PDB not found: {args.dimer}")
        sys.exit(1)

    # 0. Get State Info
    requested_workdir = args.workdir
    args.workdir, candidates, tried_dirs = autodetect_workdir_and_candidates(args.workdir)
    if args.workdir != requested_workdir:
        print(f"    - Auto-selected workdir: {args.workdir}")

    if not candidates:
        tried_msg = ", ".join(f"{wd} ({n} states)" for wd, n in tried_dirs)
        print(f"[!] Target state not found. Tried: {tried_msg}")
        sys.exit(1)

    print_excited_state_table(candidates)

    requested_root = args.root
    try:
        target_state, density_file, density_kind, root_select_mode = select_target_state_and_density(
            candidates, args.workdir, args.density_mode, requested_root=requested_root
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"[!] {exc}")
        sys.exit(1)
    args.root = target_state['root']
    if requested_root is None:
        print(
            f"    - Auto-selected root {args.root} "
            f"(mode={root_select_mode}, lambda={target_state['nm']:.1f} nm, f={target_state['osc']:.4f})"
        )

    mu_target = oscillator_to_dipole_au(target_state['ev'], target_state['osc'])
    if not np.isfinite(mu_target):
        print(
            f"[!] Could not compute dipole magnitude for root {args.root} "
            f"(ev={target_state['ev']}, f={target_state['osc']})."
        )
        sys.exit(1)
    print(f"[*] Root {args.root}: f={target_state['osc']:.4f}, Target Dipole={mu_target:.4f} a.u.")

    # 1. Build transforms with PyMOL
    matrix_A, matrix_B, aln_A, aln_B, pymol_err = get_super_matrices_with_pymol(args.monomer, args.dimer)
    if pymol_err:
        print(f"[!] PyMOL error: {pymol_err}")
        sys.exit(1)
    
    print("    - Using PyMOL super transforms.")
    print(f"    - super A: RMS={aln_A[0]:.4f}, atoms={int(aln_A[1])}")
    print(f"    - super B: RMS={aln_B[0]:.4f}, atoms={int(aln_B[1])}")
    print_transform_matrix("Monomer A -> Dimer A", matrix_A)
    print_transform_matrix("Monomer A -> Dimer B", matrix_B)

    # 2. Load full transition density
    print("    - Loading transition density...")
    pts_opt, q_opt = read_dx(density_file, threshold=args.thresh, stride=args.stride)

    # 3. Renormalize
    local_origin = np.mean(pts_opt, axis=0)
    pts_local = pts_opt - local_origin
    dip_vec = np.dot(q_opt, pts_local)
    dip_mag = np.linalg.norm(dip_vec) / BOHR_TO_ANGSTROM
    if dip_mag > 1e-6:
        scale = mu_target / dip_mag
        print(f"    - Renormalizing by {scale:.2f} (Local Grid Dipole={dip_mag:.4f} a.u.)")
        q_opt *= scale

    # 4. Transform and Coupling
    pts_A = apply_pymol_matrix(pts_opt, matrix_A)
    pts_B = apply_pymol_matrix(pts_opt, matrix_B)
    print(f"    - Dimer Separation (Density Centroids): {np.linalg.norm(np.mean(pts_A, axis=0) - np.mean(pts_B, axis=0)):.2f} A")
    
    # CR2 check
    cr2_A = parse_pdb_residue_atoms(args.dimer, res_names=["CR2"], chain_id="A")
    cr2_B = parse_pdb_residue_atoms(args.dimer, res_names=["CR2"], chain_id="B")
    if cr2_A and cr2_B:
        pts_cr2_A = np.array([a['xyz'] for a in cr2_A])
        pts_cr2_B = np.array([a['xyz'] for a in cr2_B])
        diffs = pts_cr2_A[:, None, :] - pts_cr2_B[None, :, :]
        dists = np.sqrt(np.sum(diffs**2, axis=2))
        print(f"    - CR2-CR2 Distances: min={np.min(dists):.2f} A, max={np.max(dists):.2f} A")

    try:
        J = calculate_coupling(
            pts_A,
            q_opt,
            pts_B,
            q_opt,
            backend=args.backend,
            gpu_chunk=args.gpu_chunk,
            opencl_chunk=args.opencl_chunk,
            opencl_platform=args.opencl_platform,
            opencl_device=args.opencl_device,
        ) / args.epsilon
    except Exception as exc:
        print(f"[!] Coupling integral failed: {exc}")
        sys.exit(1)

    # Dipole diagnostics
    origin_A = np.mean(pts_A, axis=0)
    origin_B = np.mean(pts_B, axis=0)
    muA = transition_dipole_au(pts_A, q_opt, origin_angstrom=origin_A)
    muB = transition_dipole_au(pts_B, q_opt, origin_angstrom=origin_B)
    muA_mag = np.linalg.norm(muA)
    muB_mag = np.linalg.norm(muB)
    
    cosang = np.dot(muA, muB) / (muA_mag * muB_mag + 1e-30)
    angle = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
    
    print("\n--- DIPOLE ORIENTATION CHECK ---")
    print(f"muA_local (a.u.): {muA} |muA|={muA_mag:.6f}")
    print(f"muB_local (a.u.): {muB} |muB|={muB_mag:.6f}")
    print(f"muA·muB angle = {angle:.2f} deg")
    
    # Far-field sign
    RA = origin_A * ANGSTROM_TO_BOHR
    RB = origin_B * ANGSTROM_TO_BOHR
    Rvec = RB - RA
    R = np.linalg.norm(Rvec)
    Rhat = Rvec / (R + 1e-30)
    Jdd_num = np.dot(muA, muB) - 3.0 * np.dot(muA, Rhat) * np.dot(muB, Rhat)
    Vdd = Jdd_num / (R**3 * args.epsilon)
    print("\n--- FAR-FIELD DIPOLE-DIPOLE ESTIMATE ---")
    print(f"Vdd: {Vdd*HARTREE_TO_CM:.2f} cm^-1")
        
    print("\n" + "="*40 + "\nRESULTS\n" + "="*40)
    print(f"J: {J:.8f} Hartree ({J*HARTREE_TO_CM:.2f} cm^-1)")
    print(f"Splitting: {2*abs(J*HARTREE_TO_CM):.2f} cm^-1\n" + "="*40)


# ===== Pipeline Orchestrator =====
import shlex
import types


def pipeline_parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run embedded simple_old -> tddft_old_current -> davydov_coupling_old_current"
    )
    parser.add_argument("--cwd", default=".", help="Working directory for the pipeline")
    parser.add_argument("--simple-args", default="", help="Extra args for Stage 1 as one shell-style string")
    parser.add_argument(
        "--tddft-args",
        default="",
        help="Stage 2 overrides as shell-style string (supports: --input-dir, --workdir-prefix)",
    )
    parser.add_argument("--coupling-args", default="", help="Extra args for Stage 3 as one shell-style string")
    parser.add_argument(
        "--overwrite-workdir",
        action="store_true",
        help="If set, appends --overwrite-workdir to Stage 1 unless already provided",
    )
    parser.add_argument(
        "--tddft-workdir-prefix",
        default="tc_tddft_old_current",
        help="Fallback prefix used to infer latest Stage 2 workdir",
    )
    parser.add_argument("--skip-simple", action="store_true", help="Skip Stage 1")
    parser.add_argument("--skip-tddft", action="store_true", help="Skip Stage 2")
    parser.add_argument("--skip-coupling", action="store_true", help="Skip Stage 3")
    parser.add_argument(
        "--visualize-script",
        default="visualise_dimer.pml",
        help="PyMOL script to run at the end of the pipeline",
    )
    parser.add_argument("--skip-visualize", action="store_true", help="Skip final PyMOL visualization stage")
    return parser.parse_args(argv)


def split_arg_string(raw):
    try:
        return shlex.split(raw)
    except ValueError as exc:
        raise RuntimeError(f"Failed to parse arg string '{raw}': {exc}") from exc


def has_flag(tokens, flag):
    return any(tok == flag or tok.startswith(flag + "=") for tok in tokens)


def latest_workdir(cwd, prefix):
    candidates = [p for p in cwd.glob(f"{prefix}*") if p.is_dir()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def run_pymol_visualization(cwd, script_arg):
    script_path = Path(script_arg).expanduser()
    if not script_path.is_absolute():
        script_path = cwd / script_path
    script_path = script_path.resolve()

    if not script_path.exists():
        print(f"[!] Visualization script not found: {script_path}")
        return 1

    pymol_exe = shutil.which("pymol")
    if pymol_exe is None:
        print("[!] PyMOL executable not found in PATH.")
        return 1

    log_path = cwd / "pymol_visualise_dimer.log"
    print(f"[*] Stage 4: Run PyMOL script ({script_path.name})")
    with open(log_path, "w") as log:
        result = subprocess.run(
            [pymol_exe, "-cq", str(script_path)],
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if result.returncode != 0:
        print(f"[ERROR] PyMOL script failed ({result.returncode}). See {log_path}")
        return result.returncode

    print(f"    - PyMOL script completed: {script_path.name}")
    print(f"    - Log: {log_path}")
    return 0


def call_stage(func, *args, **kwargs):
    try:
        func(*args, **kwargs)
        return 0
    except SystemExit as exc:
        code = exc.code
        if isinstance(code, int):
            return code
        return 1 if code else 0


def parse_tddft_overrides(raw):
    tokens = split_arg_string(raw)
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--input-dir", dest="input_dir")
    parser.add_argument("--workdir-prefix", dest="workdir_prefix")
    ns, unknown = parser.parse_known_args(tokens)
    if unknown:
        raise RuntimeError(
            "Unsupported --tddft-args tokens: " + " ".join(unknown) +
            ". Supported: --input-dir, --workdir-prefix"
        )
    return ns


def build_stage2_analysis_module():
    mod = types.ModuleType("terachem_tddft_analysis_big")
    mod.WORKDIR = Path("tc_tddft_old_current")
    mod.TC_METHOD = ACTIVE_TC_METHOD
    mod.TC_BASIS = ACTIVE_TC_BASIS
    mod.TC_CHARGE = 0
    mod.TC_SPIN = 1
    mod.NUM_STATES = 20
    mod.TC_PATH = os.environ.get("TC_PATH", "terachem")

    def parse_excitation_energies(out_file):
        states = {}
        out_file = Path(out_file)
        if not out_file.exists():
            return []

        text = out_file.read_text(errors="replace")
        in_table = False
        for line in text.splitlines():
            if "Final Excited State Results" in line:
                in_table = True
                continue
            if not in_table:
                continue
            if "Printing MM field" in line:
                break
            if "---" in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                try:
                    root = int(parts[0])
                    ev = float(parts[2])
                    osc = float(parts[3])
                    nm = 1239.84193 / ev if ev > 0 else 0.0
                    states[root] = {"root": root, "nm": nm, "osc": osc, "ev": ev}
                except Exception:
                    continue
        return sorted(states.values(), key=lambda item: item["root"])

    def select_brightest_state(states, nm_min=450, nm_max=600):
        candidates = [state for state in states if nm_min <= state["nm"] <= nm_max]
        if not candidates:
            candidates = list(states)
        if not candidates:
            return None
        best = max(candidates, key=lambda state: state["osc"])
        return best["root"]

    def generate_densities(geom_path, charges_path, root):
        workdir = Path(mod.WORKDIR)
        inp_file = workdir / "plot.in"
        scr_dir = workdir / "scr_plot"
        out_file = workdir / "plot.out"

        geom_path = Path(geom_path)
        charges_path = Path(charges_path)
        target_root = max(int(root), 1)
        cis_states = target_root

        with open(inp_file, "w") as handle:
            handle.write(f"coordinates {geom_path.name}\n")
            handle.write("run energy\n")
            handle.write("cis yes\n")
            handle.write(f"basis {mod.TC_BASIS}\n")
            handle.write(f"method {mod.TC_METHOD}\n")
            handle.write(f"charge {mod.TC_CHARGE}\n")
            handle.write(f"spinmult {mod.TC_SPIN}\n")
            if charges_path.exists():
                handle.write(f"pointcharges {charges_path.name}\n")
            handle.write(f"scrdir {scr_dir.name}\n")
            handle.write(f"cisnumstates {cis_states}\n")
            handle.write("cismaxiter 200\n")
            handle.write("cismax 500\n")
            handle.write("scf diis+a\n")
            handle.write("threall 1.0e-13\n")
            handle.write("cisdiffdensity yes\n")
            handle.write("cistransdensity yes\n")
            handle.write(f"cistarget {target_root}\n")
            handle.write("end\n")

        with open(out_file, "w") as log:
            result = subprocess.run(
                [str(mod.TC_PATH), inp_file.name],
                cwd=workdir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if result.returncode != 0:
            raise RuntimeError(f"TeraChem plot run failed ({result.returncode}). See {out_file}")

        diff_src = scr_dir / f"diffdens_{target_root}.dx"
        trans_src = scr_dir / f"transdens_{target_root}.dx"
        if diff_src.exists():
            shutil.copy(diff_src, workdir / f"abs_diffdens_{target_root}.dx")
        if trans_src.exists():
            shutil.copy(trans_src, workdir / f"abs_transdens_{target_root}.dx")

    mod.parse_excitation_energies = parse_excitation_energies
    mod.select_brightest_state = select_brightest_state
    mod.generate_densities = generate_densities
    return mod


def resolve_pipeline_cwd(cwd_arg):
    cwd = Path(cwd_arg).resolve()
    if not cwd.exists():
        raise FileNotFoundError(f"--cwd does not exist: {cwd}")
    if not cwd.is_dir():
        raise NotADirectoryError(f"--cwd is not a directory: {cwd}")
    return cwd


def run_stage1(args):
    if args.skip_simple:
        print("[*] Skipping Stage 1")
        return 0

    stage1_tokens = split_arg_string(args.simple_args)
    if args.overwrite_workdir and not has_flag(stage1_tokens, "--overwrite-workdir"):
        stage1_tokens.append("--overwrite-workdir")

    print("[*] Stage 1: terachem_simple_old")
    rc = call_stage(stage1_main, stage1_tokens)
    if rc != 0:
        print(f"[ERROR] Stage 1 failed with exit code {rc}")
    return rc


def _with_stage2_module(stage2_callable):
    prev_mod = sys.modules.get("terachem_tddft_analysis_big")
    sys.modules["terachem_tddft_analysis_big"] = build_stage2_analysis_module()
    try:
        return stage2_callable()
    finally:
        if prev_mod is None:
            sys.modules.pop("terachem_tddft_analysis_big", None)
        else:
            sys.modules["terachem_tddft_analysis_big"] = prev_mod


def run_stage2(args, cwd):
    if args.skip_tddft:
        print("[*] Skipping Stage 2")
        return 0, latest_workdir(cwd, args.tddft_workdir_prefix)

    overrides = parse_tddft_overrides(args.tddft_args)
    prefix = overrides.workdir_prefix or args.tddft_workdir_prefix

    def _invoke_stage2():
        print("[*] Stage 2: terachem_tddft_old_current")
        return call_stage(
            stage2_main,
            input_dir=overrides.input_dir,
            workdir_prefix=prefix,
            output_dir=cwd,
        )

    rc = _with_stage2_module(_invoke_stage2)
    if rc != 0:
        print(f"[ERROR] Stage 2 failed with exit code {rc}")
        return rc, None
    return 0, latest_workdir(cwd, prefix)


def run_stage3(args, selected_tddft_workdir):
    if args.skip_coupling:
        print("[*] Skipping Stage 3")
        return 0

    coupling_tokens = split_arg_string(args.coupling_args)
    if not has_flag(coupling_tokens, "--workdir"):
        if selected_tddft_workdir is None:
            print(
                "[ERROR] Stage 3 requires a TDDFT workdir, but none could be inferred. "
                'Pass one via --coupling-args "--workdir <dir>".'
            )
            return 1
        coupling_tokens = ["--workdir", str(selected_tddft_workdir), *coupling_tokens]
        print(f"[*] Stage 3 workdir: {selected_tddft_workdir}")

    print("[*] Stage 3: terachem_davydov_coupling_old_current")
    rc = call_stage(stage3_main, coupling_tokens)
    if rc != 0:
        print(f"[ERROR] Stage 3 failed with exit code {rc}")
    return rc


def run_stage4(args, cwd):
    if args.skip_visualize:
        print("[*] Skipping Stage 4")
        return 0

    rc = run_pymol_visualization(cwd, args.visualize_script)
    if rc != 0:
        print(f"[ERROR] Stage 4 failed with exit code {rc}")
    return rc


def pipeline_main(argv=None):
    args = pipeline_parse_args(argv)
    cwd = resolve_pipeline_cwd(args.cwd)
    original_cwd = Path.cwd()

    try:
        os.chdir(cwd)

        rc = run_stage1(args)
        if rc != 0:
            return rc

        rc, selected_tddft_workdir = run_stage2(args, cwd)
        if rc != 0:
            return rc

        rc = run_stage3(args, selected_tddft_workdir)
        if rc != 0:
            return rc

        rc = run_stage4(args, cwd)
        if rc != 0:
            return rc

        print("[Done] Full pipeline completed.")
        return 0
    finally:
        os.chdir(original_cwd)


if __name__ == "__main__":
    try:
        sys.exit(pipeline_main())
    except Exception as exc:
        print(f"[ERROR] {exc}")
        sys.exit(1)

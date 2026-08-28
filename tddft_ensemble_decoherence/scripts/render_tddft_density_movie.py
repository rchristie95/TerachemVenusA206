#!/usr/bin/env python3
"""Render the tracked 100 fs TDDFT transition densities in the paper-video style.

The selected compact density for each site/frame is read from the completed
global state-tracking audit.  Every compact density is already phase-aligned to
the retained reference and excludes the three link atoms.  The full protein
trajectory is stabilized on barrel A with PyMOL ``intra_fit``; a Kabsch map from
the compact CR2 coordinates to the corresponding transformed PyMOL state then
places the density without altering it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DEFAULT_TRAJECTORY = (
    REPO
    / "tc_tandem_nvt_100fs_50frames_20260723"
    / "tandem_nvt_100fs_50frames.dcd"
)
DEFAULT_TOPOLOGY = (
    REPO
    / "tc_tandem_nvt_100fs_50frames_20260723"
    / "solvated_protonated.pdb"
)
DEFAULT_TRACKING = ROOT / "rapid_100fs50_fixedmm/global_state_tracking.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, default=DEFAULT_TRAJECTORY)
    parser.add_argument("--topology", type=Path, default=DEFAULT_TOPOLOGY)
    parser.add_argument("--tracking", type=Path, default=DEFAULT_TRACKING)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "rapid_100fs50_fixedmm/tddft_density_100fs.mp4",
    )
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--fps", type=float, default=10.0)
    parser.add_argument("--rotate-deg", type=float, default=90.0)
    parser.add_argument(
        "--density-representation",
        choices=("surface", "spheres"),
        default="surface",
        help="render genuine volumetric isosurfaces (default) or the legacy point spheres",
    )
    parser.add_argument("--n-density", type=int, default=3600)
    parser.add_argument("--density-radius", type=float, default=0.18)
    parser.add_argument("--density-transparency", type=float, default=0.30)
    parser.add_argument(
        "--surface-spacing",
        type=float,
        default=0.30,
        help="visualisation-grid spacing in Angstrom",
    )
    parser.add_argument(
        "--surface-smear",
        type=float,
        default=0.40,
        help="Gaussian smoothing width in Angstrom",
    )
    parser.add_argument(
        "--surface-pad",
        type=float,
        default=1.0,
        help="padding around retained density points in Angstrom",
    )
    parser.add_argument(
        "--iso-level",
        type=float,
        default=0.06,
        help="positive/negative isosurface level after max-absolute normalisation",
    )
    parser.add_argument("--cartoon-transparency", type=float, default=0.90)
    parser.add_argument(
        "--ribbon-clearance",
        type=float,
        default=7.0,
        help=(
            "hide cartoon residues within this distance (Angstrom) of either "
            "chromophore so ribbon geometry cannot occlude the density; 0 disables"
        ),
    )
    parser.add_argument("--transparency-mode", type=int, default=2)
    parser.add_argument("--max-frames", type=int, default=0)
    parser.add_argument("--keep-frames", action="store_true")
    return parser.parse_args()


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    centered_source = source - source_center
    centered_target = target - target_center
    u, _, vt = np.linalg.svd(centered_source.T @ centered_target)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target_center - source_center @ rotation.T
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_lobes(q: np.ndarray, total: int) -> tuple[np.ndarray, np.ndarray]:
    each = max(1, total // 2)
    positive = np.flatnonzero(q > 0.0)
    negative = np.flatnonzero(q < 0.0)
    positive = positive[np.argsort(q[positive])[-min(each, len(positive)) :]]
    negative = negative[np.argsort(np.abs(q[negative]))[-min(each, len(negative)) :]]
    return positive, negative


def density_grid(
    points: np.ndarray,
    q: np.ndarray,
    *,
    spacing: float,
    smear: float,
    pad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Rasterise signed transition charge onto a smooth regular grid.

    This follows ``voxelize_density.py``, which is used by the repository's
    paper-style STEOM isosurface movies.  The compact transition charges have
    already had link-atom Voronoi cells removed upstream.
    """
    from scipy.ndimage import gaussian_filter

    lower = points.min(axis=0) - pad
    upper = points.max(axis=0) + pad
    dimensions = np.ceil((upper - lower) / spacing).astype(int) + 1
    grid = np.zeros(tuple(dimensions), dtype=np.float64)
    fractional = (points - lower) / spacing
    base = np.floor(fractional).astype(int)
    remainder = fractional - base
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                weight = (
                    np.where(dx, remainder[:, 0], 1.0 - remainder[:, 0])
                    * np.where(dy, remainder[:, 1], 1.0 - remainder[:, 1])
                    * np.where(dz, remainder[:, 2], 1.0 - remainder[:, 2])
                )
                indices = (
                    np.clip(base[:, 0] + dx, 0, dimensions[0] - 1),
                    np.clip(base[:, 1] + dy, 0, dimensions[1] - 1),
                    np.clip(base[:, 2] + dz, 0, dimensions[2] - 1),
                )
                np.add.at(grid, indices, weight * q)
    grid = gaussian_filter(grid, sigma=smear / spacing)
    maximum = float(np.max(np.abs(grid)))
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise RuntimeError("Transition-density grid is empty or non-finite")
    return grid / maximum, lower


def density_meshes(
    points: np.ndarray,
    q: np.ndarray,
    *,
    spacing: float,
    smear: float,
    pad: float,
    level: float,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Extract positive and negative triangle meshes by marching cubes."""
    from skimage import measure

    grid, origin = density_grid(
        points,
        q,
        spacing=spacing,
        smear=smear,
        pad=pad,
    )
    meshes = {}
    for phase, sign in (("positive", 1.0), ("negative", -1.0)):
        vertices, faces, normals, _ = measure.marching_cubes(
            sign * grid,
            level=level,
            spacing=(spacing, spacing, spacing),
        )
        meshes[phase] = (
            origin + vertices,
            np.asarray(faces, dtype=np.int32),
            sign * np.asarray(normals, dtype=float),
        )
    return meshes


def main() -> int:
    args = parse_args()
    try:
        from pymol import cmd, util
        cmd.get_version()
    except Exception:
        import pymol

        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd, util
    from pymol.cgo import BEGIN, COLOR, END, NORMAL, SPHERE, TRIANGLES, VERTEX

    tracking = json.loads(args.tracking.read_text())
    if tracking["status"] != "complete":
        raise RuntimeError("Global state tracking is not complete")
    compact_paths = {
        site: [
            Path(frame["compact_path"])
            for frame in tracking["sites"][site]["frames"]
        ]
        for site in ("A", "B")
    }
    frame_count = len(compact_paths["A"])
    if frame_count <= 0 or len(compact_paths["B"]) != frame_count:
        raise RuntimeError(
            "Tracking file does not contain matching nonempty A/B density paths"
        )
    if args.max_frames:
        frame_count = min(frame_count, args.max_frames)

    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    frames_dir = out.parent / f"{out.stem}_frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    for stale in frames_dir.glob("frm_*.png"):
        stale.unlink()

    cmd.reinitialize()
    cmd.load(str(args.topology), "traj")
    cmd.load_traj(str(args.trajectory), "traj", state=1)
    if cmd.count_states("traj") < frame_count:
        raise RuntimeError("PyMOL loaded fewer trajectory states than requested")

    # Match the existing transition-density movie.
    cmd.bg_color("black")
    for name, value in (
        ("ray_opaque_background", 1),
        ("depth_cue", 0),
        ("ray_trace_mode", 0),
        ("antialias", 2),
        ("ray_shadows", 0),
        ("transparency_mode", args.transparency_mode),
        ("stick_radius", 0.14),
        ("sphere_scale", 0.22),
        ("specular", 0.2),
        ("ambient", 0.32),
        ("hash_max", 240),
        ("auto_zoom", 0),
        ("cartoon_gap_cutoff", 0),
        ("cartoon_fancy_helices", 1),
        ("cartoon_transparency", args.cartoon_transparency),
        ("two_sided_lighting", 1),
    ):
        cmd.set(name, value)
    cmd.hide("everything")
    protein = "traj and polymer.protein"
    cmd.dss("traj")
    cmd.show("cartoon", protein)
    cmd.spectrum("count", "rainbow", protein + " and name CA")
    cmd.spectrum("count", "rainbow", protein)
    cr2 = "traj and resn CR2"
    cmd.show("sticks", cr2)
    cmd.show("spheres", cr2)
    util.cbaw(cr2)
    cmd.color("gray80", cr2)
    cmd.color("yellow", f"{cr2} and elem H")
    clearance_selection = "none"
    clearance_atom_count = 0
    if args.ribbon_clearance > 0.0:
        clearance_selection = "density_ribbon_clearance"
        cmd.select(
            clearance_selection,
            (
                f"byres (({protein}) within {args.ribbon_clearance:.3f} "
                f"of ({cr2}))"
            ),
        )
        clearance_atom_count = int(cmd.count_atoms(clearance_selection, state=1))
        cmd.hide("cartoon", clearance_selection)

    # Same stabilizer used by the preceding NVT/density movies.
    align_selection = "traj and resi 1-229 and name CA"
    cmd.intra_fit(align_selection, 1)
    cmd.orient("traj", state=1)
    cmd.zoom(protein, buffer=3.0, state=1)

    first_cr2 = cmd.get_model(cr2, state=1)
    residue_keys: list[tuple[str, str]] = []
    for atom in first_cr2.atom:
        key = ((atom.chain or "").strip(), str(atom.resi).strip())
        if key not in residue_keys:
            residue_keys.append(key)
    if len(residue_keys) != 2:
        raise RuntimeError(f"Expected two CR2 residues, found {residue_keys}")

    # The compact physical atoms start with the complete 29-atom CR2 in the
    # reference AMBER atom-name order.
    from prepare_surrogate_pilot import reference_mapping

    cr2_names = reference_mapping()[0]
    if len(cr2_names) != 29:
        raise RuntimeError("Reference mapping does not contain 29 CR2 atoms")

    colors = {
        "positive": (0.90, 0.12, 0.12),
        "negative": (0.15, 0.35, 0.95),
    }

    def cloud(points: np.ndarray, rgb: tuple[float, float, float]) -> list[float]:
        result: list[float] = [COLOR, *rgb]
        radius = float(args.density_radius)
        for point in points:
            result.extend(
                [
                    SPHERE,
                    float(point[0]),
                    float(point[1]),
                    float(point[2]),
                    radius,
                ]
            )
        return result

    def surface(
        vertices: np.ndarray,
        normals: np.ndarray,
        faces: np.ndarray,
        rgb: tuple[float, float, float],
    ) -> list[float]:
        """Build a smooth triangle CGO with per-vertex normals."""
        flattened = faces.reshape(-1)
        blocks = np.empty((flattened.size, 8), dtype=float)
        blocks[:, 0] = NORMAL
        blocks[:, 1:4] = normals[flattened]
        blocks[:, 4] = VERTEX
        blocks[:, 5:8] = vertices[flattened]
        return [
            BEGIN,
            TRIANGLES,
            COLOR,
            *rgb,
            *blocks.reshape(-1).tolist(),
            END,
        ]

    placement_rmsd: dict[str, list[float]] = {"A": [], "B": []}
    surface_triangle_counts: dict[str, list[int]] = {"A": [], "B": []}
    degrees_per_frame = float(args.rotate_deg) / max(1, frame_count - 1)
    for frame in range(frame_count):
        state = frame + 1
        cmd.set("state", state)
        model = cmd.get_model(cr2, state=state)
        transformed_by_residue = {key: {} for key in residue_keys}
        for atom in model.atom:
            key = ((atom.chain or "").strip(), str(atom.resi).strip())
            if key in transformed_by_residue:
                transformed_by_residue[key][atom.name.strip()] = np.asarray(
                    atom.coord, dtype=float
                )

        for site, residue_key in zip(("A", "B"), residue_keys):
            compact = np.load(compact_paths[site][frame])
            source_cr2 = np.asarray(compact["physical_atom_coords_ang"][:29], float)
            target_by_name = transformed_by_residue[residue_key]
            if any(name not in target_by_name for name in cr2_names):
                missing = [name for name in cr2_names if name not in target_by_name]
                raise RuntimeError(
                    f"Frame {frame} site {site}: missing PyMOL CR2 atoms {missing}"
                )
            target_cr2 = np.asarray([target_by_name[name] for name in cr2_names])
            rotation, translation, rmsd = kabsch(source_cr2, target_cr2)
            placement_rmsd[site].append(rmsd)
            points = np.asarray(compact["pts_ang"], float)
            q = np.asarray(compact["q"], float)
            if args.density_representation == "surface":
                meshes = density_meshes(
                    points,
                    q,
                    spacing=args.surface_spacing,
                    smear=args.surface_smear,
                    pad=args.surface_pad,
                    level=args.iso_level,
                )
                triangles = 0
                for phase, (vertices, faces, normals) in meshes.items():
                    name = f"td_{site}_{phase}"
                    cmd.delete(name)
                    placed_vertices = vertices @ rotation.T + translation
                    placed_normals = normals @ rotation.T
                    cmd.load_cgo(
                        surface(
                            placed_vertices,
                            placed_normals,
                            faces,
                            colors[phase],
                        ),
                        name,
                    )
                    cmd.set(
                        "cgo_transparency",
                        args.density_transparency,
                        name,
                    )
                    triangles += int(len(faces))
                surface_triangle_counts[site].append(triangles)
            else:
                points = points @ rotation.T + translation
                positive, negative = select_lobes(q, args.n_density)
                for phase, indices in (
                    ("positive", positive),
                    ("negative", negative),
                ):
                    name = f"td_{site}_{phase}"
                    cmd.delete(name)
                    cmd.load_cgo(cloud(points[indices], colors[phase]), name)
                    cmd.set(
                        "cgo_transparency",
                        args.density_transparency,
                        name,
                    )

        if frame > 0 and degrees_per_frame:
            cmd.turn("y", degrees_per_frame)
        cmd.ray(args.width, args.height)
        frame_path = frames_dir / f"frm_{frame:04d}.png"
        cmd.png(str(frame_path), dpi=220)
        if frame % 5 == 0 or frame + 1 == frame_count:
            print(f"[render] frame={frame + 1}/{frame_count}", flush=True)

    pngs = sorted(frames_dir.glob("frm_*.png"))
    if len(pngs) != frame_count:
        raise RuntimeError(f"Rendered {len(pngs)} PNGs, expected {frame_count}")
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("ffmpeg is required to encode the requested MP4")
    command = [
        ffmpeg,
        "-y",
        "-loglevel",
        "error",
        "-framerate",
        str(args.fps),
        "-pattern_type",
        "glob",
        "-i",
        str(frames_dir / "frm_*.png"),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "18",
        "-preset",
        "slow",
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-movflags",
        "+faststart",
        str(out),
    ]
    subprocess.run(command, check=True)
    provenance = {
        "status": "complete",
        "output": str(out),
        "output_sha256": sha256(out),
        "frames": frame_count,
        "fps": args.fps,
        "duration_seconds": frame_count / args.fps,
        "resolution": [args.width, args.height],
        "orbit_degrees": args.rotate_deg,
        "density_representation": args.density_representation,
        "density_points_per_site": (
            args.n_density if args.density_representation == "spheres" else None
        ),
        "density_radius_A": (
            args.density_radius
            if args.density_representation == "spheres"
            else None
        ),
        "density_transparency": args.density_transparency,
        "surface_grid_spacing_A": (
            args.surface_spacing
            if args.density_representation == "surface"
            else None
        ),
        "surface_gaussian_smear_A": (
            args.surface_smear
            if args.density_representation == "surface"
            else None
        ),
        "surface_padding_A": (
            args.surface_pad
            if args.density_representation == "surface"
            else None
        ),
        "surface_normalized_iso_level": (
            args.iso_level
            if args.density_representation == "surface"
            else None
        ),
        "surface_triangle_counts_by_site": {
            site: {
                "mean": float(np.mean(values)),
                "minimum": int(np.min(values)),
                "maximum": int(np.max(values)),
            }
            for site, values in surface_triangle_counts.items()
            if values
        },
        "cartoon_transparency": args.cartoon_transparency,
        "ribbon_clearance_A": args.ribbon_clearance,
        "ribbon_clearance_selection": clearance_selection,
        "ribbon_clearance_atom_count_state_1": clearance_atom_count,
        "transparency_mode": args.transparency_mode,
        "trajectory": str(args.trajectory.resolve()),
        "trajectory_sha256": sha256(args.trajectory),
        "topology": str(args.topology.resolve()),
        "topology_sha256": sha256(args.topology),
        "tracking": str(args.tracking.resolve()),
        "tracking_sha256": sha256(args.tracking),
        "density_phase": (
            "Every compact density is phase-aligned to the retained reference; "
            "red is positive and blue is negative."
        ),
        "link_atoms": "Excluded upstream from every displayed compact density.",
        "placement_rmsd_by_site_A": {
            site: {
                "mean": float(np.mean(values)),
                "maximum": float(np.max(values)),
            }
            for site, values in placement_rmsd.items()
        },
        "camera": {
            "reference_state": 1,
            "protein_alignment": align_selection,
            "alignment_method": "PyMOL intra_fit to state 1",
            "fixed": bool(args.rotate_deg == 0.0),
            "applied_y_orbit_degrees": args.rotate_deg,
        },
        "render_style": (
            "Existing paper transition-density movie: black background, dim "
            "high-transparency rainbow cartoon, gray/white CR2 ball-and-stick, "
            "slightly transparent red/blue volumetric isosurfaces visible "
            "through explicit local clearance windows in the beta-barrel ribbon, "
            "barrel-A stabilization, "
            + (
                "fixed reference-state camera."
                if args.rotate_deg == 0.0
                else "slow y orbit."
            )
        ),
    }
    provenance_path = out.with_suffix(".provenance.json")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"[complete] {out}", flush=True)
    if not args.keep_frames:
        for path in pngs:
            path.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

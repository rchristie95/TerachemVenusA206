#!/usr/bin/env python3
"""
render_nvt_movie.py -- PyMOL movie of an NVT trajectory, styled to match the
published Venus-dimer figures (visualise_dimer.pml).

Run it THROUGH PyMOL (the rendering is done by PyMOL's own ray tracer):
    pymol -cq render_nvt_movie.py -- --traj dimer_nvt_restrained_clean.pdb \
        --out videos/nvt_restrained.mp4 --width 1600 --height 1200 \
        --fps 20 --rotate-deg 150

Visual scheme (mirrors visualise_dimer.pml):
  - white background, depth_cue off, fancy helices, plain ray-trace, antialias;
  - whole protein shown as an all-atom rainbow wireframe (the detailed "MM
    backbone rainbow" of the figures) plus a faint rainbow cartoon for shape;
  - only the two CR2 chromophores drawn as sticks + small spheres with white
    carbons (util.cbaw) and yellow hydrogens.

The transition-density isosurfaces are intentionally omitted: that density was
computed at one QM geometry, so animating it on a classical trajectory would
imply a per-frame recomputation that was not done.

Camera: stabilised on monomer A (intra_fit), then orbits by --rotate-deg over
the whole trajectory.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args():
    argv = sys.argv[1:]
    if "--" in argv:                       # when launched as `pymol -cq script -- ...`
        argv = argv[argv.index("--") + 1:]
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=1200)
    ap.add_argument("--fps", type=float, default=20)
    ap.add_argument("--rotate-deg", type=float, default=150.0)
    ap.add_argument("--rotate-axis", default="y", choices=["x", "y", "z"])
    ap.add_argument("--qm-shell", type=float, default=4.0,
                    help="Deprecated compatibility option; only CR2 is displayed as QM")
    ap.add_argument("--rep", choices=["lines", "cartoon", "density"], default="lines",
                    help="what to highlight on the SAME dimer/camera/trajectory: all-atom "
                         "rainbow 'lines' (default), a rainbow 'cartoon' ribbon, or the "
                         "translucent STEOM transition-density isosurface ('density').")
    ap.add_argument("--cr2-rep", choices=["sticks", "cartoon"], default="sticks",
                    help="draw CR2 as highlighted sticks/spheres (default) or only as cartoon")
    ap.add_argument("--iso-level", type=float, default=0.006,
                    help="transition-density isosurface level (max|f|=1 normalised map)")
    ap.add_argument("--iso-carve", type=float, default=5.0,
                    help="carve the density isosurface within this distance of the chromophore (A)")
    ap.add_argument("--iso-transparency", type=float, default=0.5,
                    help="surface transparency for the density isosurface (0 opaque .. 1 clear)")
    ap.add_argument("--density-dx",
                    default="neo_model/orca_steom/steom_transdens_oldframe.dx",
                    help="axis-aligned old-frame transition-density grid (from voxelize_density.py)")
    ap.add_argument("--monomer-pdb", default="tc_simple_old/classical_relaxed.pdb",
                    help="old-monomer reference the density is registered to (super target ref)")
    ap.add_argument("--stride", type=int, default=1,
                    help="render every Nth state (subsample long trajectories for the movie)")
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--start-frame", type=int, default=1,
                    help="first sampled movie frame to render (1-based; supports parallel chunks)")
    ap.add_argument("--end-frame", type=int, default=0,
                    help="last sampled movie frame to render (inclusive; 0 means final frame)")
    ap.add_argument("--keep-frames", action="store_true")
    return ap.parse_args(argv)


def _load_dx(path):
    """Parse an axis-aligned OpenDX scalar grid -> (grid[nx,ny,nz], origin[3], spacing)."""
    import numpy as np
    nx = ny = nz = None
    origin = None
    h = None
    vals = []
    for line in open(path):
        s = line.split()
        if not s:
            continue
        if line.startswith("object 1"):
            nx, ny, nz = int(s[-3]), int(s[-2]), int(s[-1])
        elif line.startswith("origin"):
            origin = [float(x) for x in s[1:4]]
        elif line.startswith("delta"):
            if h is None:
                h = max(float(x) for x in s[1:4])   # diagonal spacing
        elif s[0][0].isdigit() or s[0][0] in "-.":
            vals.extend(float(x) for x in s)
    grid = np.array(vals).reshape(nx, ny, nz)        # C order (z fastest), matches writer
    return grid, np.array(origin), h


def _parse_barrel_frames(traj):
    """Per MODEL: chain-A CA coords (for the intra_fit stabiliser) and CR2 atom
    dicts {name: xyz} for the two barrels (chain A = FP1, chain C = FP2)."""
    import numpy as np
    frames, cur = [], None
    for l in open(traj):
        if l.startswith("MODEL"):
            cur = {"A_ca": [], "A_cr2": {}, "C_cr2": {}}
        elif l.startswith("ENDMDL"):
            cur["A_ca"] = np.array(cur["A_ca"])
            frames.append(cur); cur = None
        elif l[:6] in ("ATOM  ", "HETATM") and cur is not None:
            ch, resn, an = l[21], l[17:20].strip(), l[12:16].strip()
            xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            if ch == "A" and an == "CA":
                cur["A_ca"].append(xyz)
            if resn == "CR2" and ch in ("A", "C"):
                cur[f"{ch}_cr2"][an] = np.array(xyz)
    return frames


def main():
    args = parse_args()

    # Use PyMOL's renderer. If not already inside PyMOL, launch it headless.
    try:
        from pymol import cmd, util
        cmd.get_version()
    except Exception:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd, util

    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    frames_dir = out.parent / (out.stem + "_frames")
    frames_dir.mkdir(exist_ok=True)
    for stale in frames_dir.glob("frm_*.png"):
        stale.unlink()

    cmd.reinitialize()
    cmd.load(args.traj, "traj")
    n_states = cmd.count_states("traj")
    if args.max_frames and args.max_frames > 0:
        n_states = min(n_states, args.max_frames)
    print(f"[*] {n_states} states from {args.traj}")

    # ---- published-figure render settings (black bg, all-atom rainbow lines,
    #      QM region as sticks/spheres; no cartoon -- matches dimer_view_*.png) ----
    cmd.bg_color("black")
    cmd.set("ray_opaque_background", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    # Shadows are useful on density surfaces, but add substantial ray-tracing
    # cost and no useful depth information to the all-line skeleton.
    cmd.set("ray_shadows", 0 if args.rep == "lines" else 1)
    cmd.set("line_width", 1.4)
    cmd.set("line_smooth", 1)
    cmd.set("stick_radius", 0.14)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.3)
    cmd.set("hash_max", 240)
    # per-frame isosurface (density mode) creates objects every frame; without this
    # PyMOL would re-zoom the camera onto each new surface and lose the framing.
    cmd.set("auto_zoom", 0)
    cmd.set("surface_quality", 1)
    cmd.set("transparency", args.iso_transparency)

    cmd.hide("everything")
    protein = "traj and not resn HOH+WAT+SOL+NA+CL"

    if args.rep in ("cartoon", "density"):
        # Rainbow cartoon ribbon (warm FP1 -> cool FP2), no dashes across the CR2 gap.
        # In density mode the ribbon is faded right back so it is only spatial context
        # for the transition-density isosurface (the highlighted aspect).
        cmd.set("cartoon_gap_cutoff", 0)
        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_transparency", 0.85 if args.rep == "density" else 0.0)
        cmd.dss("traj")
        cmd.show("cartoon", protein)
        cmd.spectrum("count", "rainbow", protein + " and name CA")
        cmd.spectrum("count", "rainbow", protein)
    else:
        # All-atom rainbow wireframe; single continuous spectrum across the whole
        # dimer (gives the warm monomer-A / cool monomer-B split of the figures).
        cmd.show("lines", protein)
        cmd.spectrum("count", "rainbow", protein + " and name CA")
        cmd.spectrum("count", "rainbow", protein)

    # Optionally highlight only CR2.  Cartoon mode intentionally adds no
    # ball-and-stick representation, leaving the STEOM density as the sole
    # highlighted object in a density movie.
    cr2 = "traj and resn CR2"
    if args.cr2_rep == "sticks":
        cmd.show("sticks", cr2)
        cmd.show("spheres", cr2)
        cmd.set("sphere_scale", 0.24 if args.rep == "density" else 0.22, cr2)
        util.cbaw(cr2)
        cmd.color("yellow", f"{cr2} and elem H")
    else:
        cmd.hide("sticks", cr2)
        cmd.hide("spheres", cr2)
        cmd.show("cartoon", cr2)

    try:
        cmd.intra_fit("traj and chain A and name CA", 1)
    except Exception as exc:
        print(f"[!] intra_fit skipped: {exc}")
    cmd.orient("traj")
    cmd.zoom("traj", buffer=5)

    all_states = list(range(1, n_states + 1, max(1, args.stride)))
    start_index = max(0, args.start_frame - 1)
    end_index = len(all_states) if args.end_frame <= 0 else min(len(all_states), args.end_frame)
    states = all_states[start_index:end_index]
    if not states:
        raise ValueError("Selected movie-frame range is empty")

    # ---- density mode: preload the STEOM transition density as a multi-state CGO ----
    # The lobes are triangulated ONCE in the old frame (marching cubes), then each
    # rendered state gets its own rigidly-placed copy loaded into a single multi-state
    # CGO. All placement maths is numpy Kabsch (no PyMOL super/contour/create in the
    # loop -- those segfault after ~30 calls in this build), and the render loop below
    # is byte-for-byte the same stable set-state/turn/ray loop the lines & cartoon
    # videos use. The density is the exact one entering J(t): spec-normalised, in the
    # old frame, stabilised on chain A and seated on each barrel's CR2 chromophore.
    if args.rep == "density":
        import numpy as np
        from skimage import measure
        from pymol.cgo import BEGIN, TRIANGLES, COLOR, NORMAL, VERTEX, END
        from align_steom_density import kabsch, cr2_atoms

        grid, origin, h = _load_dx(args.density_dx)
        cmd.set("two_sided_lighting", 1)

        def _mesh(sign):
            v, f, n, _ = measure.marching_cubes(sign * grid, args.iso_level)
            return (origin + v * h), f.reshape(-1), (sign * n)
        pos_v, pos_fv, pos_n = _mesh(+1.0)
        neg_v, neg_fv, neg_n = _mesh(-1.0)

        def _cgo(V, N, fv, rgb):
            blk = np.empty((fv.size, 8))
            blk[:, 0] = NORMAL; blk[:, 1:4] = N[fv]
            blk[:, 4] = VERTEX; blk[:, 5:8] = V[fv]
            return [BEGIN, TRIANGLES, COLOR, *rgb] + blk.reshape(-1).tolist() + [END]

        mono_cr2 = cr2_atoms(args.monomer_pdb)
        frames = _parse_barrel_frames(args.traj)
        ref_ca = frames[0]["A_ca"]

        def _placement(fr, cr2_key):
            Tr, Tt, _ = kabsch(fr["A_ca"], ref_ca)                 # intra_fit on chain A
            common = sorted(set(mono_cr2) & set(fr[cr2_key]))
            P = np.array([mono_cr2[n] for n in common])            # old-frame monomer CR2
            Q = np.array([fr[cr2_key][n] for n in common])         # this barrel's CR2
            Pr, Pt, _ = kabsch(P, Q)                               # old -> raw barrel
            return Tr @ Pr, Tr @ Pt + Tt                           # compose -> stabilised barrel

        for s in states:
            fr = frames[s - 1]
            for cr2_key, ch in (("A_cr2", "A"), ("C_cr2", "C")):
                R, t = _placement(fr, cr2_key)
                cmd.load_cgo(_cgo(pos_v @ R.T + t, pos_n @ R.T, pos_fv, (0.90, 0.12, 0.12)),
                             f"dens_{ch}_pos", state=s)
                cmd.load_cgo(_cgo(neg_v @ R.T + t, neg_n @ R.T, neg_fv, (0.15, 0.35, 0.95)),
                             f"dens_{ch}_neg", state=s)
        for o in ("dens_A_pos", "dens_A_neg", "dens_C_pos", "dens_C_neg"):
            cmd.set("cgo_transparency", args.iso_transparency, o)

    # Base the turn increment on the complete movie, not the selected chunk, so
    # independently rendered chunks concatenate into one continuous orbit.
    deg = float(args.rotate_deg) / max(1, len(all_states) - 1)
    if start_index > 0 and deg:
        cmd.turn(args.rotate_axis, deg * start_index)
    pad = len(str(n_states))
    for i, s in enumerate(states):
        cmd.set("state", s)
        if i > 0 and deg:
            cmd.turn(args.rotate_axis, deg)
        cmd.ray(args.width, args.height)
        cmd.png(str(frames_dir / f"frm_{s:0{pad}d}.png"), dpi=300)
        if i % 25 == 0 or i == 0:
            print(f"    frame {i+1}/{len(states)} (state {s})")

    pngs = sorted(frames_dir.glob("frm_*.png"))
    if not pngs:
        print("[ERROR] no frames rendered")
        return 1

    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        cmd_ff = [ffmpeg, "-y", "-framerate", str(args.fps),
                  "-pattern_type", "glob", "-i", str(frames_dir / "frm_*.png"),
                  "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
                  "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(out)]
        r = subprocess.run(cmd_ff, capture_output=True, text=True)
        if r.returncode != 0:
            print("[ERROR] ffmpeg:\n" + r.stderr[-1200:])
            return 1
        print(f"[*] Wrote {out} ({len(pngs)} frames @ {args.fps} fps, {args.width}x{args.height})")
        if not args.keep_frames:
            for p in pngs:
                p.unlink()
    else:
        print(f"[!] ffmpeg not found; frames in {frames_dir}")
    return 0


# Run unconditionally: under `pymol -cq script.py` the module name is not
# "__main__", so a normal guard would never fire.
_rc = main()
try:
    from pymol import cmd as _cmd
    _cmd.quit(_rc or 0)          # clean exit when launched via the pymol binary
except Exception:
    sys.exit(_rc or 0)

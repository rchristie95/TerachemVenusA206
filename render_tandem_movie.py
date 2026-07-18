#!/usr/bin/env python3
r"""
render_tandem_movie.py -- PyMOL movie of the unrestrained VenusA206 tandem-dimer
NVT dynamics, showing the A206 interface forming and breaking.

Run through PyMOL (or via the TeraChem-env python, which finish_launches it):
    python render_tandem_movie.py --traj tandem_nvt_whole.pdb \
        --out videos/tandem_unrestrained.mp4

Scheme: white background, cartoon barrels (FP1 orange / FP2 marine), the flexible
linker as a green tube, and the two CR2 chromophores as yellow sticks+spheres so
the coupling geometry is visible. Camera is stabilised on FP1 (chain A) so FP2 is
seen swinging out (undocking) and back (re-docking); a slow orbit adds depth.
Expects the PBC-whole, relabelled trajectory from tandem_unwrap.py
(FP1 -> chain A, linker -> chain B, FP2 -> chain C).
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args():
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=int, default=1280)
    ap.add_argument("--height", type=int, default=960)
    ap.add_argument("--fps", type=float, default=12.5)
    ap.add_argument("--stride", type=int, default=2,
                    help="render every Nth trajectory state (default: 2)")
    ap.add_argument("--rotate-deg", type=float, default=70.0)
    ap.add_argument("--keep-frames", action="store_true")
    return ap.parse_args(argv)


def main():
    args = parse_args()
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
    n = cmd.count_states("traj")
    print(f"[*] {n} states from {args.traj}")

    cmd.bg_color("white")
    for k, v in [("ray_opaque_background", 1), ("depth_cue", 0), ("antialias", 2),
                 ("ray_shadows", 0), ("cartoon_fancy_helices", 1),
                 ("cartoon_transparency", 0.0), ("specular", 0.25), ("ambient", 0.35),
                 ("stick_radius", 0.16), ("hash_max", 240)]:
        cmd.set(k, v)

    cmd.set("cartoon_gap_cutoff", 0)                  # no dashes across the CR2 gap
    cmd.hide("everything")
    cmd.dss("traj")
    cmd.show("cartoon", "traj")
    cmd.color("orange", "traj and chain A")           # FP1 barrel
    cmd.color("marine", "traj and chain C")           # FP2 barrel
    cmd.color("forest", "traj and chain B")           # flexible linker
    cmd.set("cartoon_loop_radius", 0.35)              # fat loops -> linker reads as a tube

    # chromophores: highlight both CR2 as sticks + spheres
    cmd.show("sticks", "traj and resn CR2")
    cmd.show("spheres", "traj and resn CR2")
    cmd.set("sphere_scale", 0.30, "traj and resn CR2")
    util.cbay("traj and resn CR2")                    # yellow carbons

    # stabilise the camera on FP1 so FP2 visibly docks/undocks
    try:
        cmd.intra_fit("traj and chain A and name CA", 1)
    except Exception as exc:
        print(f"[!] intra_fit skipped: {exc}")
    cmd.orient("traj and (chain A or chain C)")
    cmd.zoom("traj", buffer=8)

    states = list(range(1, n + 1, max(1, args.stride)))
    deg = float(args.rotate_deg) / max(1, len(states) - 1)
    pad = len(str(n))
    for frame_no, s in enumerate(states, 1):
        cmd.set("state", s)
        if frame_no > 1 and deg:
            cmd.turn("y", deg)
        cmd.ray(args.width, args.height)
        cmd.png(str(frames_dir / f"frm_{s:0{pad}d}.png"), dpi=200)
        if frame_no % 10 == 0 or frame_no == 1:
            print(f"    frame {frame_no}/{len(states)} (trajectory state {s}/{n})")

    pngs = sorted(frames_dir.glob("frm_*.png"))
    if not pngs:
        print("[ERROR] no frames rendered")
        return 1
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        r = subprocess.run(
            [ffmpeg, "-y", "-framerate", str(args.fps), "-pattern_type", "glob",
             "-i", str(frames_dir / "frm_*.png"), "-c:v", "libx264",
             "-pix_fmt", "yuv420p", "-crf", "18",
             "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(out)],
            capture_output=True, text=True)
        if r.returncode != 0:
            print("[ERROR] ffmpeg:\n" + r.stderr[-1200:])
            return 1
        print(f"[*] Wrote {out} ({len(pngs)} frames @ {args.fps} fps)")
        if not args.keep_frames:
            for p in pngs:
                p.unlink()
    else:
        print(f"[!] ffmpeg not found; frames in {frames_dir}")
    return 0


_rc = main()
try:
    from pymol import cmd as _cmd
    _cmd.quit(_rc or 0)
except Exception:
    sys.exit(_rc or 0)

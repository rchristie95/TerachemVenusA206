#!/usr/bin/env python3
r"""
render_qm_density_movie.py -- Movie of the two CR2 chromophores + transition
density during an NVT trajectory, showing the coupling geometry.

The transition density is placed onto each chromophore every frame using the
SAME code the coupling average uses (coupling_core.get_super_matrices_with_pymol
for the per-frame superposition transform + apply_pymol_matrix to move the
density points) -- so what is drawn is exactly the density that enters J(t),
snapped identically. The density is rendered as a coloured point cloud (red/blue
= +/- transition-density phase for barrel A, orange/cyan for barrel B), avoiding
the per-frame isosurface regeneration that intermittently dropped a chromophore.
Only the 29-atom CR2 part of the 44-atom STEOM QM geometry is transported as
ball-and-stick.  Tyr203 and the three QM/MM link hydrogens are intentionally not
drawn, so capped boundary fragments cannot be mistaken for free molecules.

Input: the 2-chain coupling trajectory (FP1 -> chain A, FP2 -> chain B).

    python render_qm_density_movie.py --traj tandem_nvt_1000_clean.pdb \
        --out videos/tandem_qm_density.mp4 --stride 2 --fps 12.5
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

B = "/home/robson/PetaChem"
MONO = f"{B}/tc_simple_old/classical_relaxed.pdb"                 # coupling's monomer ref
DENSITY = f"{B}/neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz"
QMXYZ = f"{B}/neo_model/orca_steom/steom_qm.xyz"
CR2_ATOM_COUNT = 29

from coupling_core import get_super_matrices_with_pymol, apply_pymol_matrix


def parse_args():
    argv = sys.argv[1:]
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--n-density", type=int, default=4000,
                    help="number of highest-|q| density points to draw per chromophore")
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--height", type=int, default=1100)
    ap.add_argument("--fps", type=float, default=12.5)
    ap.add_argument("--rotate-deg", type=float, default=90.0)
    ap.add_argument("--max-states", type=int, default=0)
    ap.add_argument("--keep-frames", action="store_true")
    return ap.parse_args()


def split_frames(traj):
    """Yield (index, pdb_text) for each MODEL of a multi-model 2-chain PDB."""
    cryst = ""
    cur, idx = [], 0
    with open(traj) as f:
        for l in f:
            if l.startswith("CRYST1"):
                cryst = l
            elif l.startswith("MODEL"):
                cur = []
            elif l.startswith("ENDMDL"):
                idx += 1
                yield idx, cryst + "".join(cur) + "END\n"
            elif l[:6] in ("ATOM  ", "HETATM"):
                cur.append(l)


def main():
    args = parse_args()
    from pymol import cmd, util
    try:
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

    # density points (highest |q| -> the visible lobes)
    d = np.load(DENSITY)
    pts, q = np.asarray(d["pts_ang"], float), np.asarray(d["q"], float)
    keep = np.argsort(np.abs(q))[-args.n_density:]
    pts, q = pts[keep], q[keep]
    pos_mask = q > 0
    print(f"[*] density: drawing {len(pts)} points ({pos_mask.sum()} +, {(~pos_mask).sum()} -)")

    # per-frame superposition matrices via the COUPLING's code
    frames = [(i, txt) for i, txt in split_frames(args.traj)]
    frames = frames[::max(1, args.stride)]
    if args.max_states:
        frames = frames[:args.max_states]
    mats = []
    with tempfile.TemporaryDirectory() as td:
        fp = Path(td) / "frame.pdb"
        for i, txt in frames:
            fp.write_text(txt)
            mA, mB, _, _, err = get_super_matrices_with_pymol(MONO, str(fp))
            mats.append(None if err else (mA, mB))
        print(f"[*] computed {sum(m is not None for m in mats)}/{len(mats)} frame transforms")

    # ---- render session ----
    cmd.reinitialize()
    for k, v in [("ray_shadows", 0), ("antialias", 2), ("depth_cue", 0),
                 ("ray_opaque_background", 1), ("stick_radius", 0.16), ("sphere_scale", 0.24),
                 ("valence", 0), ("specular", 0.2), ("ambient", 0.4), ("hash_max", 240),
                 # CRITICAL: without this, every per-frame load_cgo re-zooms the camera
                 # onto the last-loaded density object (barrel B), losing barrel A.
                 ("auto_zoom", 0)]:
        cmd.set(k, v)
    cmd.bg_color("black")
    cmd.load(QMXYZ, "qmA"); cmd.load(QMXYZ, "qmB")
    # steom_qm.xyz is ordered as the complete 29-atom CR2 followed by Tyr203
    # and three QM/MM link hydrogens.  Keep only CR2 in the displayed geometry.
    cmd.remove(f"qmA and index {CR2_ATOM_COUNT + 1}-44")
    cmd.remove(f"qmB and index {CR2_ATOM_COUNT + 1}-44")
    for q_ in ("qmA", "qmB"):
        cmd.show("sticks", q_); cmd.show("spheres", q_)
    util.cbaw("qmA or qmB"); cmd.color("yellow", "(qmA or qmB) and elem H")
    cmd.set("sphere_scale", 0.22, "qmA or qmB")
    # the QM geometry is in the anionic frame but the density is in the OLD-monomer
    # frame (align_steom_density.py). Bring the QM into the old frame with the SAME
    # Kabsch fit so it co-locates with the density.
    qm_ref = _qm_in_old_frame(np.asarray(cmd.get_coords("qmA"), float))
    cmd.load_coords(qm_ref, "qmA"); cmd.load_coords(qm_ref, "qmB")

    from pymol.cgo import COLOR, SPHERE
    RAD = 0.16
    cols = {"Ap": (0.90, 0.10, 0.10), "An": (0.15, 0.30, 0.95),
            "Bp": (1.00, 0.55, 0.00), "Bn": (0.10, 0.85, 0.90)}

    def cloud(placed, mask, rgb):
        cgo = [COLOR, *rgb]
        for p in placed[mask]:
            cgo += [SPHERE, float(p[0]), float(p[1]), float(p[2]), RAD]
        return cgo

    from align_steom_density import kabsch

    def build(mA, mB, qA1):
        # Place each chromophore + its density with the coupling's OWN single-matrix
        # transform (individually correct -- this is exactly what enters J(t)).
        qA = apply_pymol_matrix(qm_ref, mA)
        qB = apply_pymol_matrix(qm_ref, mB)
        pA = apply_pymol_matrix(pts, mA)
        pB = apply_pymol_matrix(pts, mB)
        # Stabilise the view like the other videos' `intra_fit` on chain A: rigidly
        # fit THIS frame's barrel-A QM onto the frame-1 barrel-A QM, then apply that
        # same transform to everything. Barrel A is frozen; barrel B carries its true
        # relative (hinge) motion. Done in point space so no 4x4 super-matrix
        # composition is involved (those PyMOL object matrices don't compose cleanly).
        R, t, _ = kabsch(qA, qA1)
        stab = lambda X: (R @ np.asarray(X, float).T).T + t
        qA, qB, pA, pB = stab(qA), stab(qB), stab(pA), stab(pB)
        cmd.load_coords(qA, "qmA"); cmd.load_coords(qB, "qmB")
        for nm in ("densAp", "densAn", "densBp", "densBn"):
            cmd.delete(nm)
        cmd.load_cgo(cloud(pA, pos_mask, cols["Ap"]), "densAp")
        cmd.load_cgo(cloud(pA, ~pos_mask, cols["An"]), "densAn")
        cmd.load_cgo(cloud(pB, pos_mask, cols["Bp"]), "densBp")
        cmd.load_cgo(cloud(pB, ~pos_mask, cols["Bn"]), "densBn")
        return qA.mean(0), qB.mean(0)

    valid = [(i, m) for (i, _), m in zip(frames, mats) if m is not None]
    qA1 = apply_pymol_matrix(qm_ref, valid[0][1][0])   # frame-1 barrel-A QM (stab. ref)
    deg = float(args.rotate_deg) / max(1, len(valid) - 1)
    pad = 5
    for k, (idx, (mA, mB)) in enumerate(valid):
        build(mA, mB, qA1)
        if k == 0:
            # same camera as the other videos: orient on frame 1, then a modest
            # zoom buffer -- barrel A is frozen and barrel B's hinge only shifts its
            # centroid ~1.5 A over the whole run, so both stay comfortably framed.
            cmd.orient("qmA or qmB")
            cmd.zoom("qmA or qmB", 5)
        elif deg:
            cmd.turn("y", deg)
        cmd.ray(args.width, args.height)
        cmd.png(str(frames_dir / f"frm_{idx:0{pad}d}.png"), dpi=200)
        if k % 20 == 0 or k == 0:
            print(f"    frame {k+1}/{len(valid)}")

    pngs = sorted(frames_dir.glob("frm_*.png"))
    if not pngs:
        print("[ERROR] no frames"); return 1
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        r = subprocess.run(
            [ffmpeg, "-y", "-framerate", str(args.fps), "-pattern_type", "glob",
             "-i", str(frames_dir / "frm_*.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
             "-crf", "18", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(out)],
            capture_output=True, text=True)
        if r.returncode != 0:
            print("[ERROR] ffmpeg:\n" + r.stderr[-1200:]); return 1
        print(f"[*] Wrote {out} ({len(pngs)} frames @ {args.fps} fps, "
              f"{len(pngs)/args.fps:.0f} s)")
        if not args.keep_frames:
            for p in pngs:
                p.unlink()
    return 0


def _qm_in_old_frame(qm):
    """Kabsch-fit the anionic-frame QM into the old-monomer frame (matching
    align_steom_density.py), so the QM co-locates with the old-frame density."""
    from align_steom_density import cr2_atoms, kabsch
    anion = cr2_atoms(f"{B}/tc_simple_anionic/monomer_relaxed.pdb")
    old = cr2_atoms(f"{B}/tc_simple_old/classical_relaxed.pdb")
    common = sorted(set(anion) & set(old))
    P = np.array([anion[n] for n in common])
    Q = np.array([old[n] for n in common])
    R, t, _ = kabsch(P, Q)
    return (R @ qm.T).T + t


_rc = main()
try:
    from pymol import cmd as _cmd
    _cmd.quit(_rc or 0)
except Exception:
    sys.exit(_rc or 0)

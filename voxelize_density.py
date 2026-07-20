#!/usr/bin/env python3
r"""
voxelize_density.py -- Turn the coupling's spec-normalised STEOM transition
density point cloud (steom_transdens_specnorm_oldframe.npz, already in the "old"
monomer frame the NVT dimer chains were built from) into an axis-aligned OpenDX
volumetric grid, so it can be rendered as a translucent isosurface and placed on
each trajectory barrel with the SAME super/matrix_copy mechanism the published
figures use (visualise_dimer.pml). This is the exact density that enters J(t),
just rasterised for display.

    python voxelize_density.py \
        --npz neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz \
        --out neo_model/orca_steom/steom_transdens_oldframe.dx
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz")
    ap.add_argument("--out", default="neo_model/orca_steom/steom_transdens_oldframe.dx")
    ap.add_argument("--spacing", type=float, default=0.30, help="grid spacing (Angstrom)")
    ap.add_argument("--pad", type=float, default=2.5, help="padding around points (Angstrom)")
    ap.add_argument("--smear", type=float, default=0.40, help="Gaussian smear width (Angstrom)")
    args = ap.parse_args()

    d = np.load(args.npz)
    pts, q = np.asarray(d["pts_ang"], float), np.asarray(d["q"], float)
    h = args.spacing
    lo = pts.min(0) - args.pad
    hi = pts.max(0) + args.pad
    dims = np.ceil((hi - lo) / h).astype(int) + 1
    grid = np.zeros(tuple(dims), float)

    # trilinear deposit of the signed transition charge onto the grid
    f = (pts - lo) / h
    i0 = np.floor(f).astype(int)
    fr = f - i0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = (np.where(dx, fr[:, 0], 1 - fr[:, 0]) *
                     np.where(dy, fr[:, 1], 1 - fr[:, 1]) *
                     np.where(dz, fr[:, 2], 1 - fr[:, 2]))
                idx = (np.clip(i0[:, 0] + dx, 0, dims[0] - 1),
                       np.clip(i0[:, 1] + dy, 0, dims[1] - 1),
                       np.clip(i0[:, 2] + dz, 0, dims[2] - 1))
                np.add.at(grid, idx, w * q)

    grid = gaussian_filter(grid, sigma=args.smear / h)
    grid /= np.abs(grid).max()   # normalise to max|f| = 1 (matches the cube convention)

    # ---- write OpenDX (Angstrom, z fastest = C order) ----
    nx, ny, nz = dims
    flat = grid.reshape(-1)          # C order: ((ix*ny)+iy)*nz+iz -> z fastest
    lines = [
        f"object 1 class gridpositions counts {nx} {ny} {nz}",
        f"origin {lo[0]:.6f} {lo[1]:.6f} {lo[2]:.6f}",
        f"delta {h:.6f} 0 0",
        f"delta 0 {h:.6f} 0",
        f"delta 0 0 {h:.6f}",
        f"object 2 class gridconnections counts {nx} {ny} {nz}",
        f"object 3 class array type double rank 0 items {flat.size} data follows",
    ]
    for i in range(0, flat.size, 3):
        lines.append(" ".join(f"{x:.6e}" for x in flat[i:i + 3]))
    lines.append('object "density" class field')
    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"[voxelize] wrote {args.out}: grid {nx}x{ny}x{nz} @ {h} A, "
          f"origin {lo.round(2)}, |v|max=1 (smear {args.smear} A)")


if __name__ == "__main__":
    main()

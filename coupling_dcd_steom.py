#!/usr/bin/env python3
"""Evaluate cap-masked STEOM transition-density coupling for every DCD frame.

The rigid CR2 transforms are extracted separately with extract_cr2_transforms.py.
For tractable ensemble evaluation, nearby same-sign voxels are charge-conserving
aggregated: each cell retains separate positive and negative charges at their
charge-weighted centroids.  This preserves total charge and the full transition
dipole exactly; convergence is checked by repeating representative frames at
successively finer cell sizes.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
import pyopencl as cl

BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.63

KERNEL = r"""
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void coulomb_row_sum_f64(
    __global const double *pts1,
    __global const double *q1,
    __global const double *pts2,
    __global const double *q2,
    const int n,
    __global double *out
){
    int i = get_global_id(0);
    if (i >= n) return;
    double xi = pts1[3*i], yi = pts1[3*i+1], zi = pts1[3*i+2];
    double qi = q1[i], acc = 0.0;
    for (int j = 0; j < n; ++j) {
        double dx = xi-pts2[3*j], dy = yi-pts2[3*j+1], dz = zi-pts2[3*j+2];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.1) r = 0.1;
        acc += qi*q2[j]/r;
    }
    out[i] = acc;
}

__kernel void coulomb_row_sum_f32(
    __global const float *pts1,
    __global const float *q1,
    __global const float *pts2,
    __global const float *q2,
    const int n,
    __global float *out
){
    int i = get_global_id(0);
    if (i >= n) return;
    float xi = pts1[3*i], yi = pts1[3*i+1], zi = pts1[3*i+2];
    float qi = q1[i], acc = 0.0f;
    for (int j = 0; j < n; ++j) {
        float dx = xi-pts2[3*j], dy = yi-pts2[3*j+1], dz = zi-pts2[3*j+2];
        float r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.1f) r = 0.1f;
        acc += qi*q2[j]/r;
    }
    out[i] = acc;
}
"""


def norm(vector: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.asarray(vector) ** 2)))


def rotate(points: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    """Return points @ rotation.T without dispatching a BLAS matmul."""
    p = np.asarray(points)
    return np.stack(
        (
            p[..., 0] * rotation[0, 0] + p[..., 1] * rotation[0, 1] + p[..., 2] * rotation[0, 2],
            p[..., 0] * rotation[1, 0] + p[..., 1] * rotation[1, 1] + p[..., 2] * rotation[1, 2],
            p[..., 0] * rotation[2, 0] + p[..., 1] * rotation[2, 1] + p[..., 2] * rotation[2, 2],
        ),
        axis=-1,
    )


def aggregate_density(points: np.ndarray, charges: np.ndarray, width: float):
    if width <= 0:
        return np.asarray(points, dtype=np.float64), np.asarray(charges, dtype=np.float64)
    grid_origin = points.min(axis=0) - 1.0e-9
    cell = np.floor((points - grid_origin) / width).astype(np.int32)
    output_points, output_charges = [], []
    for sign_mask in (charges > 0, charges < 0):
        keys = cell[sign_mask]
        q = charges[sign_mask]
        p = points[sign_mask]
        _, inverse = np.unique(keys, axis=0, return_inverse=True)
        q_sum = np.bincount(inverse, weights=q)
        p_sum = np.stack([np.bincount(inverse, weights=q * p[:, axis]) for axis in range(3)], axis=1)
        keep = np.abs(q_sum) > 1.0e-18
        output_charges.append(q_sum[keep])
        output_points.append(p_sum[keep] / q_sum[keep, None])
    return (
        np.ascontiguousarray(np.concatenate(output_points), dtype=np.float64),
        np.ascontiguousarray(np.concatenate(output_charges), dtype=np.float64),
    )


class CoulombOpenCL:
    def __init__(self, charges: np.ndarray, precision: str = "f64"):
        platforms = cl.get_platforms()
        devices = []
        for platform in platforms:
            devices.extend(platform.get_devices(device_type=cl.device_type.GPU))
        if not devices:
            devices = platforms[0].get_devices()
        self.device = devices[0]
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        self.program = cl.Program(self.context, KERNEL).build()
        self.dtype = np.float32 if precision == "f32" else np.float64
        self.kernel = self.program.coulomb_row_sum_f32 if precision == "f32" else self.program.coulomb_row_sum_f64
        self.precision = precision
        self.n = int(len(charges))
        self.local = min(256, int(self.device.max_work_group_size))
        self.global_size = ((self.n + self.local - 1) // self.local) * self.local
        flags = cl.mem_flags
        q = np.ascontiguousarray(charges, dtype=self.dtype)
        itemsize = int(np.dtype(self.dtype).itemsize)
        self.q_buffer = cl.Buffer(self.context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=q)
        self.a_buffer = cl.Buffer(self.context, flags.READ_ONLY, size=self.n * 3 * itemsize)
        self.b_buffer = cl.Buffer(self.context, flags.READ_ONLY, size=self.n * 3 * itemsize)
        self.out_buffer = cl.Buffer(self.context, flags.WRITE_ONLY, size=self.n * itemsize)
        self.out = np.empty(self.n, dtype=self.dtype)
        print(f"[opencl] {self.device.name}; {self.n} signed cell charges; precision={precision}")

    def calculate(self, points_a: np.ndarray, points_b: np.ndarray) -> float:
        a = np.ascontiguousarray(points_a, dtype=self.dtype).reshape(-1)
        b = np.ascontiguousarray(points_b, dtype=self.dtype).reshape(-1)
        cl.enqueue_copy(self.queue, self.a_buffer, a, is_blocking=False)
        cl.enqueue_copy(self.queue, self.b_buffer, b, is_blocking=False)
        self.kernel(
            self.queue,
            (self.global_size,),
            (self.local,),
            self.a_buffer,
            self.q_buffer,
            self.b_buffer,
            self.q_buffer,
            np.int32(self.n),
            self.out_buffer,
        )
        cl.enqueue_copy(self.queue, self.out, self.out_buffer)
        self.queue.finish()
        return float(np.sum(self.out, dtype=np.float64))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--density", type=Path, required=True)
    ap.add_argument("--transforms", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--epsilon", type=float, default=1.77)
    ap.add_argument("--bin-width", type=float, default=0.8, help="Charge-conserving cell width in Angstrom")
    ap.add_argument("--precision", choices=("f32", "f64"), default="f64",
                    help="OpenCL arithmetic precision (validate f32 against f64 before production use)")
    ap.add_argument("--start", type=int, default=0, help="First zero-based frame")
    ap.add_argument("--end", type=int, default=-1, help="Exclusive zero-based end frame")
    ap.add_argument("--stride", type=int, default=1)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    density = np.load(args.density)
    points_full = np.asarray(density["pts_ang"], dtype=np.float64)
    charges_full = np.asarray(density["q"], dtype=np.float64)
    density_origin = points_full.mean(axis=0)
    mu_local = np.sum(charges_full[:, None] * (points_full - density_origin), axis=0) * ANGSTROM_TO_BOHR

    points, charges = aggregate_density(points_full, charges_full, args.bin_width)
    mu_aggregated = np.sum(charges[:, None] * (points - density_origin), axis=0) * ANGSTROM_TO_BOHR
    charge_error = float(charges.sum() - charges_full.sum())
    dipole_error = norm(mu_aggregated - mu_local)
    print(
        f"[density] {len(points_full)} voxels -> {len(points)} signed cells at {args.bin_width:.3f} A; "
        f"DeltaQ={charge_error:+.3e}, |DeltaMu|={dipole_error:.3e} au"
    )

    transforms = np.load(args.transforms)
    rotations = np.asarray(transforms["rotation"], dtype=np.float64)
    translations = np.asarray(transforms["translation"], dtype=np.float64)
    fit_rmsd = np.asarray(transforms["rmsd_A"], dtype=np.float64)
    n_total = rotations.shape[0]
    end = n_total if args.end < 0 else min(args.end, n_total)
    frame_indices = list(range(max(0, args.start), end, max(1, args.stride)))
    if not frame_indices:
        raise ValueError("No frames selected")

    engine = CoulombOpenCL(charges, precision=args.precision)
    rows = []
    started = time.time()
    for count, frame in enumerate(frame_indices, 1):
        ra, rb = rotations[frame, 0], rotations[frame, 1]
        ta, tb = translations[frame, 0], translations[frame, 1]
        points_a = rotate(points, ra) + ta
        points_b = rotate(points, rb) + tb
        raw = engine.calculate(points_a, points_b)
        # raw has units of inverse Angstrom. Convert 1/A to 1/bohr by
        # multiplying by 0.529177, not the reciprocal coordinate factor.
        j_cm = raw * BOHR_TO_ANGSTROM * HARTREE_TO_CM / args.epsilon

        origin_a = rotate(density_origin, ra) + ta
        origin_b = rotate(density_origin, rb) + tb
        mu_a, mu_b = rotate(mu_local, ra), rotate(mu_local, rb)
        r_vec_a = origin_b - origin_a
        separation = norm(r_vec_a)
        r_vec = r_vec_a * ANGSTROM_TO_BOHR
        r = norm(r_vec)
        r_hat = r_vec / r
        mu_a_mag, mu_b_mag = norm(mu_a), norm(mu_b)
        cosine = float(np.sum(mu_a * mu_b) / (mu_a_mag * mu_b_mag + 1.0e-30))
        angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
        jdd = float(np.sum(mu_a * mu_b) - 3.0 * np.sum(mu_a * r_hat) * np.sum(mu_b * r_hat))
        j_pda_cm = jdd / (r ** 3 * args.epsilon) * HARTREE_TO_CM
        rows.append(
            {
                "frame": frame,
                "J_cm": j_cm,
                "J_pda_cm": j_pda_cm,
                "angle_deg": angle,
                "separation_A": separation,
                "aln_A_rms": float(fit_rmsd[frame, 0]),
                "aln_B_rms": float(fit_rmsd[frame, 1]),
                "mu_A": mu_a,
                "mu_B": mu_b,
                "r_A": origin_a,
                "r_B": origin_b,
            }
        )
        if count == 1 or count % 25 == 0 or count == len(frame_indices):
            elapsed = time.time() - started
            print(
                f"[coupling] {count}/{len(frame_indices)} frame {frame}: J={j_cm:.4f} cm^-1, "
                f"PDA={j_pda_cm:.4f}, R={separation:.3f} A, elapsed={elapsed:.1f}s"
            )

    j = np.asarray([row["J_cm"] for row in rows])
    stats = {
        "n": int(len(j)),
        "mean": float(j.mean()),
        "std": float(j.std(ddof=1)) if len(j) > 1 else 0.0,
        "min": float(j.min()),
        "max": float(j.max()),
        "median": float(np.median(j)),
        "epsilon": float(args.epsilon),
        "bin_width_A": float(args.bin_width),
        "opencl_precision": args.precision,
        "original_points": int(len(points_full)),
        "aggregated_signed_cells": int(len(points)),
        "charge_error_au": charge_error,
        "dipole_error_au": dipole_error,
        "dipole_norm_au": norm(mu_local),
        "fit_rmsd_mean_A": float(fit_rmsd[frame_indices].mean()),
        "fit_rmsd_max_A": float(fit_rmsd[frame_indices].max()),
        "samples": j.tolist(),
        "tdc_units": {
            "status": "corrected",
            "pair_distance_unit": "angstrom",
            "reciprocal_distance_to_atomic_units": BOHR_TO_ANGSTROM,
        },
    }
    columns = ["frame", "J_cm", "J_pda_cm", "angle_deg", "separation_A", "aln_A_rms", "aln_B_rms"]
    with open(args.out / "coupling_samples.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with open(args.out / "coupling_distribution.json", "w") as handle:
        json.dump(stats, handle, indent=2)
    np.savez(
        args.out / "coupling_geometry.npz",
        frame=np.asarray([row["frame"] for row in rows], dtype=int),
        J_cm=j,
        mu_A=np.asarray([row["mu_A"] for row in rows]),
        mu_B=np.asarray([row["mu_B"] for row in rows]),
        r_A=np.asarray([row["r_A"] for row in rows]),
        r_B=np.asarray([row["r_B"] for row in rows]),
        epsilon=float(args.epsilon),
    )
    print(
        f"[result] J = {stats['mean']:.4f} +/- {stats['std']:.4f} cm^-1 "
        f"(sample SD, n={stats['n']}); range {stats['min']:.4f}..{stats['max']:.4f}"
    )


if __name__ == "__main__":
    main()

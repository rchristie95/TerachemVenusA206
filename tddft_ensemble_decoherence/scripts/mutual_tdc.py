#!/usr/bin/env python3
"""Signed rectangular GPU TDC for independently computed site densities."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

cl = None

BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.6314

KERNEL = r"""
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void mutual_rows(
    __global const double *p1, __global const double *q1, const int n1,
    __global const double *p2, __global const double *q2, const int n2,
    __global double *out)
{
    int i = get_global_id(0);
    if (i >= n1) return;
    double xi=p1[3*i], yi=p1[3*i+1], zi=p1[3*i+2], acc=0.0;
    for (int j=0; j<n2; ++j) {
        double dx=xi-p2[3*j], dy=yi-p2[3*j+1], dz=zi-p2[3*j+2];
        double r=sqrt(dx*dx+dy*dy+dz*dz);
        acc += q1[i]*q2[j]/r;
    }
    out[i]=acc;
}
"""


class Engine:
    def __init__(self) -> None:
        global cl
        try:
            import pyopencl as opencl
        except ModuleNotFoundError as error:
            raise RuntimeError(
                "pyopencl is required for a TDC calculation; activate the "
                "validated TeraChem Conda environment"
            ) from error
        cl = opencl
        devices = [
            device for platform in cl.get_platforms()
            for device in platform.get_devices(device_type=cl.device_type.GPU)
        ]
        if not devices:
            raise RuntimeError("No OpenCL GPU found")
        self.device = devices[0]
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        self.program = cl.Program(self.context, KERNEL).build()
        self.kernel = cl.Kernel(self.program, "mutual_rows")

    def calculate(self, p1: np.ndarray, q1: np.ndarray, p2: np.ndarray, q2: np.ndarray) -> float:
        p1 = np.ascontiguousarray(p1, dtype=np.float64).reshape(-1)
        p2 = np.ascontiguousarray(p2, dtype=np.float64).reshape(-1)
        q1 = np.ascontiguousarray(q1, dtype=np.float64)
        q2 = np.ascontiguousarray(q2, dtype=np.float64)
        n1, n2 = len(q1), len(q2)
        output = np.empty(n1, dtype=np.float64)
        flags = cl.mem_flags
        buffers = [cl.Buffer(self.context, flags.READ_ONLY | flags.COPY_HOST_PTR, hostbuf=array) for array in (p1, q1, p2, q2)]
        out_buffer = cl.Buffer(self.context, flags.WRITE_ONLY, output.nbytes)
        local = min(256, int(self.device.max_work_group_size))
        global_size = ((n1 + local - 1) // local) * local
        self.kernel(
            self.queue, (global_size,), (local,), *buffers[:2], np.int32(n1),
            *buffers[2:], np.int32(n2), out_buffer,
        )
        cl.enqueue_copy(self.queue, output, out_buffer)
        self.queue.finish()
        return float(output.sum(dtype=np.float64))


def cpu_subset(p1: np.ndarray, q1: np.ndarray, p2: np.ndarray, q2: np.ndarray, count: int) -> float:
    n1, n2 = min(count, len(q1)), min(count, len(q2))
    total = 0.0
    for start in range(0, n1, 128):
        delta = p1[start : start + 128, None, :] - p2[None, :n2, :]
        distance = np.sqrt(np.sum(delta * delta, axis=2))
        total += float(np.sum(q1[start : start + 128, None] * q2[None, :n2] / distance))
    return total


def dipole(points: np.ndarray, charges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    origin = np.average(points, axis=0, weights=np.abs(charges))
    mu = np.sum(charges[:, None] * (points - origin), axis=0) * ANGSTROM_TO_BOHR
    return origin, mu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-a", type=Path, required=True)
    parser.add_argument("--site-b", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--epsilon", type=float, default=1.77)
    parser.add_argument("--cpu-subset", type=int, default=1000)
    parser.add_argument(
        "--box-lengths", nargs=3, type=float, metavar=("LX", "LY", "LZ"),
        help="Orthorhombic box; translate site B as a whole to the nearest image of site A",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    a, b = np.load(args.site_a), np.load(args.site_b)
    pa, qa = np.asarray(a["pts_ang"], float), np.asarray(a["q"], float)
    pb, qb = np.asarray(b["pts_ang"], float), np.asarray(b["q"], float)
    lattice_shift_b = np.zeros(3, dtype=float)
    if args.box_lengths is not None:
        box = np.asarray(args.box_lengths, dtype=float)
        if np.any(box <= 0.0):
            raise ValueError("Box lengths must be positive")
        center_a = np.asarray(a["physical_atom_coords_ang"], float)[:29].mean(axis=0)
        center_b = np.asarray(b["physical_atom_coords_ang"], float)[:29].mean(axis=0)
        lattice_shift_b = -box * np.round((center_b - center_a) / box)
        pb = pb + lattice_shift_b
    engine = Engine()
    started = time.time()
    raw_ab = engine.calculate(pa, qa, pb, qb)
    raw_ba = engine.calculate(pb, qb, pa, qa)
    runtime = time.time() - started
    symmetry_abs = abs(raw_ab - raw_ba)
    symmetry_rel = symmetry_abs / max(abs(raw_ab), abs(raw_ba), 1.0e-30)

    subset = max(1, args.cpu_subset)
    count_a, count_b = min(subset, len(qa)), min(subset, len(qb))
    index_a = np.argpartition(np.abs(qa), -count_a)[-count_a:]
    index_b = np.argpartition(np.abs(qb), -count_b)[-count_b:]
    pa_sub, qa_sub = pa[index_a], qa[index_a]
    pb_sub, qb_sub = pb[index_b], qb[index_b]
    cpu = cpu_subset(pa_sub, qa_sub, pb_sub, qb_sub, subset)
    gpu_subset = engine.calculate(pa_sub, qa_sub, pb_sub, qb_sub)
    cpu_gpu_abs = abs(cpu - gpu_subset)
    cpu_gpu_rel = cpu_gpu_abs / max(abs(cpu), abs(gpu_subset), 1.0e-30)

    oa, mua = dipole(pa, qa)
    ob, mub = dipole(pb, qb)
    r_ang = ob - oa
    r_bohr = r_ang * ANGSTROM_TO_BOHR
    r = float(np.linalg.norm(r_bohr))
    rhat = r_bohr / r
    pda_hartree = float((np.dot(mua, mub) - 3 * np.dot(mua, rhat) * np.dot(mub, rhat)) / r**3)
    angle = float(np.degrees(np.arccos(np.clip(np.dot(mua, mub) / (np.linalg.norm(mua) * np.linalg.norm(mub)), -1, 1))))
    # Kernel sums q_i q_j / r with r in Angstrom: a reciprocal distance
    # converts to atomic units with BOHR_TO_ANGSTROM, not ANGSTROM_TO_BOHR.
    # The latter inflates J by 1/0.529177**2 = 3.5711.
    j_vacuum = raw_ab * BOHR_TO_ANGSTROM * HARTREE_TO_CM
    result = {
        "status": "complete" if symmetry_rel < 1.0e-10 and cpu_gpu_rel < 1.0e-10 else "validation_failed",
        "opencl_device": engine.device.name.strip(),
        "site_A_points": int(len(qa)),
        "site_B_points": int(len(qb)),
        "J_vacuum_cm-1": j_vacuum,
        "J_screened_cm-1": j_vacuum / args.epsilon,
        "epsilon_opt": args.epsilon,
        "J_pda_vacuum_cm-1": pda_hartree * HARTREE_TO_CM,
        "J_pda_screened_cm-1": pda_hartree * HARTREE_TO_CM / args.epsilon,
        "separation_A": float(np.linalg.norm(r_ang)),
        "dipole_angle_deg": angle,
        "symmetry_AB_raw": raw_ab,
        "symmetry_BA_raw": raw_ba,
        "symmetry_absolute_error": symmetry_abs,
        "symmetry_relative_error": symmetry_rel,
        "cpu_subset_points": int(min(subset, len(qa), len(qb))),
        "cpu_subset_raw": cpu,
        "gpu_subset_raw": gpu_subset,
        "cpu_gpu_absolute_error": cpu_gpu_abs,
        "cpu_gpu_relative_error": cpu_gpu_rel,
        "runtime_seconds": runtime,
        "site_B_lattice_shift_A": lattice_shift_b.tolist(),
        "box_lengths_A": list(args.box_lengths) if args.box_lengths is not None else None,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    with (args.out.parent / "coupling_samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)
    print(json.dumps(result, indent=2))
    if result["status"] != "complete":
        raise SystemExit(2)


if __name__ == "__main__":
    main()

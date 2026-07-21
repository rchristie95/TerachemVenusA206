#!/usr/bin/env python3
"""Estimate tandem-Venus solvation and pure dephasing from a QM/MM trajectory.

This is a deliberately lightweight numerical test.  A quantitative STEOM
excited-minus-ground atom-centred charge probe is rigidly fitted to each CR2
site in every frame.  Its Coulomb interaction with the explicit AMBER/TIP3P
MM charges gives an environmental contribution to each vertical site energy.
The two site-energy traces are then converted into solvation response
functions and a classical second-cumulant estimate of the inter-site
coherence decay.

Limitations: the probe is fixed, the environment is non-polarizable, the
chromophore/Tyr atoms belonging to the probed site are excluded, and no
intrachromophore excited-state energy fluctuation is included.  The result is
therefore an electrostatic nuclear-solvation/pure-dephasing test, not a final
quantum spectral density.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pyopencl as cl
from openmm import NonbondedForce, unit
from openmm.app import ForceField, HBonds, PME, PDBFile
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit

import run_nvt


ANGSTROM_TO_BOHR = 1.0 / 0.529177210903
HARTREE_TO_CM = 219474.6313632
LIGHT_CM_PER_S = 2.99792458e10

WATER_NAMES = {"HOH", "WAT", "TIP3", "TIP3P", "SOL"}
ION_NAMES = {
    "NA", "CL", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "F", "BR", "I",
    "Na+", "Cl-", "K+", "Ca2+", "Mg2+",
}

KERNEL = r"""
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define CHUNK 256

__kernel void electrostatic_chunks_f32(
    __global const float *probe_xyz,
    __global const float *probe_q,
    const int nprobe,
    __global const float *mm_xyz,
    __global const float *mm_q,
    __global const char *group_code,
    __global const char *exclude_code,
    const int nmm,
    const int nchunks,
    const int site,
    const float lx,
    const float ly,
    const float lz,
    const float cutoff,
    __global float *out
){
    int gid = get_global_id(0);
    int total = nprobe*nchunks;
    if (gid >= total) return;
    int ip = gid/nchunks;
    int chunk = gid - ip*nchunks;
    int start = chunk*CHUNK;
    int stop = min(start + CHUNK, nmm);
    float px = probe_xyz[3*ip];
    float py = probe_xyz[3*ip+1];
    float pz = probe_xyz[3*ip+2];
    float qp = probe_q[ip];
    float a0 = 0.0f, a1 = 0.0f, a2 = 0.0f;
    for (int j=start; j<stop; ++j) {
        if ((int)exclude_code[j] == site) continue;
        float dx = px-mm_xyz[3*j];
        float dy = py-mm_xyz[3*j+1];
        float dz = pz-mm_xyz[3*j+2];
        dx -= rint(dx/lx)*lx;
        dy -= rint(dy/ly)*ly;
        dz -= rint(dz/lz)*lz;
        float r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.25f) continue;
        if (cutoff > 0.0f && r > cutoff) continue;
        float value = qp*mm_q[j]/r;
        int g = (int)group_code[j];
        if (g == 1) a1 += value;
        else if (g == 2) a2 += value;
        else a0 += value;
    }
    int base = 3*gid;
    out[base] = a0;
    out[base+1] = a1;
    out[base+2] = a2;
}
"""


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    source_c = source - source_center
    target_c = target - target_center
    u, _, vt = np.linalg.svd(source_c.T @ target_c)
    handedness = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, handedness]) @ u.T
    translation = target_center - rotation @ source_center
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def build_exact_charges(
    topology_pdb: Path,
    generic_xml: Path,
    amber_cr2_prmtop: Path,
) -> tuple[np.ndarray, list, object]:
    pdb = PDBFile(str(topology_pdb))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml", str(generic_xml))
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=run_nvt.choose_cutoff_from_box(pdb.topology, default_nm=1.0),
        constraints=HBonds,
        ignoreExternalBonds=True,
    )
    run_nvt.apply_amber_cr2_parameters(system, pdb.topology, amber_cr2_prmtop)
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    charges = np.empty(system.getNumParticles(), dtype=np.float32)
    for index in range(system.getNumParticles()):
        charge, _, _ = nonbonded.getParticleParameters(index)
        charges[index] = charge.value_in_unit(unit.elementary_charge)
    return charges, list(pdb.topology.atoms()), pdb.topology


def classify_atoms(openmm_atoms: list) -> tuple[np.ndarray, np.ndarray, dict]:
    group_code = np.zeros(len(openmm_atoms), dtype=np.int8)
    exclude_code = np.full(len(openmm_atoms), -1, dtype=np.int8)
    residue_rows = {}
    cr2_residues = []
    for atom in openmm_atoms:
        residue = atom.residue
        key = ((residue.chain.id or "").strip(), str(residue.id).strip(), residue.name)
        residue_rows.setdefault(key, []).append(atom.index)
        name_upper = residue.name.upper()
        if name_upper in WATER_NAMES:
            group_code[atom.index] = 1
        elif name_upper in ION_NAMES:
            group_code[atom.index] = 2
        if residue.name == "CR2" and key not in cr2_residues:
            cr2_residues.append(key)
    if len(cr2_residues) != 2:
        raise RuntimeError(f"Expected two CR2 sites, found {cr2_residues}")

    # The tandem topology merges each three-residue chromophore into CR2, so
    # paper Tyr203 appears as residue 202 and its second-repeat partner as 464.
    qm_keys = []
    for site, cr2_key in enumerate(cr2_residues):
        chain, resid, _ = cr2_key
        target_tyr = "202" if site == 0 else "464"
        selected = [cr2_key]
        tyr_candidates = [key for key in residue_rows if key[0] == chain and key[1] == target_tyr and key[2] == "TYR"]
        if len(tyr_candidates) != 1:
            raise RuntimeError(f"Could not identify stacked Tyr for site {site}: {tyr_candidates}")
        selected.extend(tyr_candidates)
        for key in selected:
            exclude_code[residue_rows[key]] = site
        qm_keys.append(selected)
    return group_code, exclude_code, {"cr2_residues": cr2_residues, "excluded_qm_residues": qm_keys}


class ElectrostaticEngine:
    def __init__(
        self,
        probe_q: np.ndarray,
        mm_q: np.ndarray,
        group_code: np.ndarray,
        exclude_code: np.ndarray,
        electrostatic_cutoff_a: float = 0.0,
    ) -> None:
        platforms = cl.get_platforms()
        devices = [
            device
            for platform in platforms
            for device in platform.get_devices(device_type=cl.device_type.GPU)
        ]
        if not devices:
            raise RuntimeError("No OpenCL GPU available")
        self.device = devices[0]
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context)
        self.program = cl.Program(self.context, KERNEL).build()
        self.kernel = cl.Kernel(self.program, "electrostatic_chunks_f32")
        flags = cl.mem_flags
        self.nprobe = int(len(probe_q))
        self.nmm = int(len(mm_q))
        self.nchunks = (self.nmm + 255) // 256
        self.global_size = self.nprobe * self.nchunks
        self.electrostatic_cutoff_a = float(electrostatic_cutoff_a)
        self.probe_q_buffer = cl.Buffer(
            self.context, flags.READ_ONLY | flags.COPY_HOST_PTR,
            hostbuf=np.ascontiguousarray(probe_q, dtype=np.float32),
        )
        self.mm_q_buffer = cl.Buffer(
            self.context, flags.READ_ONLY | flags.COPY_HOST_PTR,
            hostbuf=np.ascontiguousarray(mm_q, dtype=np.float32),
        )
        self.group_buffer = cl.Buffer(
            self.context, flags.READ_ONLY | flags.COPY_HOST_PTR,
            hostbuf=np.ascontiguousarray(group_code, dtype=np.int8),
        )
        self.exclude_buffer = cl.Buffer(
            self.context, flags.READ_ONLY | flags.COPY_HOST_PTR,
            hostbuf=np.ascontiguousarray(exclude_code, dtype=np.int8),
        )
        self.probe_xyz_buffer = cl.Buffer(self.context, flags.READ_ONLY, size=self.nprobe * 3 * 4)
        self.mm_xyz_buffer = cl.Buffer(self.context, flags.READ_ONLY, size=self.nmm * 3 * 4)
        self.out = np.empty((self.global_size, 3), dtype=np.float32)
        self.out_buffer = cl.Buffer(self.context, flags.WRITE_ONLY, size=self.out.nbytes)

    def set_mm_positions(self, positions_ang: np.ndarray) -> None:
        cl.enqueue_copy(
            self.queue,
            self.mm_xyz_buffer,
            np.ascontiguousarray(positions_ang, dtype=np.float32).reshape(-1),
            is_blocking=False,
        )

    def calculate(self, probe_xyz_ang: np.ndarray, box_ang: np.ndarray, site: int) -> np.ndarray:
        cl.enqueue_copy(
            self.queue,
            self.probe_xyz_buffer,
            np.ascontiguousarray(probe_xyz_ang, dtype=np.float32).reshape(-1),
            is_blocking=False,
        )
        self.kernel(
            self.queue,
            (self.global_size,),
            None,
            self.probe_xyz_buffer,
            self.probe_q_buffer,
            np.int32(self.nprobe),
            self.mm_xyz_buffer,
            self.mm_q_buffer,
            self.group_buffer,
            self.exclude_buffer,
            np.int32(self.nmm),
            np.int32(self.nchunks),
            np.int32(site),
            np.float32(box_ang[0]),
            np.float32(box_ang[1]),
            np.float32(box_ang[2]),
            np.float32(self.electrostatic_cutoff_a),
            self.out_buffer,
        )
        cl.enqueue_copy(self.queue, self.out, self.out_buffer)
        self.queue.finish()
        # Kernel uses 1/r_A. Convert to 1/r_bohr, then Hartree to cm^-1.
        return self.out.sum(axis=0, dtype=np.float64) / ANGSTROM_TO_BOHR * HARTREE_TO_CM


def correlation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    n = len(a)
    return np.correlate(a, b, mode="full")[n - 1 :] / np.arange(n, 0, -1)


def first_crossing(time_fs: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    hits = np.where(values <= threshold)[0]
    if hits.size == 0:
        return None
    index = int(hits[0])
    if index == 0:
        return float(time_fs[0])
    x0, x1 = time_fs[index - 1], time_fs[index]
    y0, y1 = values[index - 1], values[index]
    if y1 == y0:
        return float(x1)
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def positive_integral_time(time_fs: np.ndarray, normalized: np.ndarray) -> float:
    zeros = np.where(normalized <= 0.0)[0]
    stop = int(zeros[0]) + 1 if zeros.size else len(normalized)
    return float(np.trapezoid(normalized[:stop], time_fs[:stop]))


def biexponential_fit(time_fs: np.ndarray, normalized: np.ndarray, max_fs: float = 1000.0) -> dict | None:
    mask = (time_fs <= max_fs) & np.isfinite(normalized)
    t = time_fs[mask]
    y = normalized[mask]
    if len(t) < 12:
        return None

    def model(x, amplitude, tau_fast, tau_slow, offset):
        return offset + (1.0 - offset) * (
            amplitude * np.exp(-x / tau_fast)
            + (1.0 - amplitude) * np.exp(-x / tau_slow)
        )

    try:
        params, _ = curve_fit(
            model,
            t,
            y,
            p0=(0.6, 40.0, 400.0, 0.0),
            bounds=((0.0, 1.0, 10.0, -0.5), (1.0, 400.0, 5000.0, 0.5)),
            maxfev=30000,
        )
    except Exception:
        return None
    fitted = model(t, *params)
    rmse = float(np.sqrt(np.mean((fitted - y) ** 2)))
    return {
        "fast_amplitude": float(params[0]),
        "tau_fast_fs": float(params[1]),
        "tau_slow_fs": float(params[2]),
        "offset": float(params[3]),
        "fit_rmse": rmse,
        "fit_window_fs": float(max_fs),
    }


def coherence_from_correlation(cdiff_cm2: np.ndarray, dt_fs: float) -> tuple[np.ndarray, np.ndarray]:
    c_omega = (2.0 * np.pi * LIGHT_CM_PER_S) ** 2 * np.asarray(cdiff_cm2, dtype=float)
    time_s = np.arange(len(c_omega), dtype=float) * dt_fs * 1.0e-15
    first_integral = cumulative_trapezoid(c_omega, time_s, initial=0.0)
    g = cumulative_trapezoid(first_integral, time_s, initial=0.0)
    # Sampling noise can make the long-time cumulant non-monotonic; the first
    # decay is the meaningful quantity in this short test.
    coherence = np.exp(-np.clip(g, 0.0, 700.0))
    return g, coherence


def analyze_series(time_fs: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[dict, dict[str, np.ndarray]]:
    caa = correlation(a, a)
    cbb = correlation(b, b)
    cab = correlation(a, b)
    cba = correlation(b, a)
    difference = a - b
    cdiff = correlation(difference, difference)
    sa = caa / caa[0]
    sb = cbb / cbb[0]
    sdiff = cdiff / cdiff[0]
    lag_fs = time_fs - time_fs[0]
    g, coherence = coherence_from_correlation(cdiff, float(np.median(np.diff(time_fs))))
    rho0 = float(0.5 * (cab[0] + cba[0]) / np.sqrt(caa[0] * cbb[0]))
    summary = {
        "site_A_mean_cm-1": float(np.mean(a)),
        "site_B_mean_cm-1": float(np.mean(b)),
        "site_A_sigma_cm-1": float(np.std(a, ddof=1)),
        "site_B_sigma_cm-1": float(np.std(b, ddof=1)),
        "difference_sigma_cm-1": float(np.std(difference, ddof=1)),
        "site_cross_correlation_rho0": rho0,
        "site_A_solvation_1e_fs": first_crossing(lag_fs, sa, 1.0 / np.e),
        "site_B_solvation_1e_fs": first_crossing(lag_fs, sb, 1.0 / np.e),
        "difference_correlation_1e_fs": first_crossing(lag_fs, sdiff, 1.0 / np.e),
        "site_A_positive_integral_fs": positive_integral_time(lag_fs, sa),
        "site_B_positive_integral_fs": positive_integral_time(lag_fs, sb),
        "difference_positive_integral_fs": positive_integral_time(lag_fs, sdiff),
        "classical_cumulant_T2_1e_fs": first_crossing(lag_fs, coherence, 1.0 / np.e),
        "site_A_biexponential": biexponential_fit(lag_fs, sa),
        "site_B_biexponential": biexponential_fit(lag_fs, sb),
        "difference_biexponential": biexponential_fit(lag_fs, sdiff),
    }
    arrays = {
        "lag_fs": lag_fs,
        "C_AA_cm-2": caa,
        "C_BB_cm-2": cbb,
        "C_AB_cm-2": cab,
        "C_BA_cm-2": cba,
        "C_difference_cm-2": cdiff,
        "S_A": sa,
        "S_B": sb,
        "S_difference": sdiff,
        "g_classical": g,
        "coherence_classical": coherence,
    }
    return summary, arrays


def write_csv(path: Path, columns: dict[str, np.ndarray]) -> None:
    names = list(columns)
    rows = zip(*(np.asarray(columns[name]).tolist() for name in names))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(rows)


def make_plot(path: Path, time_fs: np.ndarray, energies: np.ndarray, arrays: dict[str, np.ndarray]) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    axes[0, 0].plot(time_fs, energies[:, 0, :].sum(axis=1), label="site A", lw=1)
    axes[0, 0].plot(time_fs, energies[:, 1, :].sum(axis=1), label="site B", lw=1)
    axes[0, 0].set(xlabel="time (fs)", ylabel=r"environmental gap shift (cm$^{-1}$)")
    axes[0, 0].legend()

    lag = arrays["lag_fs"]
    axes[0, 1].plot(lag, arrays["S_A"], label="A")
    axes[0, 1].plot(lag, arrays["S_B"], label="B")
    axes[0, 1].plot(lag, arrays["S_difference"], label="A-B", alpha=0.8)
    axes[0, 1].axhline(0, color="0.6", lw=0.7)
    axes[0, 1].set(xlim=(0, min(1000, lag[-1])), ylim=(-0.5, 1.05), xlabel="lag (fs)", ylabel="normalized correlation")
    axes[0, 1].legend()

    labels = ["protein", "water", "ions"]
    for group, label in enumerate(labels):
        differential = energies[:, 0, group] - energies[:, 1, group]
        axes[1, 0].plot(time_fs, differential - differential.mean(), label=label, lw=1)
    axes[1, 0].set(xlabel="time (fs)", ylabel=r"differential fluctuation (cm$^{-1}$)")
    axes[1, 0].legend()

    axes[1, 1].plot(lag, arrays["coherence_classical"], color="black")
    axes[1, 1].axhline(1 / np.e, color="tab:red", ls="--", lw=0.8)
    axes[1, 1].set(xlim=(0, min(500, lag[-1])), ylim=(0, 1.02), xlabel="time (fs)", ylabel=r"$|L(t)|$ (classical cumulant)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--generic-xml", type=Path, required=True)
    parser.add_argument("--amber-cr2-prmtop", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--dt-fs", type=float, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=-1)
    parser.add_argument(
        "--electrostatic-cutoff-a",
        type=float,
        default=0.0,
        help="Optional real-space probe-environment cutoff in Angstrom (0 uses the full minimum-image box)",
    )
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    probe = np.load(args.probe)
    probe_coords = np.asarray(probe["atom_coords_ang"], dtype=float)
    probe_q = np.asarray(probe["atom_delta_q_e"], dtype=float)
    atomic_numbers = np.asarray(probe["atomic_numbers"], dtype=int)
    cr2_names = [str(value) for value in probe["cr2_atom_names"]]
    cr2_count = int(probe["cr2_atom_count"])
    heavy_names = [cr2_names[index] for index in range(cr2_count) if atomic_numbers[index] != 1]
    heavy_source = np.asarray(
        [probe_coords[cr2_names.index(name)] for name in heavy_names],
        dtype=float,
    )

    mm_q, openmm_atoms, _ = build_exact_charges(
        args.topology, args.generic_xml, args.amber_cr2_prmtop
    )
    group_code, exclude_code, topology_audit = classify_atoms(openmm_atoms)
    engine = ElectrostaticEngine(
        probe_q,
        mm_q,
        group_code,
        exclude_code,
        electrostatic_cutoff_a=args.electrostatic_cutoff_a,
    )

    md_topology = md.load_topology(str(args.topology))
    cr2_sites = [residue for residue in md_topology.residues if residue.name == "CR2"]
    if len(cr2_sites) != 2:
        raise RuntimeError(f"Expected two CR2 sites in MDTraj topology, found {len(cr2_sites)}")
    site_indices = []
    for residue in cr2_sites:
        atom_map = {atom.name: atom.index for atom in residue.atoms}
        missing = set(heavy_names) - set(atom_map)
        if missing:
            raise RuntimeError(f"Missing CR2 fitting atoms at {residue}: {sorted(missing)}")
        site_indices.append(np.asarray([atom_map[name] for name in heavy_names], dtype=int))

    energy_rows = []
    rmsd_rows = []
    frame_count = 0
    for chunk in md.iterload(str(args.trajectory), top=str(args.topology), chunk=25, stride=max(1, args.stride)):
        xyz_ang = np.asarray(chunk.xyz, dtype=np.float32) * 10.0
        if chunk.unitcell_lengths is None:
            raise RuntimeError("Trajectory lacks periodic unit-cell lengths")
        boxes_ang = np.asarray(chunk.unitcell_lengths, dtype=np.float32) * 10.0
        for local_frame in range(chunk.n_frames):
            if args.max_frames >= 0 and frame_count >= args.max_frames:
                break
            positions = xyz_ang[local_frame]
            engine.set_mm_positions(positions)
            frame_energies = np.empty((2, 3), dtype=float)
            frame_rmsd = np.empty(2, dtype=float)
            for site in range(2):
                target = positions[site_indices[site]]
                rotation, translation, rmsd = kabsch(heavy_source, target)
                placed_probe = probe_coords @ rotation.T + translation
                frame_energies[site] = engine.calculate(placed_probe, boxes_ang[local_frame], site)
                frame_rmsd[site] = rmsd
            energy_rows.append(frame_energies)
            rmsd_rows.append(frame_rmsd)
            frame_count += 1
            if frame_count == 1 or frame_count % 100 == 0:
                print(f"[gap] {frame_count} frames")
        if args.max_frames >= 0 and frame_count >= args.max_frames:
            break

    energies = np.asarray(energy_rows, dtype=float)
    rmsd = np.asarray(rmsd_rows, dtype=float)
    time_fs = np.arange(frame_count, dtype=float) * args.dt_fs * max(1, args.stride)
    total_a = energies[:, 0, :].sum(axis=1)
    total_b = energies[:, 1, :].sum(axis=1)
    total_summary, arrays = analyze_series(time_fs, total_a, total_b)

    group_names = ["protein", "water", "ions"]
    group_summaries = {}
    for group, name in enumerate(group_names):
        group_summaries[name], _ = analyze_series(
            time_fs, energies[:, 0, group], energies[:, 1, group]
        )

    summary = {
        "method": "fixed STEOM multipair difference-density probe against explicit AMBER/TIP3P MM charges",
        "trajectory": f"external/{args.trajectory.name}",
        "frames": frame_count,
        "dt_fs": float(args.dt_fs * max(1, args.stride)),
        "duration_fs": float(time_fs[-1] - time_fs[0]) if frame_count > 1 else 0.0,
        "opencl_device": engine.device.name.strip(),
        "electrostatic_cutoff_A": float(args.electrostatic_cutoff_a),
        "probe_net_delta_charge_e": float(probe_q.sum()),
        "probe_delta_mu_e_bohr": np.asarray(probe["delta_mu_atom_e_bohr"]).tolist(),
        "topology": topology_audit,
        "charge_counts": {
            "total_atoms": int(len(mm_q)),
            "protein_atoms": int(np.sum(group_code == 0)),
            "water_atoms": int(np.sum(group_code == 1)),
            "ion_atoms": int(np.sum(group_code == 2)),
            "net_system_charge_e": float(mm_q.sum(dtype=np.float64)),
        },
        "fit_rmsd_A": {
            "site_A_mean": float(rmsd[:, 0].mean()),
            "site_B_mean": float(rmsd[:, 1].mean()),
            "maximum": float(rmsd.max()),
        },
        "total": total_summary,
        "components": group_summaries,
        "limitations": [
            "fixed STEOM difference density; no geometry-dependent QM excitation calculation",
            "TIP3P/AMBER fixed-charge environment; no induced electronic polarization",
            "classical high-temperature second cumulant; no quantum correction or detailed balance",
            "short trajectory intended for ultrafast numerical testing, not converged slow protein dynamics",
        ],
    }

    np.savez_compressed(
        args.out / "gap_timeseries.npz",
        time_fs=time_fs,
        energies_cm=energies,
        fit_rmsd_A=rmsd,
        group_names=np.asarray(group_names),
    )
    write_csv(
        args.out / "gap_timeseries.csv",
        {
            "time_fs": time_fs,
            **{
                f"site_{site_name}_{group_name}_cm-1": energies[:, site, group]
                for site, site_name in enumerate(("A", "B"))
                for group, group_name in enumerate(group_names)
            },
            "site_A_total_cm-1": total_a,
            "site_B_total_cm-1": total_b,
            "difference_total_cm-1": total_a - total_b,
            "fit_rmsd_A_site_A": rmsd[:, 0],
            "fit_rmsd_A_site_B": rmsd[:, 1],
        },
    )
    write_csv(args.out / "correlations.csv", arrays)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    make_plot(args.out / "solvation_decoherence_test.png", time_fs, energies, arrays)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

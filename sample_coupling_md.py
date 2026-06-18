#!/usr/bin/env python3
"""
sample_coupling_md.py  --  Conformational sampling of the Davydov coupling J.

Reviewer item 1 (the criticism that sank the JPCL submission): the reported
J = 74.38 cm^-1 comes from a single minimised geometry. A protein environment
makes a single frame indefensible. This script samples J over an ensemble of
MD snapshots and reports J as mean +/- std with a histogram, converting the
weakness into a result (the spread is the static + dynamic disorder that feeds
the dephasing / lineshape analysis in lineshape_cd.py).

Two modes (configurable via --mode):

  rigid  (default, fast):
        Run TDDFT once (on the relaxed geometry; reuse an existing transition
        density). For every MD frame, re-derive the monomer->siteA / siteB
        PyMOL `super` transforms from THAT frame's chain A / chain B
        coordinates, re-map the fixed transition density onto both sites, and
        recompute J with the existing GPU TDC routine. Captures the
        orientational / positional disorder of the dimer. No per-frame QM.

  full   (rigorous, GPU + TeraChem heavy):
        For a --subset of frames, re-run the whole QM/MM TDDFT->TDC pipeline
        (terachem_full_pipeline.py) on the frame geometry. Used to validate the
        rigid approximation on a handful of frames. Requires TeraChem + OpenMM.

All heavy numerics reuse coupling_core.py (read_dx, calculate_coupling,
get_super_matrices_with_pymol, apply_pymol_matrix, transition_dipole_au, ...),
identical to terachem_full_pipeline.py Stage 3, so a single-frame run of this
script reproduces the pipeline's J.

Outputs (in --out, default `coupling_sampling_out/`):
  - coupling_samples.csv         per-frame J, J_PDA, dipole angle, separation
  - coupling_distribution.json   {mean, std, n, samples, ...}  (read by lineshape_cd.py)
  - Fig_Coupling_Histogram.pdf   histogram of J with mean +/- std

Example:
  python sample_coupling_md.py --traj tc_dimer_nvt/dimer_nvt_trajectory.pdb \\
      --workdir tc_tddft_old_current --monomer tc_simple_old/classical_relaxed.pdb \\
      --n-frames 100 --mode rigid
"""

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from coupling_core import (
    HARTREE_TO_CM,
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    read_dx,
    calculate_coupling,
    get_super_matrices_with_pymol,
    apply_pymol_matrix,
    transition_dipole_au,
    oscillator_to_dipole_au,
    autodetect_workdir_and_candidates,
    select_target_state_and_density,
    print_excited_state_table,
)


# --------------------------------------------------------------------------- #
# Trajectory handling
# --------------------------------------------------------------------------- #
def split_pdb_models(traj_path):
    """
    Split a multi-MODEL PDB trajectory into a list of frame line-blocks.

    Each returned element is a list of PDB lines (ATOM/HETATM/TER/etc.) for one
    MODEL. If the file has no MODEL records it is treated as a single frame.
    """
    frames = []
    current = []
    saw_model = False
    with open(traj_path, "r") as f:
        for line in f:
            tag = line[:6].strip()
            if tag == "MODEL":
                saw_model = True
                current = []
                continue
            if tag == "ENDMDL":
                if current:
                    frames.append(current)
                current = []
                continue
            if tag in ("ATOM", "HETATM", "TER", "HETERO", "ANISOU"):
                current.append(line)
    if not saw_model:
        # Single-structure PDB: one frame.
        return [current] if current else []
    if current:  # trailing model without ENDMDL
        frames.append(current)
    return frames


def select_frame_indices(n_total, n_frames, stride, seed, randomize):
    """Choose which frame indices to evaluate."""
    all_idx = list(range(0, n_total, max(1, stride)))
    if n_frames is None or n_frames >= len(all_idx):
        return all_idx
    if randomize:
        rng = np.random.default_rng(seed)
        return sorted(rng.choice(all_idx, size=n_frames, replace=False).tolist())
    # Evenly spaced subset.
    picks = np.linspace(0, len(all_idx) - 1, n_frames).round().astype(int)
    return sorted({all_idx[i] for i in picks})


def write_frame_pdb(frame_lines, path):
    with open(path, "w") as f:
        f.writelines(frame_lines)
        f.write("END\n")


def frame_has_chains(frame_lines, chains=("A", "B")):
    seen = {line[21] for line in frame_lines if line[:6].strip() in ("ATOM", "HETATM")}
    return all(c in seen for c in chains)


# --------------------------------------------------------------------------- #
# Coupling for one frame (rigid mode)
# --------------------------------------------------------------------------- #
def coupling_for_frame(frame_dimer_pdb, monomer_pdb, pts_opt, q_opt, epsilon, backend_kwargs):
    """
    Compute J (TDC) and the point-dipole estimate J_PDA for one dimer frame,
    re-using the fixed monomer transition density (rigid-density approximation).

    Mirrors terachem_full_pipeline.stage3_main lines ~3275-3331.

    Returns a dict of per-frame observables, or None if the PyMOL alignment failed.
    """
    matrix_A, matrix_B, aln_A, aln_B, err = get_super_matrices_with_pymol(monomer_pdb, frame_dimer_pdb)
    if err:
        return {"error": err}

    pts_A = apply_pymol_matrix(pts_opt, matrix_A)
    pts_B = apply_pymol_matrix(pts_opt, matrix_B)

    origin_A = np.mean(pts_A, axis=0)
    origin_B = np.mean(pts_B, axis=0)
    separation = float(np.linalg.norm(origin_A - origin_B))

    # Full transition-density coupling.
    J_ha = calculate_coupling(pts_A, q_opt, pts_B, q_opt, **backend_kwargs) / epsilon
    J_cm = J_ha * HARTREE_TO_CM

    # Point-dipole approximation (far-field), same construction as Stage 3.
    muA = transition_dipole_au(pts_A, q_opt, origin_angstrom=origin_A)
    muB = transition_dipole_au(pts_B, q_opt, origin_angstrom=origin_B)
    muA_mag = float(np.linalg.norm(muA))
    muB_mag = float(np.linalg.norm(muB))
    cosang = np.dot(muA, muB) / (muA_mag * muB_mag + 1e-30)
    angle = float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))

    Rvec = (origin_B - origin_A) * ANGSTROM_TO_BOHR
    R = np.linalg.norm(Rvec)
    Rhat = Rvec / (R + 1e-30)
    Jdd = np.dot(muA, muB) - 3.0 * np.dot(muA, Rhat) * np.dot(muB, Rhat)
    Vdd = Jdd / (R ** 3 * epsilon)
    J_pda_cm = float(Vdd * HARTREE_TO_CM)

    return {
        "J_cm": float(J_cm),
        "J_pda_cm": J_pda_cm,
        "angle_deg": angle,
        "separation_A": separation,
        "aln_A_rms": float(aln_A[0]) if aln_A else float("nan"),
        "aln_B_rms": float(aln_B[0]) if aln_B else float("nan"),
    }


# --------------------------------------------------------------------------- #
# Full-pipeline mode (subprocess per frame)
# --------------------------------------------------------------------------- #
J_LINE_RE = re.compile(r"J:\s*[-\d.eE+]+\s*Hartree\s*\(\s*([-\d.eE+]+)\s*cm")


def coupling_for_frame_full(frame_dimer_pdb, out_dir, pipeline_args, python_exe):
    """
    Run the complete QM/MM TDDFT->TDC pipeline on a single dimer frame by
    invoking terachem_full_pipeline.py as a subprocess in its own working dir,
    then parse the printed `J: ... cm^-1` line. Returns J in cm^-1 or None.
    """
    cmd = [
        python_exe,
        str(Path(__file__).with_name("terachem_full_pipeline.py")),
        "--cwd", str(out_dir),
    ]
    cmd += pipeline_args
    # The pipeline's Stage 1 reads --pdb via --simple-args; expose the frame.
    cmd += ["--simple-args", f"--pdb {frame_dimer_pdb}"]
    print(f"    [full] running: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout[-2000:])
    if proc.returncode != 0:
        print(f"    [full] pipeline failed (rc={proc.returncode}):\n{proc.stderr[-1500:]}")
        return None
    matches = J_LINE_RE.findall(proc.stdout)
    if not matches:
        print("    [full] could not parse J from pipeline output.")
        return None
    return float(matches[-1])


# --------------------------------------------------------------------------- #
# Aggregation & plotting
# --------------------------------------------------------------------------- #
def summarize(j_values):
    arr = np.asarray([v for v in j_values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return {"n": 0, "mean": float("nan"), "std": float("nan")}
    return {
        "n": int(arr.size),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "median": float(np.median(arr)),
    }


def plot_histogram(j_values, stats, out_pdf, single_frame_ref=74.38):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    arr = np.asarray([v for v in j_values if np.isfinite(v)], dtype=float)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    nbins = max(8, min(40, int(np.sqrt(arr.size)) * 2)) if arr.size else 8
    ax.hist(arr, bins=nbins, color="#4C72B0", alpha=0.78, edgecolor="white")
    mean, std = stats["mean"], stats["std"]
    ax.axvline(mean, color="k", lw=2, label=fr"mean $= {mean:.1f}$ cm$^{{-1}}$")
    ax.axvspan(mean - std, mean + std, color="k", alpha=0.10,
               label=fr"$\pm\sigma = {std:.1f}$ cm$^{{-1}}$")
    if single_frame_ref is not None:
        ax.axvline(single_frame_ref, color="#C44E52", ls="--", lw=2,
                   label=fr"single frame $= {single_frame_ref:.2f}$ cm$^{{-1}}$")
    ax.set_xlabel(r"Davydov coupling $J$ (cm$^{-1}$)")
    ax.set_ylabel("MD snapshots")
    ax.set_title(fr"$J = {mean:.1f} \pm {std:.1f}$ cm$^{{-1}}$ ($n={stats['n']}$)")
    ax.legend(fontsize=9, frameon=True)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def write_csv(rows, path):
    cols = ["frame", "J_cm", "J_pda_cm", "angle_deg", "separation_A", "aln_A_rms", "aln_B_rms"]
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in cols) + "\n")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", type=Path, required=True,
                   help="Multi-MODEL dimer trajectory PDB (from run_dimer_nvt.py).")
    p.add_argument("--monomer", type=Path, default=Path("tc_simple_old/classical_relaxed.pdb"),
                   help="Monomer reference PDB (chain A) for the super alignment.")
    p.add_argument("--workdir", type=Path, default=Path("tc_tddft_old_current"),
                   help="TDDFT workdir to pull the transition density / state from.")
    p.add_argument("--density", type=Path, default=None,
                   help="Explicit transition density .dx (overrides --workdir autodetect).")
    p.add_argument("--root", type=int, default=None, help="Force a specific excited-state root.")
    p.add_argument("--density-mode", choices=["signed", "auto", "abs"], default="auto")
    p.add_argument("--mode", choices=["rigid", "full"], default="rigid")
    p.add_argument("--n-frames", type=int, default=100, help="Number of frames to sample (rigid).")
    p.add_argument("--subset", type=int, default=5, help="Number of frames for full TDDFT mode.")
    p.add_argument("--stride", type=int, default=1, help="Trajectory frame stride before sampling.")
    p.add_argument("--random", action="store_true", help="Random frame selection (else evenly spaced).")
    p.add_argument("--seed", type=int, default=20260618, help="RNG seed for random frame selection.")
    p.add_argument("--epsilon", type=float, default=1.77, help="Dielectric screening (optical limit).")
    p.add_argument("--thresh", type=float, default=1e-7, help="Transition-density grid threshold.")
    p.add_argument("--grid-stride", type=int, default=1, help="Density grid subsampling stride.")
    p.add_argument("--backend", default="auto", choices=["auto", "gpu", "opencl"])
    p.add_argument("--gpu-chunk", type=int, default=10000)
    p.add_argument("--opencl-chunk", type=int, default=20000)
    p.add_argument("--opencl-platform", type=int, default=None)
    p.add_argument("--opencl-device", type=int, default=None)
    p.add_argument("--out", type=Path, default=Path("coupling_sampling_out"))
    p.add_argument("--pipeline-args", default="", help="Extra args forwarded to the pipeline in full mode.")
    p.add_argument("--python", default=sys.executable, help="Python used for full-mode subprocesses.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    if not args.traj.exists():
        print(f"[!] Trajectory not found: {args.traj}")
        sys.exit(1)

    print(f"[*] Reading trajectory {args.traj} ...")
    frames = split_pdb_models(args.traj)
    if not frames:
        print("[!] No frames parsed from trajectory.")
        sys.exit(1)
    print(f"    - {len(frames)} frames in trajectory.")

    n_pick = args.subset if args.mode == "full" else args.n_frames
    idxs = select_frame_indices(len(frames), n_pick, args.stride, args.seed, args.random)
    print(f"    - Sampling {len(idxs)} frames (mode={args.mode}).")

    backend_kwargs = dict(
        backend=args.backend,
        gpu_chunk=args.gpu_chunk,
        opencl_chunk=args.opencl_chunk,
        opencl_platform=args.opencl_platform,
        opencl_device=args.opencl_device,
    )

    rows = []

    if args.mode == "rigid":
        # Resolve the (fixed) transition density + target dipole once.
        if args.density is not None:
            density_file = args.density
            mu_target = None
            print(f"[*] Using explicit density {density_file} (no oscillator renormalisation).")
        else:
            workdir, candidates, tried = autodetect_workdir_and_candidates(args.workdir)
            if not candidates:
                print(f"[!] No excited-state candidates found. Tried: {tried}")
                sys.exit(1)
            print_excited_state_table(candidates)
            state, density_file, _, mode = select_target_state_and_density(
                candidates, workdir, args.density_mode, requested_root=args.root)
            mu_target = oscillator_to_dipole_au(state["ev"], state["osc"])
            print(f"    - Root {state['root']} ({mode}): f={state['osc']:.4f}, "
                  f"|mu|_target={mu_target:.4f} a.u., density={density_file}")

        print("[*] Loading fixed transition density ...")
        pts_opt, q_opt = read_dx(density_file, threshold=args.thresh, stride=args.grid_stride)
        if pts_opt.size == 0:
            print("[!] Transition density empty after thresholding.")
            sys.exit(1)

        # Renormalise to the oscillator-strength dipole, exactly as Stage 3 does.
        if mu_target is not None and np.isfinite(mu_target):
            local_origin = np.mean(pts_opt, axis=0)
            dip_vec = np.dot(q_opt, pts_opt - local_origin)
            dip_mag = np.linalg.norm(dip_vec) / BOHR_TO_ANGSTROM
            if dip_mag > 1e-6:
                q_opt = q_opt * (mu_target / dip_mag)
                print(f"    - Renormalised grid dipole {dip_mag:.4f} -> {mu_target:.4f} a.u.")

        with tempfile.TemporaryDirectory() as td:
            frame_pdb = Path(td) / "frame_dimer.pdb"
            for n, fi in enumerate(idxs, 1):
                flines = frames[fi]
                if not frame_has_chains(flines):
                    print(f"    - frame {fi}: missing chain A/B, skipped.")
                    continue
                write_frame_pdb(flines, frame_pdb)
                res = coupling_for_frame(frame_pdb, args.monomer, pts_opt, q_opt,
                                         args.epsilon, backend_kwargs)
                if res is None or "error" in res:
                    msg = res.get("error") if res else "unknown"
                    print(f"[!] frame {fi}: {msg}")
                    if "PyMOL" in str(msg):
                        sys.exit(1)  # no point continuing without PyMOL
                    continue
                res["frame"] = fi
                rows.append(res)
                print(f"    [{n}/{len(idxs)}] frame {fi}: J={res['J_cm']:.2f} cm^-1 "
                      f"(PDA {res['J_pda_cm']:.2f}, sep {res['separation_A']:.2f} A)")

    else:  # full mode
        pipeline_args = args.pipeline_args.split() if args.pipeline_args else []
        with tempfile.TemporaryDirectory() as td:
            for n, fi in enumerate(idxs, 1):
                flines = frames[fi]
                if not frame_has_chains(flines):
                    print(f"    - frame {fi}: missing chain A/B, skipped.")
                    continue
                frame_dir = args.out / f"full_frame_{fi:04d}"
                frame_dir.mkdir(parents=True, exist_ok=True)
                frame_pdb = (frame_dir / "frame_dimer.pdb").resolve()
                write_frame_pdb(flines, frame_pdb)
                J = coupling_for_frame_full(frame_pdb, frame_dir, pipeline_args, args.python)
                if J is None:
                    continue
                rows.append({"frame": fi, "J_cm": J})
                print(f"    [{n}/{len(idxs)}] frame {fi}: J={J:.2f} cm^-1 (full pipeline)")

    if not rows:
        print("[!] No successful frames; nothing to summarise.")
        sys.exit(1)

    j_values = [r["J_cm"] for r in rows]
    stats = summarize(j_values)
    stats["mode"] = args.mode
    stats["epsilon"] = args.epsilon
    stats["samples"] = j_values

    csv_path = args.out / "coupling_samples.csv"
    json_path = args.out / "coupling_distribution.json"
    fig_path = args.out / "Fig_Coupling_Histogram.pdf"
    write_csv(rows, csv_path)
    with open(json_path, "w") as f:
        json.dump(stats, f, indent=2)
    try:
        plot_histogram(j_values, stats, fig_path)
    except Exception as exc:  # plotting is non-fatal
        print(f"[!] Histogram plotting failed: {exc}")
        fig_path = None

    print("\n" + "=" * 52)
    print(f"  J = {stats['mean']:.2f} +/- {stats['std']:.2f} cm^-1  (n={stats['n']})")
    print(f"  range [{stats.get('min', float('nan')):.2f}, {stats.get('max', float('nan')):.2f}] cm^-1")
    print(f"  Davydov splitting 2|J| = {2*abs(stats['mean']):.2f} cm^-1")
    print("=" * 52)
    print(f"  samples : {csv_path}")
    print(f"  summary : {json_path}")
    if fig_path:
        print(f"  figure  : {fig_path}")


if __name__ == "__main__":
    main()

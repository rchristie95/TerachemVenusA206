#!/usr/bin/env python3
"""
coupling_core.py

Pure (numpy / numba / pyopencl / pymol) building blocks for the Transition
Density Coupling (TDC) analysis, factored out of terachem_full_pipeline.py so
that lightweight analysis scripts can reuse them WITHOUT importing OpenMM /
PDBFixer (which terachem_full_pipeline.py loads at module import time).

Nothing in this module imports openmm. The numba/pyopencl/pymol dependencies are
all imported lazily or behind try/except, so importing this module is cheap and
safe on machines without a GPU.

terachem_full_pipeline.py imports every public name defined here, so the
Davydov-coupling pipeline behaviour is unchanged by this extraction.
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np

# Suppress Numba performance warnings
try:
    from numba.core.errors import NumbaPerformanceWarning
    warnings.simplefilter('ignore', category=NumbaPerformanceWarning)
except ImportError:
    pass

try:
    from numba import cuda
    NUMBA_CUDA_AVAILABLE = True
    NUMBA_CUDA_IMPORT_ERROR = None
except Exception as exc:
    cuda = None
    NUMBA_CUDA_AVAILABLE = False
    NUMBA_CUDA_IMPORT_ERROR = f"numba.cuda import failed: {exc}"

try:
    import pyopencl as cl
    OPENCL_AVAILABLE = True
    OPENCL_IMPORT_ERROR = None
except Exception as exc:
    cl = None
    OPENCL_AVAILABLE = False
    OPENCL_IMPORT_ERROR = f"pyopencl import failed: {exc}"

# Constants
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_CM = 219474.63


def oscillator_to_dipole_au(ev, osc):
    """
    Estimate transition-dipole magnitude (a.u.) from transition energy (eV)
    and oscillator strength.
    """
    try:
        ev = float(ev)
        osc = float(osc)
    except (TypeError, ValueError):
        return np.nan
    if ev <= 0.0 or osc < 0.0:
        return np.nan
    return np.sqrt(1.5 * osc / (ev / HARTREE_TO_EV))


def print_excited_state_table(states, indent="    - "):
    """
    Print all parsed excited-state roots with wavelength and dipole magnitude.
    """
    if not states:
        print(f"{indent}No excited-state roots parsed.")
        return
    print(f"{indent}Excited-state roots (lambda in nm, |mu| in a.u.):")
    for state in sorted(states, key=lambda item: item.get("root", 0)):
        root = state.get("root", "?")
        nm = float(state.get("nm", np.nan))
        osc = float(state.get("osc", np.nan))
        mu = oscillator_to_dipole_au(state.get("ev", np.nan), osc)
        nm_text = f"{nm:.1f}" if np.isfinite(nm) else "nan"
        mu_text = f"{mu:.4f}" if np.isfinite(mu) else "nan"
        osc_text = f"{osc:.4f}" if np.isfinite(osc) else "nan"
        print(f"      Root {root}: lambda={nm_text} nm, |mu|={mu_text} a.u., f={osc_text}")


OPENCL_KERNEL_FP64 = r"""
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
__kernel void coulomb_row_sum_f64(
    __global const double *pts1,
    __global const double *q1,
    __global const double *pts2,
    __global const double *q2,
    const int n1,
    const int n2,
    __global double *out
){
    int i = get_global_id(0);
    if (i >= n1) return;

    double xi = pts1[3*i + 0];
    double yi = pts1[3*i + 1];
    double zi = pts1[3*i + 2];
    double qi = q1[i];
    double acc = 0.0;

    for (int j = 0; j < n2; ++j){
        double dx = xi - pts2[3*j + 0];
        double dy = yi - pts2[3*j + 1];
        double dz = zi - pts2[3*j + 2];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.1) r = 0.1;
        acc += (qi * q2[j]) / r;
    }
    out[i] = acc;
}
"""

OPENCL_KERNEL_FP32 = r"""
__kernel void coulomb_row_sum_f32(
    __global const float *pts1,
    __global const float *q1,
    __global const float *pts2,
    __global const float *q2,
    const int n1,
    const int n2,
    __global float *out
){
    int i = get_global_id(0);
    if (i >= n1) return;

    float xi = pts1[3*i + 0];
    float yi = pts1[3*i + 1];
    float zi = pts1[3*i + 2];
    float qi = q1[i];
    float acc = 0.0f;

    for (int j = 0; j < n2; ++j){
        float dx = xi - pts2[3*j + 0];
        float dy = yi - pts2[3*j + 1];
        float dz = zi - pts2[3*j + 2];
        float r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.1f) r = 0.1f;
        acc += (qi * q2[j]) / r;
    }
    out[i] = acc;
}
"""

def parse_pdb_residue_atoms(pdb_path, res_names=None, res_ids=None, chain_id=None, heavy_only=False):
    """Extracts atoms for specific residues from PDB."""
    atoms = []
    if not pdb_path.exists(): return atoms
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                if chain_id and chain != chain_id:
                    continue
                resName = line[17:20].strip()
                resSeq = line[22:26].strip()
                match = False
                if res_names and resName in res_names: match = True
                if res_ids and resSeq in res_ids: match = True

                if match:
                    if resName in ["HOH", "WAT"]: continue
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    elem = line[76:78].strip()
                    if not elem: elem = line[12:16].strip()[0]
                    if heavy_only and elem.upper() == "H":
                        continue
                    atoms.append({'elem': elem, 'xyz': np.array([x,y,z]), 'resSeq': resSeq})
    return atoms

def get_super_matrices_with_pymol(monomer_pdb, dimer_pdb):
    """
    Reproduce the exact transform workflow from visualise_dimer_big.pml:
      super siteA and chain A, dimer_ref and chain A
      super siteB and chain A, dimer_ref and chain B
    """
    try:
        from pymol import cmd
    except Exception as exc:
        return None, None, None, None, f"PyMOL unavailable: {exc}"

    cmd.reinitialize()
    cmd.load(str(dimer_pdb), "dimer_ref")
    cmd.load(str(monomer_pdb), "siteA")
    cmd.load(str(monomer_pdb), "siteB")
    aln_a = cmd.super("siteA and chain A", "dimer_ref and chain A", quiet=1)
    aln_b = cmd.super("siteB and chain A", "dimer_ref and chain B", quiet=1)
    mat_a = np.array(cmd.get_object_matrix("siteA"), dtype=float).reshape(4, 4)
    mat_b = np.array(cmd.get_object_matrix("siteB"), dtype=float).reshape(4, 4)
    return mat_a, mat_b, aln_a, aln_b, None

def apply_pymol_matrix(points, matrix4):
    """Apply a PyMOL object matrix to an (N,3) point array."""
    rot = matrix4[:3, :3]
    trans = matrix4[:3, 3]
    return np.dot(points, rot.T) + trans

def transition_dipole_au(points_angstrom, charges, origin_angstrom=None):
    """
    Compute transition dipole in atomic units from point charges.

    If origin_angstrom is provided, positions are shifted by that origin first.
    This keeps dipole diagnostics robust when residual net charge is non-zero.
    """
    pts_bohr = points_angstrom * ANGSTROM_TO_BOHR
    if origin_angstrom is not None:
        pts_bohr = pts_bohr - (origin_angstrom * ANGSTROM_TO_BOHR)
    return np.sum(charges[:, None] * pts_bohr, axis=0)

def print_transform_matrix(name, matrix4):
    print(f"    - {name} transform matrix:")
    for row in matrix4:
        print("      " + " ".join(f"{v: .8f}" for v in row))

def read_dx(filename, threshold=1e-7, stride=1):
    with open(filename, 'r') as f:
        counts = None
        origin = None
        deltas = []
        # Read header
        for line in f:
            if line.startswith('object 1'):
                counts = np.array([int(x) for x in line.split()[5:8]], dtype=int)
            elif line.startswith('origin'):
                origin = np.array([float(x) for x in line.split()[1:4]], dtype=float)
            elif line.startswith('delta'):
                deltas.append([float(x) for x in line.split()[1:4]])
            elif line.startswith('object 3'):
                break
        raw = f.read()
    if counts is None or origin is None or len(deltas) != 3:
        raise ValueError(f"DX header parse failed for {filename}")

    # Trim any trailing attribute blocks
    if 'attribute' in raw:
        raw = raw.split('attribute')[0]

    # Parse numeric payload
    vals = np.fromstring(raw, sep=' ')
    nx, ny, nz = counts
    expected = nx * ny * nz
    if vals.size < expected:
        raise ValueError(
            f"DX payload too small in {filename}: got {vals.size} values, expected {expected}"
        )
    if vals.size > expected:
        # Some DX files include trailing numeric blocks; keep only grid payload.
        vals = vals[:expected]

    vals = vals.reshape((nx, ny, nz))

    deltas = np.array(deltas, dtype=float)
    d_vol = abs(np.linalg.det(deltas))  # voxel volume

    # 3D striding
    if stride is None or stride < 1:
        stride = 1

    # Subsample the grid
    vals_s = vals[::stride, ::stride, ::stride]
    scale = float(stride**3)

    # Thresholding: keep signed density
    mask = (vals_s > threshold) | (vals_s < -threshold)
    if not np.any(mask):
        return np.array([]), np.array([])

    # Get indices of points above threshold in the subsampled grid
    ix_s, iy_s, iz_s = np.where(mask)

    # Map subsampled indices back to full-grid coordinates
    # DX convention: origin + i*dx + j*dy + k*dz
    pts = origin + np.outer(ix_s * stride, deltas[0]) + \
                   np.outer(iy_s * stride, deltas[1]) + \
                   np.outer(iz_s * stride, deltas[2])

    vals_kept = vals_s[ix_s, iy_s, iz_s]
    # The values in TeraChem DX files are often already pre-multiplied by d_vol.
    # We use the raw values as discrete charges to avoid redundant scaling.
    charges = vals_kept * scale
    return pts, charges

if NUMBA_CUDA_AVAILABLE:
    @cuda.jit
    def coupling_cuda_kernel(pts1, q1, pts2, q2, result):
        i = cuda.grid(1)
        if i >= pts1.shape[0]:
            return

        xi = pts1[i, 0]
        yi = pts1[i, 1]
        zi = pts1[i, 2]
        qi = q1[i]

        acc = 0.0
        n2 = pts2.shape[0]
        for j in range(n2):
            dx = xi - pts2[j, 0]
            dy = yi - pts2[j, 1]
            dz = zi - pts2[j, 2]
            r = (dx * dx + dy * dy + dz * dz) ** 0.5
            if r < 0.1:
                r = 0.1
            acc += (qi * q2[j]) / r

        cuda.atomic.add(result, 0, acc)

def _is_cuda_ready():
    if not NUMBA_CUDA_AVAILABLE:
        return False
    try:
        return cuda.is_available()
    except Exception:
        return False

def _cuda_unavailable_reason():
    if not NUMBA_CUDA_AVAILABLE:
        return NUMBA_CUDA_IMPORT_ERROR
    try:
        if cuda.is_available():
            return None
    except Exception as exc:
        return f"cuda.is_available() error: {exc}"

    # cuda.is_available() is False: report likely driver/runtime cause.
    try:
        _ = len(cuda.gpus)
    except Exception as exc:
        return f"no CUDA device/driver visible ({exc})"
    return "no CUDA-capable device visible to this Python environment"

def calculate_coupling_gpu(pts1, q1, pts2, q2, gpu_chunk=10000):
    if not _is_cuda_ready():
        reason = _cuda_unavailable_reason() or "unknown CUDA initialization issue"
        raise RuntimeError(f"CUDA backend is not available: {reason}")

    pts1 = np.ascontiguousarray(pts1, dtype=np.float64)
    q1 = np.ascontiguousarray(q1, dtype=np.float64)
    pts2 = np.ascontiguousarray(pts2, dtype=np.float64)
    q2 = np.ascontiguousarray(q2, dtype=np.float64)

    if len(q1) == 0 or len(q2) == 0:
        return 0.0

    device = cuda.get_current_device()
    dev_name = getattr(device, "name", "CUDA device")
    if isinstance(dev_name, bytes):
        dev_name = dev_name.decode(errors="ignore")

    threads_per_block = 256
    print(f"    - Calculating J on GPU ({dev_name})...")

    start_t = time.time()
    j_sum = 0.0
    n_chunks = (len(q1) + gpu_chunk - 1) // gpu_chunk

    d_pts2 = cuda.to_device(pts2)
    d_q2 = cuda.to_device(q2)

    for chunk_idx, i_start in enumerate(range(0, len(q1), gpu_chunk), start=1):
        i_end = min(i_start + gpu_chunk, len(q1))
        pts1_chunk = pts1[i_start:i_end]
        q1_chunk = q1[i_start:i_end]

        d_pts1 = cuda.to_device(pts1_chunk)
        d_q1 = cuda.to_device(q1_chunk)
        d_result = cuda.to_device(np.array([0.0], dtype=np.float64))

        blocks = (len(pts1_chunk) + threads_per_block - 1) // threads_per_block
        coupling_cuda_kernel[blocks, threads_per_block](d_pts1, d_q1, d_pts2, d_q2, d_result)
        cuda.synchronize()
        j_sum += float(d_result.copy_to_host()[0])

        if n_chunks > 5 and (chunk_idx % max(1, n_chunks // 5) == 0 or chunk_idx == n_chunks):
            print(f"      GPU progress: {chunk_idx}/{n_chunks} chunks")

    print(f"    - Calculation took {time.time() - start_t:.2f} seconds.")
    return j_sum * ANGSTROM_TO_BOHR

def _is_opencl_ready():
    if not OPENCL_AVAILABLE:
        return False
    try:
        return len(cl.get_platforms()) > 0
    except Exception:
        return False

def _opencl_unavailable_reason():
    if not OPENCL_AVAILABLE:
        return OPENCL_IMPORT_ERROR
    try:
        platforms = cl.get_platforms()
    except Exception as exc:
        return f"OpenCL platform query failed ({exc})"
    if not platforms:
        return "no OpenCL platforms found"
    for platform in platforms:
        try:
            if platform.get_devices():
                return None
        except Exception:
            continue
    return "no OpenCL devices found"

def _create_opencl_context(opencl_platform=None, opencl_device=None):
    platforms = cl.get_platforms()
    if not platforms:
        raise RuntimeError("no OpenCL platforms found")

    if opencl_platform is not None:
        if opencl_platform < 0 or opencl_platform >= len(platforms):
            raise ValueError(f"opencl-platform index {opencl_platform} out of range [0, {len(platforms)-1}]")
        platform = platforms[opencl_platform]
    else:
        # Prefer a platform with at least one GPU device.
        platform = None
        for p in platforms:
            try:
                if p.get_devices(device_type=cl.device_type.GPU):
                    platform = p
                    break
            except Exception:
                continue
        if platform is None:
            platform = platforms[0]

    devices = platform.get_devices()
    if not devices:
        raise RuntimeError(f"platform '{platform.name}' has no devices")

    if opencl_device is not None:
        if opencl_device < 0 or opencl_device >= len(devices):
            raise ValueError(f"opencl-device index {opencl_device} out of range [0, {len(devices)-1}]")
        device = devices[opencl_device]
    else:
        gpu_devices = [d for d in devices if d.type & cl.device_type.GPU]
        device = gpu_devices[0] if gpu_devices else devices[0]

    ctx = cl.Context([device])
    queue = cl.CommandQueue(ctx)
    return platform, device, ctx, queue

def _build_opencl_program(ctx, device):
    extensions = (getattr(device, "extensions", "") or "").lower()
    fp64_supported = ("cl_khr_fp64" in extensions) or ("cl_amd_fp64" in extensions)

    if fp64_supported:
        try:
            program = cl.Program(ctx, OPENCL_KERNEL_FP64).build()
            return program, np.float64, "f64"
        except Exception as exc:
            print(f"    - OpenCL fp64 build failed ({exc}); trying fp32 kernel.")

    try:
        program = cl.Program(ctx, OPENCL_KERNEL_FP32).build()
        return program, np.float32, "f32"
    except Exception as exc:
        raise RuntimeError(f"OpenCL kernel build failed: {exc}")

def calculate_coupling_opencl(
    pts1, q1, pts2, q2, opencl_chunk=20000, opencl_platform=None, opencl_device=None
):
    if not _is_opencl_ready():
        reason = _opencl_unavailable_reason() or "unknown OpenCL initialization issue"
        raise RuntimeError(f"OpenCL backend is not available: {reason}")

    pts1 = np.ascontiguousarray(pts1, dtype=np.float64)
    q1 = np.ascontiguousarray(q1, dtype=np.float64)
    pts2 = np.ascontiguousarray(pts2, dtype=np.float64)
    q2 = np.ascontiguousarray(q2, dtype=np.float64)
    if len(q1) == 0 or len(q2) == 0:
        return 0.0

    platform, device, ctx, queue = _create_opencl_context(opencl_platform, opencl_device)
    program, kernel_dtype, precision_label = _build_opencl_program(ctx, device)
    kernel = getattr(program, f"coulomb_row_sum_{precision_label}")
    mf = cl.mem_flags

    pts2_dev_host = np.ascontiguousarray(pts2, dtype=kernel_dtype).reshape(-1)
    q2_dev_host = np.ascontiguousarray(q2, dtype=kernel_dtype)
    d_pts2 = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pts2_dev_host)
    d_q2 = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q2_dev_host)

    max_wg = int(getattr(device, "max_work_group_size", 256))
    local_size = min(256, max_wg) if max_wg > 0 else 64

    print(f"    - Calculating J on OpenCL ({platform.name} | {device.name})...")

    start_t = time.time()
    j_sum = 0.0
    n2 = np.int32(len(q2_dev_host))
    n_chunks = (len(q1) + opencl_chunk - 1) // opencl_chunk

    for chunk_idx, i_start in enumerate(range(0, len(q1), opencl_chunk), start=1):
        i_end = min(i_start + opencl_chunk, len(q1))
        n1_chunk = i_end - i_start
        pts1_chunk = np.ascontiguousarray(pts1[i_start:i_end], dtype=kernel_dtype).reshape(-1)
        q1_chunk = np.ascontiguousarray(q1[i_start:i_end], dtype=kernel_dtype)
        out_host = np.empty(n1_chunk, dtype=kernel_dtype)

        d_pts1 = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pts1_chunk)
        d_q1 = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q1_chunk)
        d_out = cl.Buffer(ctx, mf.WRITE_ONLY, size=out_host.nbytes)

        global_size = ((n1_chunk + local_size - 1) // local_size) * local_size
        kernel(
            queue,
            (global_size,),
            (local_size,),
            d_pts1,
            d_q1,
            d_pts2,
            d_q2,
            np.int32(n1_chunk),
            n2,
            d_out,
        )
        cl.enqueue_copy(queue, out_host, d_out)
        queue.finish()
        j_sum += float(np.sum(out_host, dtype=np.float64))

        if n_chunks > 5 and (chunk_idx % max(1, n_chunks // 5) == 0 or chunk_idx == n_chunks):
            print(f"      OpenCL progress: {chunk_idx}/{n_chunks} chunks")

    print(f"    - Calculation took {time.time() - start_t:.2f} seconds.")
    return j_sum * ANGSTROM_TO_BOHR

def calculate_coupling(
    pts1,
    q1,
    pts2,
    q2,
    backend="auto",
    gpu_chunk=10000,
    opencl_chunk=20000,
    opencl_platform=None,
    opencl_device=None,
):
    backend = backend.lower()
    if backend not in {"auto", "gpu", "opencl"}:
        raise ValueError("Unknown backend '{}'. Expected one of: auto, gpu, opencl.".format(backend))

    # Auto: try CUDA, but fallback safely to OpenCL if it fails (e.g. PTX error)
    if backend == "auto":
        if _is_cuda_ready():
            try:
                return calculate_coupling_gpu(pts1, q1, pts2, q2, gpu_chunk=gpu_chunk)
            except Exception as exc:
                print(f"    - CUDA backend failed ({exc}); falling back to OpenCL.")
        else:
            reason = _cuda_unavailable_reason() or "unknown reason"
            print(f"    - CUDA backend unavailable ({reason}); trying OpenCL.")

        return calculate_coupling_opencl(
            pts1, q1, pts2, q2,
            opencl_chunk=opencl_chunk,
            opencl_platform=opencl_platform,
            opencl_device=opencl_device
        )

    # Explicit requests
    if backend == "gpu":
        try:
            return calculate_coupling_gpu(pts1, q1, pts2, q2, gpu_chunk=gpu_chunk)
        except Exception as exc:
            # Even if explicitly requested, if it fails due to the known PTX version issue,
            # we offer to fallback or just fail. Given the "fix" mandate, let's fallback with a warning.
            print(f"    - CUDA backend failed ({exc}); falling back to OpenCL.")
            return calculate_coupling_opencl(
                pts1, q1, pts2, q2,
                opencl_chunk=opencl_chunk,
                opencl_platform=opencl_platform,
                opencl_device=opencl_device
            )

    if backend == "opencl":
        return calculate_coupling_opencl(
            pts1,
            q1,
            pts2,
            q2,
            opencl_chunk=opencl_chunk,
            opencl_platform=opencl_platform,
            opencl_device=opencl_device,
        )

def parse_excited_state_candidates(energy_file):
    """Parse excited-state summary table from TeraChem energy.out."""
    candidates = []
    if not energy_file.exists():
        return candidates

    in_table = False
    with open(energy_file, 'r') as f:
        for line in f:
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
            if len(parts) < 4 or not parts[0].isdigit():
                continue
            try:
                r = int(parts[0])
                ev = float(parts[2])
                osc = float(parts[3])
            except ValueError:
                continue
            nm = 1239.84 / ev if ev != 0 else np.inf
            candidates.append({'root': r, 'ev': ev, 'osc': osc, 'nm': nm})
    return candidates

def autodetect_workdir_and_candidates(preferred_workdir):
    """
    Return (workdir, candidates) by checking preferred path first, then
    tc_tddft_old_current* directories in cwd.
    """
    search_dirs = []
    seen = set()

    def add_dir(path_obj):
        path_obj = Path(path_obj)
        key = str(path_obj.resolve()) if path_obj.exists() else str(path_obj)
        if key in seen:
            return
        seen.add(key)
        search_dirs.append(path_obj)

    add_dir(preferred_workdir)
    for p in sorted(Path.cwd().glob("tc_tddft_old_current*")):
        if p.is_dir():
            add_dir(p)

    tried = []
    fallback = None

    for idx, wd in enumerate(search_dirs):
        energy_file = wd / "energy.out"
        cand = parse_excited_state_candidates(energy_file)
        tried.append((wd, len(cand)))
        if not cand:
            continue

        if idx == 0:
            return wd, cand, tried
        if fallback is None:
            fallback = (wd, cand)

    if fallback is not None:
        return fallback[0], fallback[1], tried

    return preferred_workdir, [], tried

def find_transition_density_file(workdir, root, mode="signed"):
    """
    Resolve transition-density DX file path.
    """
    signed_candidates = [
        workdir / f"transdens_{root}.dx",
        workdir / "scr_plot" / f"transdens_{root}.dx",
    ]
    legacy_abs_candidates = [
        workdir / f"abs_transdens_{root}.dx",
    ]

    def first_existing(paths):
        for p in paths:
            if p.exists():
                return p
        return None

    if mode in ("signed", "auto"):
        p = first_existing(signed_candidates)
        if p is not None:
            return p, "signed"
        if mode == "signed":
            tried = ", ".join(str(x) for x in signed_candidates)
            raise FileNotFoundError(
                f"Signed transition density not found for root {root}. Tried: {tried}. "
                "Cannot reconstruct signed density from absolute-valued data."
            )

    p = first_existing(legacy_abs_candidates)
    if p is not None:
        label = "abs (requested)" if mode == "abs" else "legacy abs filename (auto fallback)"
        return p, label

    tried_paths = signed_candidates + legacy_abs_candidates
    tried = ", ".join(str(x) for x in tried_paths)
    raise FileNotFoundError(f"Transition density file not found for root {root}. Tried: {tried}")

def select_target_state_and_density(candidates, workdir, density_mode, requested_root=None):
    """
    Choose excited state and matching density file.
    """
    if requested_root is not None:
        state = next((c for c in candidates if c['root'] == requested_root), None)
        if state is None:
            roots = ", ".join(str(c['root']) for c in candidates)
            raise ValueError(f"Root {requested_root} not found in {workdir}. Available roots: {roots}")
        density_file, density_kind = find_transition_density_file(workdir, state['root'], mode=density_mode)
        return state, density_file, density_kind, "requested"

    visible = [c for c in candidates if 450 <= c['nm'] <= 600]
    visible_roots = {c['root'] for c in visible}
    ranked_visible = sorted(visible, key=lambda x: x['osc'], reverse=True)
    ranked_rest = sorted(
        [c for c in candidates if c['root'] not in visible_roots],
        key=lambda x: x['osc'],
        reverse=True,
    )

    for state in ranked_visible:
        try:
            density_file, density_kind = find_transition_density_file(workdir, state['root'], mode=density_mode)
            return state, density_file, density_kind, "auto-visible"
        except FileNotFoundError:
            continue

    for state in ranked_rest:
        try:
            density_file, density_kind = find_transition_density_file(workdir, state['root'], mode=density_mode)
            return state, density_file, density_kind, "auto-global"
        except FileNotFoundError:
            continue

    raise FileNotFoundError(
        "No candidate root with an available transition-density file was found "
        f"in {workdir} (mode={density_mode})."
    )

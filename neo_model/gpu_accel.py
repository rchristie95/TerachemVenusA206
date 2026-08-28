#!/usr/bin/env python3
"""
GPU acceleration (gpu4pyscf / RTX 4080) for the chromophore NEO pipeline.

Three levers, in order of how cleanly GPU applies:
  1. gpu_tddft()      : electronic TDDFT on GPU            -> the excitation step
  2. gpu_optimize()   : classical electronic SCF+grad+opt  -> the `cl` geometry opt
  3. gpu_neo_hook()   : routes the NEO *electronic* full veff (J+K+XC) to the GPU,
                        leaving the (tiny) quantum-proton SCF + e-p coupling on the
                        CPU. The electronic structure is the dominant cost of
                        NEO-DFT, so this is the real "faster NEO".

We replace the electronic ComponentSCF.get_veff with a gpu4pyscf RKS veff engine
(grids cached), re-tag ecoul/exc, and re-add NEO's cached inter-component
potential (_vint) so the multicomponent physics is identical. Always verify the
GPU NEO energy against CPU NEO before trusting a production run.
"""
import numpy as np


def _ensure_pyscf_version():
    """The editable fork installs pyscf as a namespace package (no __version__);
    gpu4pyscf._patch_pyscf needs it. Set it before any gpu4pyscf import."""
    import pyscf
    if not getattr(pyscf, "__version__", None):
        pyscf.__version__ = "2.13.1"


_ensure_pyscf_version()


_MANAGED = False


def enable_managed_memory(threshold=100_000_000):
    """Override gpu4pyscf's allocator so LARGE (>100 MB) allocations use CUDA
    managed (unified) memory, which oversubscribes VRAM and spills to host RAM;
    small allocations keep using the fast device pool. Needed for full TDDFT of
    this anion on a 16 GB card (the non-Hermitian eigensolver peaks past ~15 GB).

    gpu4pyscf installs its own conditional allocator at import (__init__.py),
    routing large allocations to direct CUDA malloc -> physical VRAM -> OOM. So
    this MUST be called *after* gpu4pyscf is imported to take effect. Idempotent.
    """
    global _MANAGED
    if _MANAGED:
        return
    import cupy
    managed = cupy.cuda.MemoryPool(cupy.cuda.malloc_managed)
    pool_malloc = cupy.get_default_memory_pool().malloc

    def malloc(size):
        return managed.malloc(size) if size >= threshold else pool_malloc(size)

    cupy.cuda.set_allocator(malloc)
    _MANAGED = True


def to_cpu(x):
    """cupy -> numpy (no-op for numpy)."""
    if x is None:
        return None
    get = getattr(x, "get", None)
    return get() if callable(get) else np.asarray(x)


# ---------------------------------------------------------------- electronic GPU
def gpu_rks(mol, xc="b3lyp"):
    from gpu4pyscf.dft import rks
    mf = rks.RKS(mol, xc=xc)
    return mf


def gpu_scf(mol, xc="b3lyp"):
    mf = gpu_rks(mol, xc)
    e = float(to_cpu(mf.kernel()))
    return mf, e


def _osc_from_xy(mf, td):
    """Oscillator strengths from TDDFT (X,Y) + dipole integrals (gpu4pyscf RKS
    TDDFT has no oscillator_strength()). Absolute scale uses the standard
    f = (2/3) E |mu|^2 with the spin-restricted sqrt(2) convention; the bright
    state is robustly the argmax regardless of small prefactor ambiguities."""
    mo = to_cpu(mf.mo_coeff); occ = to_cpu(mf.mo_occ)
    orbo = mo[:, occ > 0]; orbv = mo[:, occ == 0]
    ao_dip = mf.mol.intor("int1e_r")                      # (3, nao, nao)
    dip_ov = np.einsum("xpq,pi,qa->xia", ao_dip, orbo, orbv)
    e = np.asarray(to_cpu(td.e))
    fs = []
    for xy in td.xy:
        x = to_cpu(xy[0]); y = xy[1]
        xpy = x + (to_cpu(y) if np.ndim(y) > 0 else 0.0)
        tdip = np.einsum("xia,ia->x", dip_ov, xpy) * 2.0   # restricted; xy stored *sqrt(.5)
        fs.append((2.0 / 3.0) * e[len(fs)] * float(tdip @ tdip))
    return np.array(fs)


def gpu_tddft(mol, xc="b3lyp", nstates=8, tda=False):
    """Bright pi->pi* on GPU. Returns (energies_eV, osc, bright_index, mf).
    Frees the cupy pool first and falls back to TDA on out-of-memory."""
    import cupy
    EV = 27.211386
    mf = gpu_rks(mol, xc)          # imports gpu4pyscf (which installs its own allocator)...
    enable_managed_memory()        # ...so override AFTER: large allocs -> managed (spill to host)
    mf.kernel()
    from gpu4pyscf.tdscf import rks as gtd
    def _is_oom(ex):
        return "memory" in str(ex).lower() or "OutOfMemory" in type(ex).__name__
    try:
        td = (gtd.TDA(mf) if tda else gtd.TDDFT(mf)); td.nstates = nstates; td.kernel()
    except Exception as ex:                      # full-TDDFT dense eigensolver can exceed VRAM
        if not _is_oom(ex):
            raise
        cupy.get_default_memory_pool().free_all_blocks()
        td = gtd.TDA(mf); td.nstates = nstates; td.kernel()   # TDA ~ 1/4 the memory
        print("  [gpu_tddft] fell back to TDA (out of GPU memory for full TDDFT)")
    e = np.asarray(to_cpu(td.e)) * EV
    osc = _osc_from_xy(mf, td)
    return e, osc, int(np.argmax(osc)), mf


def gpu_optimize(mol, xc="b3lyp", maxsteps=80, callback=None):
    """Classical geometry optimization with GPU SCF+gradients."""
    from pyscf.geomopt.geometric_solver import optimize
    mf = gpu_rks(mol, xc)
    return optimize(mf, maxsteps=maxsteps, callback=callback)


# ---------------------------------------------------------------- NEO electronic J/K on GPU
def gpu_neo_hook(nmf, xc="b3lyp"):
    """Route the NEO *electronic* component's FULL veff (J + K + XC) to the GPU.

    This is the dominant cost of NEO-DFT; the quantum-proton SCF and the e-p
    Coulomb coupling are tiny and stay on CPU. The electronic veff is computed
    by a persistent gpu4pyscf RKS engine (grids built once), converted back to
    numpy, re-tagged with ecoul/exc, and the cached inter-component potential
    (_vint) is added exactly as NEO's ComponentSCF.get_veff does.

    Returns nmf (modified in place); call nmf.scf() as usual. Always verify the
    resulting energy against CPU NEO before a production run.
    """
    import cupy
    from pyscf import lib
    from gpu4pyscf.dft import rks as grks

    ecomp = nmf.components["e"]
    mol_e = ecomp.mol
    gmf = grks.RKS(mol_e, xc=xc)          # persistent GPU engine (grids cached)

    def gpu_get_veff(mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
        vg = gmf.get_veff(mol_e, cupy.asarray(np.asarray(dm)))   # tagged cupy veff
        veff = lib.tag_array(to_cpu(vg), ecoul=float(vg.ecoul), exc=float(vg.exc),
                             vj=to_cpu(vg.vj),
                             vk=(to_cpu(vg.vk) if vg.vk is not None else None))
        # add inter-component (e-p Coulomb / epc) potential, mirroring NEO
        if ecomp._vint is not None:
            vhf_self = veff.view(np.ndarray).copy()
            vhf_self = lib.tag_array(vhf_self, **veff.__dict__)
            veff = lib.tag_array(veff + ecomp._vint, **veff.__dict__,
                                 vhf_self=vhf_self, vint=ecomp._vint)
        return veff

    ecomp.get_veff = gpu_get_veff
    ecomp._gpu_engine = gmf               # keep reference alive
    nmf._gpu_neo = True
    return nmf


# ---------------------------------------------------------------- self-test / benchmark
def _read_model():
    import json
    idx = json.load(open("idx_full.json"))
    L = open("chromo_his.xyz").read().splitlines(); n = int(L[0].split()[0])
    atomstr = "\n".join(L[2:2+n])
    return atomstr, idx


if __name__ == "__main__":
    import time, json
    from pyscf import gto, dft, neo, tddft
    EV, NM = 27.211386, 1239.841984
    atomstr, idx = _read_model()
    HD1 = idx["HD1"]

    print("=== GPU sanity: cupy + gpu4pyscf ===", flush=True)
    import cupy
    print("  cupy device:", cupy.cuda.runtime.getDeviceProperties(0)["name"].decode())

    mol = gto.M(atom=atomstr, basis="6-31g*", charge=-1, spin=0, verbose=0, max_memory=30000)

    # (1) TDDFT  CPU vs GPU
    print("\n=== TDDFT  CPU vs GPU (b3lyp/6-31g*) ===", flush=True)
    t = time.time(); mfc = dft.RKS(mol); mfc.xc = "b3lyp"; mfc.kernel()
    td = tddft.TDDFT(mfc); td.nstates = 8; td.kernel()
    oc = td.oscillator_strength(); ec = np.asarray(td.e)*EV; ic = int(np.argmax(oc)); tc = time.time()-t
    print(f"  CPU: bright {ec[ic]:.3f} eV / {NM/ec[ic]:.1f} nm  ({tc:.1f}s)")
    t = time.time(); eg, og, ig, _ = gpu_tddft(mol, "b3lyp", 8); tg = time.time()-t
    print(f"  GPU: bright {eg[ig]:.3f} eV / {NM/eg[ig]:.1f} nm  ({tg:.1f}s)   speedup {tc/tg:.1f}x"
          f"   dE={(eg[ig]-ec[ic])*1000:.2f} meV")

    # (3) NEO SCF  CPU vs GPU-hooked
    print("\n=== NEO-DFT SCF  CPU vs GPU-electronic-J/K ===", flush=True)
    def neomol():
        return neo.M(atom=atomstr, basis="6-31g*", nuc_basis="pb4d", charge=-1,
                     spin=0, quantum_nuc=[HD1], verbose=0, max_memory=30000)
    t = time.time(); ncpu = neo.KS(neomol(), xc="b3lyp"); ecpu = ncpu.scf(); tcpu = time.time()-t
    print(f"  CPU NEO: E={ecpu:.6f}  ({tcpu:.1f}s)  conv={ncpu.converged}")
    t = time.time(); ngpu = gpu_neo_hook(neo.KS(neomol(), xc="b3lyp"), "b3lyp")
    egpu = ngpu.scf(); tgpu = time.time()-t
    print(f"  GPU NEO: E={egpu:.6f}  ({tgpu:.1f}s)  conv={ngpu.converged}  speedup {tcpu/tgpu:.1f}x")
    print(f"  >>> energy match: dE = {(egpu-ecpu)*1e6:.2f} micro-Eh  "
          f"({'OK' if abs(egpu-ecpu) < 1e-5 else 'MISMATCH -- do NOT trust GPU NEO'})")

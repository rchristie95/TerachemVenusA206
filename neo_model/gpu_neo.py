#!/usr/bin/env python3
"""
GPU-native quantum chemistry engine for the chromophore pipeline (RTX 4080).

Backends: gpu4pyscf (GPU electronic DFT/TDDFT) + the PySCF NEO fork (quantum
nuclei). NEO-DFT runs with its dominant electronic veff on the GPU via
gpu_accel.gpu_neo_hook; the quantum-proton SCF + e-p coupling stay on CPU.

  python gpu_neo.py dft   --geom g.xyz --charge -1 --basis 6-31g* --xc b3lyp
  python gpu_neo.py tddft --geom g.xyz --charge -1 --nstates 8
  python gpu_neo.py neo   --geom g.xyz --charge -1 --quantum 36 [--epc 17-2] [--verify]

--verify (neo) also runs CPU NEO and asserts the energy matches (<5e-5 Eh).
Run from this directory with the conda env `pyscfneo` and:
  export LD_LIBRARY_PATH=$CONDA_PREFIX/../../lib:$LD_LIBRARY_PATH  (MKL) ; MKL_THREADING_LAYER=GNU
"""
import argparse, sys, time, numpy as np
import gpu_accel as G   # runs the pyscf-version shim on import

EV, NM = 27.211386, 1239.841984


def read_geom(path):
    L = open(path).read().splitlines(); n = int(L[0].split()[0])
    return "\n".join(L[2:2+n]), n


def main():
    p = argparse.ArgumentParser(description="GPU-native DFT / TDDFT / NEO-DFT")
    p.add_argument("mode", choices=["dft", "tddft", "neo"])
    p.add_argument("--geom", required=True, help="xyz file")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--basis", default="6-31g*")
    p.add_argument("--xc", default="b3lyp")
    p.add_argument("--nstates", type=int, default=8)
    p.add_argument("--quantum", default="", help="comma-sep atom indices treated quantum (neo)")
    p.add_argument("--nuc-basis", default="pb4d")
    p.add_argument("--epc", default=None, help="e-p correlation functional, e.g. 17-2 (neo SCF)")
    p.add_argument("--verify", action="store_true", help="also run CPU and check energy (neo)")
    a = p.parse_args()
    atomstr, n = read_geom(a.geom)
    from pyscf import gto
    mol = gto.M(atom=atomstr, basis=a.basis, charge=a.charge, spin=a.spin,
                verbose=0, max_memory=30000)
    print(f"# {a.mode}: {n} atoms  basis={a.basis} xc={a.xc} charge={a.charge} (GPU: RTX 4080)")

    if a.mode == "dft":
        t = time.time(); mf, e = G.gpu_scf(mol, a.xc)
        print(f"GPU DFT  E = {e:.8f} Eh   ({time.time()-t:.1f}s)")

    elif a.mode == "tddft":
        t = time.time(); e, osc, i, _ = G.gpu_tddft(mol, a.xc, a.nstates); dt = time.time()-t
        for k in range(len(e)):
            print(f"  S{k+1}: {e[k]:.3f} eV  {NM/e[k]:6.1f} nm  f={osc[k]:.3f}" +
                  ("   <== bright" if k == i else ""))
        print(f"GPU TDDFT bright pi->pi*: {e[i]:.3f} eV = {NM/e[i]:.1f} nm (f={osc[i]:.3f})  ({dt:.1f}s)")

    elif a.mode == "neo":
        from pyscf import neo
        qn = [int(x) for x in a.quantum.split(",") if x != ""]
        if not qn:
            sys.exit("neo needs --quantum <atom indices>")
        def build():
            return neo.M(atom=atomstr, basis=a.basis, nuc_basis=a.nuc_basis, charge=a.charge,
                         spin=a.spin, quantum_nuc=qn, verbose=0, max_memory=30000)
        nmf = neo.KS(build(), xc=a.xc, epc=a.epc)
        G.gpu_neo_hook(nmf, a.xc)
        t = time.time(); e = nmf.scf(); dt = time.time()-t
        print(f"GPU NEO-DFT  E = {e:.8f} Eh   ({dt:.1f}s)  converged={nmf.converged}  "
              f"quantum_nuc={qn} epc={a.epc}")
        if a.verify:
            tc = time.time(); ncpu = neo.KS(build(), xc=a.xc, epc=a.epc); ec = ncpu.scf(); tcd = time.time()-tc
            ok = abs(e - ec) < 5e-5
            print(f"CPU NEO-DFT  E = {ec:.8f} Eh   ({tcd:.1f}s)")
            print(f">>> dE = {(e-ec)*1e6:+.1f} micro-Eh  speedup {tcd/dt:.2f}x  "
                  f"[{'VERIFIED' if ok else 'MISMATCH'}]")


if __name__ == "__main__":
    main()

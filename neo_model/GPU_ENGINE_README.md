# GPU-native DFT / TDDFT / NEO-DFT engine

GPU-accelerated electronic-structure for the Venus/A206 chromophore pipeline,
on the local **RTX 4080 (16 GB)**. Built on **gpu4pyscf** (GPU-native DFT/TDDFT,
a TeraChem-class engine) bridged to the **PySCF NEO fork** (quantum nuclei).

## What's GPU-native here
| capability | engine | status |
|---|---|---|
| DFT SCF | gpu4pyscf | exact vs CPU |
| TDDFT (excitations) | gpu4pyscf | exact vs CPU |
| **NEO-DFT (quantum proton)** | gpu4pyscf electronic veff + NEO fork | **exact vs CPU, faster** |

NEO-DFT is the novel piece: **no public GPU NEO code exists** (TeraChem has no NEO
at all). The dominant cost of NEO-DFT is the *electronic* structure; we route the
full electronic veff (J+K+XC) to the GPU via `gpu_accel.gpu_neo_hook`, while the
tiny quantum-proton SCF and electron-proton coupling stay on the CPU. The
multicomponent physics is identical (energy matches CPU NEO to <0.01 µHartree).

## Files
- `gpu_accel.py` — library: `gpu_scf`, `gpu_tddft`, `gpu_optimize`, `gpu_neo_hook`,
  `enable_managed_memory`.
- `gpu_neo.py` — CLI engine (`dft` | `tddft` | `neo`).

## Use
```bash
conda activate pyscfneo
export LD_LIBRARY_PATH=/home/robson/anaconda3/lib:$LD_LIBRARY_PATH   # MKL
export MKL_THREADING_LAYER=GNU
cd neo_model

python gpu_neo.py dft   --geom chromo_his.xyz --charge -1
python gpu_neo.py tddft --geom chromo_his.xyz --charge -1 --nstates 6
python gpu_neo.py neo   --geom chromo_his.xyz --charge -1 --quantum 36 --epc 17-2 --verify
```
`--verify` (neo) also runs CPU NEO and asserts the energy matches (<5e-5 Eh).

## Validated results (40-atom chromophore+His148, b3lyp/6-31G*, charge -1)
- GPU DFT  : E = -1118.75389839 Eh  (matches CPU)            ~25 s
- GPU TDDFT: bright pi->pi* 3.106 eV = 399.2 nm, f=1.075 (matches CPU)  ~62 s
- GPU NEO  : E = -1118.71577311 Eh  (dE vs CPU = 0.0 uEh)    2.3x faster (40.8s vs 94.5s CPU)

## Important notes / honest scope
- **This is not a from-scratch TeraChem clone.** gpu4pyscf *is* a mature, GPU-native
  DFT engine competitive with TeraChem for DFT/TDDFT; we wrap it and add the NEO
  bridge. Re-implementing GPU integral kernels from scratch would be slower and
  buggier than gpu4pyscf — the value here is the **GPU NEO-DFT bridge**.
- **GPU speedup is system-size dependent.** At 40 atoms NEO is ~2-2.5x (it has
  irreducible CPU parts: nuclear SCF, e-p coupling). The GPU advantage grows with
  the electronic system size / basis; for large QM regions it will be much larger.
- **Full TDDFT memory**: the non-Hermitian eigensolver for this anion peaks past
  the 4080's ~15 GB free. `enable_managed_memory()` (called automatically in
  `gpu_tddft`) routes >100 MB allocations to CUDA *managed* memory so they spill to
  host RAM instead of OOM-ing. On a 40/80 GB datacenter GPU this is unnecessary.
- Always run NEO with `--verify` on a new system before trusting production runs.

## Environment provenance
- pyscf 2.13.1 = theorychemyang NEO fork, built from source + installed as a real
  package in site-packages (NOT editable: editable + namespace plugins conflict).
- gpu4pyscf-cuda12x 1.7.3 + gpu4pyscf-libxc-cuda12x + cupy-cuda12x 14.1.1
  (CUDA-12 wheels run fine over the CUDA-13 driver). pyscf-dispersion layered in.

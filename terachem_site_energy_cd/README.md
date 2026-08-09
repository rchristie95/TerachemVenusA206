# Frame-resolved TeraChem site energies and nondegenerate CD

This directory is a restartable, fail-closed workflow for two equivalent
44-atom site calculations per MD frame. It does not edit either manuscript and
does not permit a compact `J_i` table to be joined to a trajectory unless the
trajectory SHA-256 matches the archived production source.

## Current gate

The available full-solvent 1,000-frame rerun is structurally usable for a smoke
test (69,456 atoms; 29 atoms per CR2; frames saved every 500 MD steps), but its
SHA-256 is `d854e18b...`, not the archived production hash `5da1f8b...`.
Consequently `results/run_manifest.json` sets `production_join_allowed=false`.
No 20-frame, 100-frame, or 1,000-frame calculation is authorized from these
inputs.

The smoke input uses wB97X-D3/6-31G*, TDA, seven roots, charge/multiplicity
`-1/1`, no PCM, a fixed frame-0 12 A embedding selection, conserved
partial-residue charges, and the complete RESP-charged partner CR2 in the MM
field. The installed TeraChem manual documents electric transition dipoles but
not magnetic transition dipoles or rotational strengths, so the spectral code
reports relative interaction-induced exciton-chirality CD only.

## Exact smoke commands

Run with the TeraChem Conda Python:

```bash
TC_PY=/home/robson/anaconda3/envs/TeraChem/bin/python
TC_BIN=/home/robson/Desktop/TeraChemPython/TeraChem/bin/terachem
TOP=tc_tandem_nvt_rerun_20260722_retry2/solvated_protonated.pdb
TRAJ=tc_tandem_nvt_rerun_20260722_retry2/tandem_nvt_1000.dcd

$TC_PY terachem_site_energy_cd/scripts/audit_inputs.py \
  --topology "$TOP" --trajectory "$TRAJ" \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --terachem "$TC_BIN" --output terachem_site_energy_cd/results/run_manifest.json
```

The audit exits 2 on the present expected hash mismatch while still writing the
manifest. Preparing frame 0 creates the fixed MM atom inventory; frame 499 then
reuses it:

```bash
$TC_PY terachem_site_energy_cd/scripts/prepare_production_frame.py \
  --topology "$TOP" --trajectory "$TRAJ" --frame-index 0 \
  --output-dir terachem_site_energy_cd/results/frame_0000_selection \
  --embedding-cache terachem_site_energy_cd/results/embedding_charges.npz \
  --fixed-embedding-selection terachem_site_energy_cd/results/fixed_embedding_selection.json \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --retain-partner-cr2-charges --conserve-boundary-residue-charge

$TC_PY terachem_site_energy_cd/scripts/prepare_production_frame.py \
  --topology "$TOP" --trajectory "$TRAJ" --frame-index 499 \
  --output-dir terachem_site_energy_cd/results/smoke_frame_0499 \
  --embedding-cache terachem_site_energy_cd/results/embedding_charges.npz \
  --fixed-embedding-selection terachem_site_energy_cd/results/fixed_embedding_selection.json \
  --amber-cr2-prmtop anionic_build/monomer_solv.prmtop \
  --retain-partner-cr2-charges --conserve-boundary-residue-charge

$TC_PY terachem_site_energy_cd/scripts/launch_jobs.py \
  terachem_site_energy_cd/results/smoke_frame_0499 --terachem "$TC_BIN" --gpu 0
```

The launcher runs one independent TeraChem process at a time on the explicitly
selected GPU, preserves incomplete logs, and reuses only outputs containing
`Job finished:`. Run one launcher process per GPU on disjoint frame directories
for bounded scale-out.

## Scale-up gate and cost

After restoring the exact archived trajectory/topology, rerun the audit and
require `production_join_allowed=true`. Then generate transition densities or
NTOs for every candidate root in the smoke frame, track the physical state by
overlap, repeat a cold-start/restart comparison, and require an unambiguous
assignment before preparing the deterministic 20-frame stratified pilot.

The two energy jobs used 85.8 GPU-seconds total, or 0.0238 GPU-hours per frame.
Energy-only lower bounds are therefore about 2.4 GPU-hours for 100 frames and
23.8 GPU-hours for 1,000 frames on one RTX 4080. Density/NTO tracking, basis
checks, cold starts, and failed retries add cost; their multiplier must be
measured in the validated 20-frame pilot rather than guessed.

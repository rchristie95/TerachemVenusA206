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
field.

The installed TeraChem *manual* documents electric transition dipoles but not
magnetic ones. The *output* does print them -- see the audit below -- so do not
repeat the inference that they are unavailable. They are simply not usable.

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

## Splitting-profile investigation (no electronic-structure launch)

The generalized two-site identifiability and line-shape sweep is independent
of the blocked TDDFT production join:

```bash
python terachem_site_energy_cd/investigate_splitting_profile.py
```

It reconstructs the published constrained dVenus Lorentzian model, maps the
exact `J/Delta` ridge, tests the corrected 1000-frame coupling/geometry
ensemble, and quantifies differential-disorder sensitivity. See
`SPLITTING_PROFILE_INVESTIGATION.md` and
`results/splitting_profile_investigation/`. It makes no absolute molar-CD
claim and does not pair the archived coupling rows with an unmatched TDDFT
trajectory.

## Magnetic transition dipoles: printed, but not usable

Every completed `results/v2_camb3lyp_frame_*/site_*/tddft.out` carries

```
Magnetic transition dipole moments and rotational strengths:
  Root   Lx   Ly   Lz   R(len)   R(vel) (a.u.)
```

`scripts/audit_magnetic_dipoles.py` decodes the conventions by numerical
identity over all 1540 roots of all 220 runs (max residual ~1e-4):

    R(len) = -mu_len . L          R(vel) = (mu_vel . L) / omega

so `L` is the raw angular-momentum matrix element with no 1/2c, both gauges are
built from the *same* `L`, and the gauge origin is molecule-centred (|L| ~ 1.4
a.u.; it would be ~20 a.u. about the coordinate origin, which sits ~112 A away).

The numbers fail three ways:

| diagnostic | value | exact value |
|---|---|---|
| \|mu_vel\| / (omega \|mu_len\|) | 0.233 | 1.0 |
| angle(mu_vel, mu_len) | 129.4 +- 6.8 deg | 180 deg |
| f_len / f_vel | 19.2 (1.13 vs 0.062) | 1.0 |
| R(len), R(vel) same sign | 12 / 220 runs | 220 / 220 |
| origin shift for a 100 % change in R(len) | 0.79 A | infinite |

The bright state carries essentially all the chromophore absorption intensity,
so f_vel = 0.06 is untenable: the velocity-gauge quantity is the broken one.
The origin dependence is the same defect seen a second way, because a gauge
shift by `a` moves the length-gauge value by `a . (mu_vel x mu_len)`, which
vanishes only when the hypervirial relation `mu_vel = -omega mu_len` holds.
Either failure alone disqualifies these numbers.

**It is not a TeraChem bug.** `scripts/pyscf_gauge_crosscheck.py` repeats the
bright state in PySCF on the same geometry and functional, also TDA, gas phase
(MM point charges are a local multiplicative potential and cannot break the
hypervirial relation):

| | TeraChem 6-311+G** | PySCF 6-31G* | PySCF 6-311+G** |
|---|---|---|---|
| \|mu_vel\| / (omega \|mu_len\|) | 0.233 | 0.229 | 0.243 |
| f_len / f_vel | 19.2 | 19.1 | 17.0 |
| f_len | 1.13 | 1.08 | 0.87 |

Two codes agree to within the ensemble scatter, and the defect does **not**
improve on adding diffuse functions. The cause is the method: TDA discards the
de-excitation block the off-diagonal hypervirial relation requires, and the
range-separated exchange operator does not commute with `r`. If usable
rotational strengths are ever wanted from these runs, the thing to change is
TDA (`cis yes`), not the code. That has not been tested here.

## Absolute CD against Nguyen 2025 Table S3

`scripts/absolute_cd_vs_experiment.py`. None of the above blocks absolute CD,
because the exciton-chirality couplet needs only the two electric transition
dipoles and their displacement. The comparison is anchored on the couplet first
moment, which is independent of detuning, homogeneous linewidth and
inhomogeneous spread because the `1/Omega` in `R_pm` cancels the `Omega` lever
arm:

    M = int (Delta_eps / nubar)(nubar - nubar_0) d nubar = -pi nubar_0 J T / 2.296e-39

Result over 95 usable production frames:

| quantity | predicted | experiment (Table S3) |
|---|---|---|
| couplet first moment | +1.51e7 | -3.70e6 |
| peak-to-peak [Theta] | 2.94e5 | 2.83e5 |
| dissymmetry factor g | 4.83e-4 | 4.65e-4 |

The magnitude is right to a factor of 4.1 on the invariant measure; the
**handedness is opposite**. The peak-to-peak agreement is partly fortuitous --
our couplet is spread over `<Omega>` = 584 cm^-1 against the observed 262, and
the lost peak amplitude offsets the fourfold excess in M.

Sign provenance, since the handedness is the headline: the structures are
natural L-proteins (93.4 % of backbone phi negative, median -97 deg), so the
lab frame is right-handed; the chirality triple product is negative in 95/95
frames here and 69.6 +- 11.7 (1000/1000 negative) in the independent
transition-density-centroid ensemble; and the stored `J_cm` is positive in
1000/1000 frames. The couplet handedness is therefore not a sign convention.

Two traps worth recording. The stored `J_cm` = 32.82 and `J_pda_cm` = 27.64 are
already screened by epsilon = 1.77 (the unscreened point-dipole value
recomputed from `coupling_geometry.npz` is 48.93 = 27.64 x 1.77). And the
experimental `[theta]` normalisation is per mole of *construct*, fixed by the
exact 2:1 ratio of eps_516 for dVenus-TD (184,400) and dVenus-TDX (92,200).

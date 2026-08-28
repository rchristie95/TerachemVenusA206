# Numerical report: tandem-Venus TDDFT ensemble decoherence

> **Superseded in part, 2026-08-27.** The 4 fs campaign is bit-reproducible; site B alone was anomalous and is corrected by the bright-moment estimator. Homogeneous HWHM is now 346 cm^-1. See `CORRECTIONS_2026-08-27.md`.

Date: 22 July 2026  
Repository: `/home/robson/PetaChem`  
Commit: `34775aa82847f654d77f8720182ddc36782095c7` (`main`)

## Outcome

The same-geometry reference pilot, calibration, restart manifests, transition-
density compaction, and mutual-TDC validation are complete. The requested
3,000-frame ensemble calculation is not complete because its exact raw
coordinate trajectories are absent from the local filesystem and deliberately
excluded from Git. No frame was substituted, interpolated, or silently marked
complete.

This is a source-data blocker, not an electronic-structure failure. Both short
DCD names are retained only in the compact solvation summaries, and the latest
long production package retains only coupling/geometry data. These compact
products cannot reconstruct instantaneous 44-atom QM coordinates, stacked-Tyr
coordinates, solvent coordinates, or the MM point-charge field.

## Repository update and preservation

The initial checkout was branch `agent/recover-orca-paper-pipeline` at
`dced07dbf2d403b7aadcf6105cac51e090539dce` with pre-existing untracked files
under `neo_model/`. The requested commands were run in order:

```text
git status --short --branch
git fetch origin
git switch main
git pull --ff-only origin main
```

`main` fast-forwarded from `b11f87c` to `34775aa`. The pre-existing untracked
files remain in place; nothing was reset, discarded, overwritten, or stashed.

## Required reading and inherited conventions

The solvation note, handoff, numerical test report, reproducibility README,
manuscript, TeraChem inputs/outputs, QM-region builders, trajectory preparation,
GPU TDC, density-masking, and analysis programs were inspected. They establish:

- two contiguous short segments of 1,000 frames each, sampled every 4 fs;
- 2,000 short frames over 0--7,996 fs;
- 1,000 long production frames at 1 ps intervals over 1--1,000 ps;
- two CR2 sites at topology residues 66 and 328 and stacked Tyr residues 202
  and 464;
- a 44-atom QM region: 29-atom anionic CR2, 12-atom Tyr phenol group, and three
  hydrogen link atoms;
- 41 physical QM atoms and explicit exclusion of link/cap-associated density;
- `wB97X-D3/6-311G**`, charge -1, singlet, point-charge embedding;
- optical coupling screening by `epsilon_opt = 1.77` in the audited production
  coupling (the manuscript describes the same convention as approximately
  1.77--1.78);
- whole-residue periodic translations only for production imaging, without
  coordinate fitting or deformation.

The short topology audit records 69,471 atoms: 7,761 protein atoms, 61,587
water atoms, and 123 ions. The two short sites are CR2 66/Tyr 202 and CR2
328/Tyr 464. Periodic box vectors are required by the retained minimum-image
and PME analyses, but their exact values cannot be re-read without the missing
DCD headers/topology.

## Rejected trajectory substitution

The local July 6 file `tc_tandem_nvt_1000/tandem_nvt_1000.pdb` contains 1,000
models and a 90.183 A cubic box, but it has 69,418 atoms, 61,524 water atoms,
and 31 atoms in each CR2. The exact July 21 production workflow uses the
retained 29-atom AMBER CR2. This older file is therefore not the requested
trajectory and was rejected rather than silently substituted.

The compact long file `coupling_geometry.npz` has transition-dipole and centroid
geometry but not the full protein/water coordinates. The protein-only clean PDB
also lacks water and ion coordinates. Neither can supply the requested
instantaneous electrostatic embedding.

## Static GPU TeraChem pilot

The reference calculation was rerun in an independent directory using the
identical retained 44-atom geometry and 2,350 point charges. The TeraChem and
ORCA point-charge numerical rows are identical. Hardware and software were:

- NVIDIA GeForce RTX 4080, GPU 0;
- driver 580.126.09, CUDA 13.0 reported by the driver;
- TeraChem 1.97B-beta-251105, Git
  `3fcf25291f21a40754a20c9c6450f71e006a2729`;
- 48 s wall time, 45.79 s TeraChem processing time;
- 1,734.15 GPU mega-words reported available by TeraChem.

SCF converged, and the five-root Davidson solve converged in 15 iterations.
The relevant roots were:

| Root | Energy (eV) | f | |mu| (a.u.) | Dominant configuration |
|---:|---:|---:|---:|---|
| 1 | 3.08312196 | 0.0065 | 0.2935 | 93 -> 94 |
| 2 | 3.43441587 | 1.0228 | 3.4864 | 93 -> 95 |
| 3 | 3.45608312 | 0.2634 | 1.7637 | 92 -> 94 |
| 4 | 4.02299611 | 0.0037 | 0.1946 | 85 -> 94 |
| 5 | 4.14137713 | 0.0682 | 0.8198 | 92 -> 95 |

Root 2 reproduces the retained bright pi-pi* assignment and reference output to
printed precision. Root 3 is close in energy and non-negligibly bright, so a
fixed-root assumption would be unsafe. Five roots are a provisional minimal
pilot choice; beginning/middle/end frames from both trajectories must confirm
that this root count avoids root loss before production.

The production state tracker is designed to use the cap-masked physical-atom
transition-density overlap first, then configuration character, oscillator
strength/dipole orientation, and finally adjacent-frame energy continuity in
the short window. The overlap also fixes the arbitrary transition-density
phase so signed couplings are comparable from frame to frame.

## STEOM calibration

The final ORCA STEOM electric-dipole reference is 2.366633 eV, 523.90 nm, and
`f = 0.885045444`. With the reproduced TeraChem reference energy of
3.43441587 eV,

```text
Delta = E_STEOM - E_TDDFT = -1.06778287 eV.
```

The initial calibration is `E_corrected = E_TDDFT + Delta`. A common constant
shift changes the absolute spectral origin but cancels from `E_A-E_B`; it must
leave the differential-gap correlation, cumulant coherence, and T2* unchanged.
That invariance cannot yet be tested on frame data. The shift is not applied to
the transition density or coupling. No affine slope was fit because only one
matched final STEOM/TDDFT geometry is available, so a slope is not identifiable.

## Mutual TDC validation

A rectangular FP64 OpenCL kernel was added because independent A and B TDDFT
densities can have different retained grid sizes. It reports vacuum and
optically screened signed TDC, PDA diagnostics, separation, dipole angle,
normalization diagnostics, A/B exchange symmetry, and a brute-force CPU subset.

On the retained static TeraChem root-2 density, the physical-density builder:

- removed grid points whose nearest QM atom was one of link atoms 42--44;
- recorded the cap-mask neutrality error;
- restored exact transition-density neutrality by minimal separate rescaling
  of the positive and negative lobes;
- restored the printed TeraChem transition-dipole projection;
- retained 268,038 compact points.

The first CPU validation used the first 1,000 file-order points. Its integral
was nearly zero through cancellation, making a relative-error criterion
ill-conditioned; that failed output is retained. Retry 1 used the 1,000 largest
absolute-density points from each site and passed:

- A->B versus B->A relative difference: `1.65e-15`;
- CPU versus GPU 1,000-by-1,000 subset difference: exactly zero at printed
  double precision;
- full FP64 GPU runtime: 9.54 s.

The resulting static diagnostic is 217.79 cm-1 in vacuum and 123.05 cm-1 after
division by 1.77; PDA values are 51.79 and 29.26 cm-1. These are validation
numbers on a retained reference placement, not short- or long-window ensemble
statistics, and they are not used to replace the retained STEOM production
coupling.

## Restart and failure handling

`manifests/short_frames.csv` has 2,000 rows and `manifests/long_frames.csv` has
1,000 rows. Every unavailable coordinate row is explicitly classified as
`blocked_missing_source` / `missing_source_trajectory`. The idempotent
initializer preserves every existing classified row, including failures and
their retry counts, and never launches work or creates duplicate frame IDs.

No SCF, excited-state, root-tracking, point-charge-format, GPU-memory, or
numerical frame failure has occurred because no requested frame job was
launched. When exact inputs are restored, original failed job directories must
be retained, retries must use minimally stronger settings, and no missing
energy or coupling may be interpolated.

## Statistical and decoherence status

No short/long energy or coupling histogram has been generated because there
are no new framewise TDDFT site energies. Likewise, no TDDFT differential-gap
correlation, cumulant T2*, direct two-site propagation, affine-calibration test,
block bootstrap, or long-window inhomogeneous width can be reported yet.

It would be scientifically invalid to reuse the approximately 1 ps long-window
cadence for a 20--60 fs decay. The retained comparison values remain only the
previous fixed-charge results: 21.8 fs minimum image, 24.0 fs PME, block median
24.5 fs with SD 3.7 fs, and the manuscript's phenomenological 60 fs sensitivity
value. There is no new geometry-dependent TDDFT T2* to say whether those values
are supported, lengthened, or shortened.

## Reduced surrogate feasibility rerun

After the requested scope was reduced, beginning/middle/end snapshots were
streamed from the only local full-solvent 1,000-model PDB. This remains a
surrogate: it has 69,418 atoms and the obsolete 31-atom CR2 rather than the
audited 69,471-atom/29-atom production system.

The preflight rejected every instantaneous surrogate chromophore. The
19-heavy-atom fit RMSDs were 2.93--3.07 A at the initial frame and increased to
7.22--9.74 A at the middle/end frames, compared with the retained production
maximum of 0.705 A. No TeraChem job used those malformed geometries.

One initial-frame electrostatic feasibility test was run by rigidly placing the
validated reference CR2/Tyr stacking geometry from the instantaneous stacked
Tyr heavy atoms. The Tyr fits were 0.056 A at site A and 0.037 A at site B.
Both five-root energy calculations converged in about 54 s. Physical-state
tracking selected different numerical roots:

| Site | Root | Energy (eV) | f | Reference-density score |
|---|---:|---:|---:|---:|
| A | 3 | 3.39905601 | 1.2440 | 0.9023 |
| B | 2 | 3.42337285 | 1.2699 | 0.9788 |

Both are dominated by the intended `93 -> 95` character. The differential gap
is -0.02431684 eV (-196.13 cm-1). Adding the common STEOM correction gives
2.33127314 and 2.35558998 eV but leaves that gap exactly unchanged.

Selected-root transition-density jobs used only three roots for site A and two
for site B and converged in 42 and 32 s. Link-cap density was excluded, density
phase was set by reference overlap, and neutrality was restored explicitly.
The first TDC integral retained the raw periodic images and exposed an invalid
68.9 A separation. That output is retained. A whole-site B translation of
`[0, 0, 90.183] A` produced a 34.61 A surrogate separation and passed:

- signed vacuum TDC: -32.57 cm-1;
- screened TDC (`epsilon = 1.77`): -18.40 cm-1;
- signed vacuum/screened PDA: -9.62/-5.43 cm-1;
- A/B symmetry relative error: `1.28e-14`;
- CPU/GPU subset relative error: `1.73e-15`.

These numbers validate convergence, state tracking, periodic-image handling,
and signed mutual integration only. The old environment, separation, and
chromophore history are not the requested production ensemble. A single rigid
surrogate frame cannot yield T2, disorder distributions, or a comparison with
the prior 21.8--24.5 fs estimates.

### Root-cause correction and valid static rerun

Continued debugging located the obsolete-trajectory failure upstream. Its
generated generic CR2 XML contains only 12 X--H bonds per site and omits the
heavy-atom covalent network. The solvated pre-minimization structure has 0.296
A CR2 heavy-atom fit RMSDs at both sites. Current topology inference adds the
40 missing bonds, giving 64 CR2-connected bonds across two sites. By contrast,
the obsolete relaxed structure retains only 24 and has 8.74--9.83 A fits.
Those dynamics cannot be repaired into an ensemble and were not rerun.

`run_nvt.py` now has a fail-closed CR2 connectivity guard. In regression tests,
the valid initial structure forms one 31-atom component with 32 internal bonds
at each site after inference. The corrupted NVT frame is rejected with 19
components before minimization or dynamics can start.

The valid pre-minimization structure was used for an improved static
diagnostic with its actual instantaneous heavy-atom coordinates. Since CR2 is
wholly QM, both CR2 residues were removed before standard AMBER14/TIP3P MM
charge construction. This charge route was checked against the retained AMBER
monomer topology: the maximum per-atom discrepancy is `2.20e-9 e`.

| Site | Root | Energy (eV) | f | Reference-density score |
|---|---:|---:|---:|---:|
| A | 3 | 3.13677405 | 0.9722 | 0.9394 |
| B | 3 | 2.97345490 | 0.8058 | 0.9454 |

Site B also has a competing root 4 at 3.00266487 eV with `f = 0.2814`. The
selected states retain `93 -> 95` character. The static differential gap is
0.16331915 eV (1317.26 cm-1); the common STEOM shift leaves it unchanged.

The signed TDDFT-density coupling is 179.42 cm-1 in vacuum and 101.37 cm-1
after optical screening. The separation is 23.39 A and the transition-dipole
angle is 106.42 degrees. A/B exchange symmetry (`2.63e-15`) and the CPU/GPU
subset check (`1.72e-16`) pass. This coupling scale is consistent with the
retained production ensemble mean of 117.19 cm-1, but the single obsolete
initial structure is not a trajectory distribution and cannot determine T2.

Spatial-density tracking of all first-three roots confirms embedded root 3 at
both sites: the mean bidirectional overlaps are 0.522 and 0.452, while all
lower embedded roots are at most 0.136. The overlap margins are 0.388 and
0.316, comfortably above the 0.10 ambiguity guard. Gas-phase state mixing is
substantial. The leading branches are site-A root 1 at 2.92354498 eV and
site-B root 3 at 3.16891755 eV, unchanged over interpolation widths from 0.08
to 0.20 A. Their margins are only 0.054 and 0.081, however, so both gas
assignments are explicitly classified as ambiguous.

Protein-only fields retain root 3 and give 3.13276511/2.97512038 eV at A/B,
versus 3.13677405/2.97345490 eV with protein and nearby waters. Thus, for this
one structure:

- gas-phase differential gap: -1979.06 cm-1;
- protein-only differential gap: +1271.49 cm-1;
- full differential gap: +1317.26 cm-1;
- protein-induced change relative to the leading gas branch: +3250.55 cm-1;
- incremental water contribution: +45.77 cm-1.

The protein field clearly drives state reordering. The quoted gas-to-protein
gap inversion is a branch sensitivity because gas assignments are ambiguous;
the six/four nearby-water increment between robust embedded root-3 states is a
smaller correction. Because this is a single pre-minimization structure, it
does not establish the corresponding trajectory averages or their role in T2.

## Exact unblock requirement

Restore or provide local paths to all of the following, with their matching
topology/system metadata:

1. `full_highres_0to4ps_4fs.dcd` (1,000 frames, 69,471 atoms);
2. `full_highres_4to8ps_4fs.dcd` (1,000 frames, 69,471 atoms);
3. the exact July 21 long full-system production DCD whose recorded SHA-256 is
   `5DA1F8B2CE814DD04467B4D121C3BA70CCC8619852D045120D2E73086DF53E5E`,
   or its coordinate-equivalent PBC-whole DCD with recorded SHA-256
   `C8285024076091EB14BCEB2C25A93CD73BAB7FB5675CBA41BF42E7C064425C93`;
4. the topology/system files matching those trajectories.

Once restored, the first action is a strict frame/atom/box/checksum preflight,
then beginning/middle/end site calculations for both windows. Production must
remain disabled until state assignment, transition-density output, memory,
runtime, and the minimal root count pass that pilot.

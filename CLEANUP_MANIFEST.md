# Cleanup manifest — 2026-08-11

Produced by a four-way read-only inventory of the 429 GB tree. Every size below
was measured directly, not estimated. Nothing in the FLAGGED section has been
deleted.

## Already done

| action | reclaimed | evidence |
|---|---|---|
| Deleted 235 ORCA `.tmp` scratch files from the crashed `steom_phenol_svpd_robust2` DLPNO-STEOM-CCSD run | **251 GB** | run crashed at `Making the (pseudo)densities` (MPI signal 11) *after* printing its spectrum; 527.3 nm / f = 0.8467 preserved in the 324 kB `.out`, bit-identical to `_robust`. 0 tracked files touched. |
| Rescued 6 harvested ensemble arrays from a `/tmp` scratchpad into `terachem_site_energy_cd/results/ensembles/` + README | — | `ens_v2.npz` (n=88), `ens_final.npz` (n=87), `ensemble_camb3lyp.npz` (n=23), `ct.npz`, `pol_v2.npz`, `surrogate_fit.npz`; md5-verified copies |
| Added a `.gitignore` negation so those arrays are committable | — | the global `*.npz` rule (line 42) was hiding them |

`neo_model/` went 252 GB -> 735 MB. Disk: 189 GB free -> 440 GB free.

## Two correctness issues found while surveying

These are not tidying items; they affect the manuscript.

1. **`manuscript/JPCB_tandem.tex` still carries the pre-fix coupling.** It states
   `J = 117.2 +/- 5.5 cm^-1`, splitting `234.4`, PDA `27.6`, for the same
   1000-frame ensemble that `README.md` and
   `coupling_nvt_production_cr2_1000_20260721/` both report as `32.82 +/- 1.55`.
   `117.2 / 32.82 = 3.570993` against the Angstrom/Bohr bug factor
   `1/0.529177^2 = 3.571068` — this file never received the round-2 correction.
   It is **tracked and modified**. Its abstract also quotes the 523.90 nm site
   energy rather than the reproducible 527.3 nm. Correct it or retire it.

   *Update 2026-08-17:* `JPCB_tandem_round_2.tex` is now the single live
   manuscript. `tandem_dimer_2.tex` (the alternative single-author framing,
   "Site-Energy Detuning Explains the Apparent Excitonic Coupling…") has been
   **deleted**; it is tracked history and recoverable with
   `git show <rev>:manuscript/tandem_dimer_2.tex`.

2. **Round-2 source tracking.** `notes/J_apparent_derivation.tex` was untracked
   as of this survey. Also untracked in
   `terachem_site_energy_cd/`: `absolute_cd.py`, `polarizable_embedding.py`,
   `vibronic_exciton.py`, `cumulant_lineshape.py`,
   `PROVENANCE_v1_DEFECTIVE_MD.md`, and 6 files under `scripts/`.

## Groups A-E — APPROVED AND EXECUTED 2026-08-11

All five groups below were approved and deleted. Verified afterwards: **0 tracked
files removed**, and every KEEP item (ORCA install, v1/v2/v3 trajectories,
`_authoritative_backup/`, `run1_nto/`, `qchem_validation/`, both nested repos,
`coupling_nvt_production_cr2_1000_20260721/`, `oqs_out/`, `multipole_out_correct/`,
`lineshape_tandem/`, the rescued ensembles, `reference_pilot/`, `calibration/`,
`short100fs50_fixedmm/`) confirmed present.

One deviation from the plan, on new evidence: `resp_cr2/` was **stripped rather
than deleted**. `cr2_resp_charges.dat` and the omega-tuning scan exist nowhere
else on the system, so the 293 MB of TeraChem `.molden` conformer scratch was
removed and the 2.4 MB parameter record kept (297 MB -> 4 MB). `cr2.prepi` turned
out to be 0 bytes — a failed output, not a deliverable.

Ordered by confidence. Group A was provably redundant; E was judgement.

### A. Verified duplicate — 11.7 GB

| path | size | why |
|---|---|---|
| `wb97xd3_100fs100_complete_package/` | 6.2 GB | Unpacked handoff bundle. Its `frames/short100fs100_wb97xd3/` is byte-identical to the live copy — 4508 files / 6.10 GB on both sides, md5 matches on spot checks. Fully reconstructible from the zip. |
| `wb97xd3_100fs100_complete_package.zip` (+`.sha256`) | 5.5 GB | Sealed archive of the above; `sha256sum -c` verified **OK**. Third copy of data already live in the tree. Consider cold storage rather than deletion if it is owed to a collaborator. |

Deleting the unpacked directory alone is the zero-risk move: the zip is intact
and checksummed.

### B. Superseded MD lineage — 17.0 GB

| path | size | why |
|---|---|---|
| `tc_dimer_nvt_restrained_1500/` | 7.5 GB | Crystal-contact **dimer** lineage, abandoned when the project moved to the tandem construct |
| `tc_dimer_nvt_long/` | 5.3 GB | same |
| `tc_dimer_nvt/`, `tc_dimer_nvt_restrained/` | 2.2 GB | same |
| `dimer_nvt_restrained_1500_clean.pdb`, `dimer_nvt_restrained_clean.pdb` | 1.0 GB | filtered text derivatives of the above |
| `tc_tandem_nvt/` + `tandem_nvt_clean.pdb` + `tandem_nvt_whole.pdb` | 0.3 GB | pre-v1 shakedown, zero-charge CR2 template |
| `tc_tandem_nvt_rerun_20260722_retry2/` | 0.8 GB | superseded by v2 |
| `tc_tandem_nvt_rerun_20260722{,_retry1}/`, `tc_tandem_nvt_highcadence_20260722/`, `tc_tandem_nvt_100ns_50frames_20260723/` | small | empty or crashed (the "100ns 50 frames" DCD holds **3** frames) |
| `tc_qmmm_opt/` | 0.2 GB | unconstrained QM/MM opt, diverged globally; superseded by `tc_qmmm_opt_constrained/` |

**Keep** `tc_tandem_nvt_1000/` (5.3 GB) and its two derived PDBs (1.2 GB)
despite v1 being defective — `reference/SHA256SUMS` tracks
`tandem_nvt_1000_whole.pdb` by hash, and `coupling_tandem_1000/` is built from
`tandem_nvt_1000_clean.pdb`.

A hypothesis worth recording as **false**: these giant `.pdb` files are *not*
redundant text copies of `.dcd` binaries. No `tc_dimer_nvt*` or pre-v2
`tc_tandem_nvt*` directory contains a `.dcd` at all — DCD only became the output
format on 2026-07-22. Each `_clean`/`_whole` file is a different atom selection
of its parent, not a format duplicate.

### C. Superseded QM ensembles in `tddft_ensemble_decoherence/` — 14.4 GB

| path | size | why |
|---|---|---|
| `frames/short100fs50/` | 13.0 GB | computed before the MM-embedding fix; redone as `frames/short100fs50_fixedmm/` (6.8 GB, kept). `finalize_fixedmm_100fs.py` exists specifically because this was redone. |
| `surrogate_pilot/` | 1.9 GB | its own README rejects it — built on the obsolete v1 CR2 topology (31 atoms/CR2), "not production ensemble data" |
| `rapid_100fs50/` | 0.5 GB | analysis products of the pre-fix window above |

### D. Superseded site-energy families — 9.8 GB

All computed on the **v1 defective trajectory** and named as superseded in
`PROVENANCE_v1_DEFECTIVE_MD.md`. Safe now that `ens_final.npz` (n=87) — the
harvest of `hi_camb3lyp_frame_*` — has been rescued into the repo.

| family | dirs | size |
|---|---|---|
| `hi_camb3lyp_frame_*` | 160 | 6.4 GB |
| `linkonly_frame_*` | 206 | 2.1 GB |
| `ct_camb3lyp_frame_*`, `ct_frame_*`, `fullsys_frame_*`, `hi_wb97xd3_frame_*`, `hi_wb97xv_frame_*`, `pilot_frame_*`, `shell{16,20,25,30}_frame_0499`, `guesstest_frame_*` | 40 | 1.0 GB |
| `production_frame_*` | 993 | 0.3 GB — **zero of 993 were ever run**; pure unrun input prep against the defective v1 trajectory |

**Keep** `v2_camb3lyp_frame_*` (121 dirs, 5.6 GB, 94% complete) — the current
campaign, harvested as `ens_v2.npz` but still live. Keep `v2_linkonly_frame_*`
(200 dirs, 1.4 GB) — queued, not yet launched.

### E. Regenerable output and transients — 0.9 GB

| path | size | why |
|---|---|---|
| `videos/` | 0.5 GB | regenerable from `render_*.py`; the manuscript links a YouTube playlist, not these files |
| `resp_cr2/` | 0.3 GB | AmberTools RESP fitting scratch; deliverables already consumed into `anionic_build/`. **Confirm nothing downstream still reads it.** |
| `steom_densities.pse` | 29 MB | saved PyMOL session, regenerable |
| `ScienceDirect_files_11Aug2026_15-31-32.068.zip` | 17 MB | today's browser download, unrelated to the pipeline |
| `__pycache__/` | 1.6 MB | bytecode cache |
| `oqs_out_new/` | 532 kB | **older than `oqs_out/` despite the name** (07-05 vs 07-06); `oqs_out/` is what the tracked manuscript figures come from |
| `multipole_out/` | 28 kB | pre-correction; superseded by `multipole_out_correct/` |
| `coupling_sampling_steom_out/` (155.82), `..._matched_out/` (96.38), `..._50rand_out/` (96.80), `..._expdip_out/`, `coupling_paper_steom_thermal_1500/`, `coupling_tandem_unrestrained/`, `lineshape_out/` | small | all carry pre-unit-fix coupling values (the 96 / 155 / 193 family) |
| 21 top-level `*.log` incl. 3 AmberTools install transcripts and 2 zero-byte files | 350 kB | transcripts of completed one-off runs |

## Stale documentation — keep for history, do not cite

Four untracked handoff notes assert pre-correction couplings and should be
marked or moved to an archive directory rather than deleted:
`STEOM_COUPLING_FINDINGS_2026-06-29.txt` (~155-160), `BRIEFING_FOR_LLM.txt`,
`NEXT_AGENT_INSTRUCTIONS.txt`, `WHY_ORCA_NOT_TERACHEM.txt` (all "J ~ 65-94"),
plus `ResearchPaper_origin.tex`. `README.md`, `ORCA_ONLY_PIPELINE.md` and
`LINUX_TERACHEM_SITE_ENERGY_CD_PROMPT.md` are current and correct (32.82 / 65.6).

## Not to be touched

- `orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/` (17 GB) — **load-bearing**.
  `reproduce_paper.py:87` resolves the ORCA binary from inside the repo, and
  `/usr/bin/orca` on this box is the GNOME screen reader, not quantum chemistry.
- `decoherence_calculation_venus/` (43 GB) — separate git repo. Its `work/` tree
  is ~42.8 GB of regenerable scratch, but it has **57 untracked source files**
  (30 tests, 22 workflow scripts, 5 modules) and 18 modified tracked files.
  Commit that work before touching anything there.
- `qchem_validation/` (1.2 MB) — expensive correlated-method benchmarks the
  manuscript cites; irreplaceable without rerunning Q-Chem.
- `UoP_Surrey_Venus_Collaboration/` — separate repo, clean, nothing unpushed.

## Totals — final

| | |
|---|---|
| Crashed ORCA scratch | 251 GB |
| A duplicate bundle | 11.7 GB |
| B superseded MD lineage | 17.0 GB |
| C superseded decoherence QM | 14.4 GB |
| D superseded v1 site-energy families (1399 dirs) | 9.8 GB |
| E transients + `resp_cr2/` strip | 0.8 GB |
| **Total reclaimed** | **~305 GB** |

Tree: 429 GB -> **124 GB**. Disk: 189 GB free -> **494 GB free** (79% -> 44% used).

## Still outstanding

1. Fix or retire `manuscript/JPCB_tandem.tex` (still `J = 117.2`).
2. Commit the untracked round-2 sources and the rescued ensembles.
3. Optional: the loose top-level scripts could be grouped into
   `scripts/{build,coupling,validation,figures,render}/`, but ~100 are tracked
   and cross-referenced by `README.md`, `ORCA_ONLY_PIPELINE.md` and each other's
   imports, so it needs the reference updates done with it — not a blind `git mv`.
4. `decoherence_calculation_venus/work/` is ~42.8 GB of regenerable scratch, but
   that repo has 57 untracked source files. Commit before cleaning it.

# T2 campaign runbook (v6t2a / v6t2b / v6t2c)

Three 100 fs windows at 1 fs cadence on the corrected 511-aa construct with
full-system embedding. 300 TDDFT frames x 2 sites, **zero new STEOM jobs**.

## Why

T2 was the weakest number in the manuscript on 2026-08-31:

- the text quoted `T2* = 15.6 fs` while `make_frame_dynamics.py` built the
  figures at `60 fs` -- a 4x internal inconsistency;
- the measured values disagreed among themselves: cumulant `11.63 fs`, direct
  propagation `11.75 fs` with J but `64.83 fs` gap-only, and `21.8-24 fs` from
  single-site `dq.V` probes;
- **all of them came from `tc_tandem_nvt_100fs_50frames_20260723`** -- the
  truncated 495-aa construct with the residue-complete 12 A cutoff embedding,
  superseded on 2026-08-27 and 2026-08-28 respectively.

And a ~12 fs T2 measured on a 4 fs grid is three points. The old stride test
(1 -> 11.63, 2 -> 11.75, 4 -> 10.69 fs) looks reassuring but carries no
information: striding a series that is already under-sampled cannot detect
under-sampling. Only a finer cadence can, which is what this is.

## Design

| window | v5 seed frame | `\|Delta\|` cm-1 | percentile | J cm-1 | velocity seed |
|---|---:|---:|---:|---:|---:|
| `v6t2a` | 900  | 826.6 | 71% | 28.84 | 20260901 |
| `v6t2b` | 1400 | 58.6  | 4%  | 28.85 | 20260902 |
| `v6t2c` | 1900 | 514.2 | 51% | 27.67 | 20260903 |

Three separated windows rather than one contiguous 300 fs block, so T2 gets a
configurational error bar and "does T2 depend on the instantaneous detuning?"
becomes measurable. The seeds span near-degenerate to strongly localised at
near-constant J, so detuning is the only variable.

**v5's frame 2400 is deliberately unused**: its `|Delta| = 110.8` duplicates
1400's, and its `J = 19.23` is a -3.1 sigma coupling outlier that would confound
the comparison.

Seeding each window on a frame whose STEOM reference already exists is what
makes the campaign free of ORCA. All 8 jobs are on disk under `steom_v5/`.

## Running it

    systemd-run --user --unit=v6t2 --same-dir \
        /home/robson/PetaChem/tddft_ensemble_decoherence/run_v6_t2_campaign.sh

Never launch from an agent shell: the harness kills its whole process tree when
the shell exits (that is how `chain_attempt7.sh` died on 2026-08-27).

Rebuild the MD first if the workdirs are absent:

    bash tandem_v6_t2_campaign/run_md_window.sh 900  -1209516.140068 20260901
    bash tandem_v6_t2_campaign/run_md_window.sh 1400 -1209328.412976 20260902
    bash tandem_v6_t2_campaign/run_md_window.sh 1900 -1209262.681868 20260903
    python scripts/initialize_v6t2_manifests.py

Then analyse:

    python scripts/summarize_100fs50.py --window v6t2a     # and b, c
    python scripts/t2_campaign_summary.py

## What is different from the v5 supervisor

- **One worker per window, not `--manifest-suffix` slices.** Required, not
  cosmetic: `run_ensemble.py` threads `previous["A"]/["B"]` transition dipoles
  frame to frame for density-overlap state tracking (`:532,:542`), so splitting a
  window across workers breaks the chain. At 1 fs the geometries barely move, so
  tracking should be far cleaner than v5's 0.5 ps case.
- **No STEOM stream**, hence no `steom_disk_guard.sh`, no `steom_pause_site.sh`,
  no 250 GB scratch, no throttle and no deadlock escape. All of that machinery
  existed only to keep ORCA from filling the disk.
- `OMP_NUM_THREADS=8` rather than 4: nothing is competing for cores now.

## Traps this campaign had to clear

**`addSolvent` is not seeded.** `--seed` reaches only the Langevin integrator
(`run_nvt.py:2012`). Re-running `--add-solvent` would place different waters, and
the atom ordering would no longer match the DCD being seeded from. The windows
therefore build from the v5 run's OWN `solvated_protonated.pdb` with
`--skip-solvation`. Verified atom-for-atom identical.

**The PDB's box is stale.** DCD frame 0 is at 92.0938 A (NPT still relaxing),
production frames at 90.6101 A. `--initial-dcd` takes positions and cell from the
DCD together, which is the documented fix for the v3 lineage's stale-box trap.

**`classical_relaxed.pdb` is the post-NVT MINIMISED structure**, not the
trajectory end. Never seed from it.

**Manifest rows must be checksum-bound before running.** `run_ensemble.py:549`
refuses any frame whose prepared coordinates do not hash to the manifest value,
and an unstamped row fails as `malformed_coordinate_or_point_charge_input`.
`initialize_v6t2_manifests.py` stamps them; `register_nvt_trajectories.register()`
does the same job for the older windows but overwrites the notes column.

**All four `run_ensemble.py` membership sets must be updated together** --
site-B offset (`:199`), charge cache (`:206`), partner-CR2 retention (`:226`) and
the full-system-embedding branch (`:236`). v5camp was omitted from the retention
set at launch and ran 5 frames with the partner chromophore zeroed. Verified for
v6t2a frame 0: `mm_net_charge_e = 0.8492` with `partner_cr2_charges_retained`
true and the restored residue at -1.000 e.

## Controls

- Velocity reseed: see `tandem_v6_t2_campaign/VELOCITY_RESEED_CONTROL.md`. No
  equilibration transient; the window mean sits 0.03 sigma from the parent
  trajectory's production distribution.
- System identity: seeded potential energies reproduce the v5 log to 0.06-0.23
  kJ/mol out of 1.21e6. Enforced at launch by
  `--expect-initial-potential-kj`.

## Blind spot, to be stated not glossed

100 fs caps the resolvable lag at ~330 cm-1. The slow and static components come
from the v5 1 ns ensemble instead. **Nothing samples 0.1-500 ps.**

---

# Chained: the 500 fs low-frequency extension (v6lfa/b/c)

Armed 2026-09-01 as `systemd --user` unit `v6lf-chain`, which waits for the
`v6t2` unit to go inactive and then runs the extension automatically.
`chain_v6lf.sh` **refuses to start** if v6t2 ended at fewer than 300/300: an
aborted or retired run needs a human, and chaining 8 h of GPU onto a broken
result is worse than stopping.

## Why

The v6t2 windows resolve only to `1/(c*100fs)` = **334 cm-1**, which is coarser
than `kT` = 208 cm-1. The part of the spectral density entitled to a classical
treatment (`hbar*w < kT`) is therefore invisible to them. 500 fs gives
**67 cm-1** and crosses that boundary, so the fast/slow split can be *measured*
across the classical/quantum line instead of inferred.

This matters because the T2 the 1 fs campaign computes is dominated by a
~1350 cm-1 mode at `hbar*w/kT = 6.5` -- quantum-frozen, and properly vibronic
rather than a dephasing bath. See the memory note
[[t2-campaign-and-sampling-floor]].

## What it does, in order

1. Waits on `systemctl --user is-active v6t2.service`.
2. Aborts unless v6t2 reached 300/300.
3. Builds three 500 fs MD windows (`run_md_ext.sh`): 2000 steps x 0.25 fs,
   DCD every 16 steps = 4 fs = 125 frames, from the SAME seed frames, velocity
   seeds and integrator seeds as v6t2. Each still verifies its seeded potential
   energy against the v5 log.
4. `check_lf_overlap.py`: hashes DCD frame 24 of the 500 fs run (t = 100 fs)
   against frame 99 of the 100 fs run. Identical seeds and system should make the
   trajectories identical, differing only in save interval -- but that is checked,
   not assumed. MATCH means the 1 fs and 4 fs records splice into one continuous
   500 fs window per seed; DIFFER means each stands alone and only the >100 fs
   frames are used. Either way the campaign proceeds.
5. `initialize_v6lf_manifests.py`: 100 rows per window, DCD frames 25-124
   (t = 104-500 fs). The first 100 fs is deliberately NOT recomputed -- v6t2
   already has it at 1 fs.
6. `run_v6_lf_campaign.sh`: same supervisor as v6t2, one worker per window,
   no STEOM. 300 frames at the measured 99 s/frame effective = **~8.2 h**.

## Cost and disk

300 frames ~ 19 GB of frame output, on top of v6t2's ~19 GB. No STEOM scratch.

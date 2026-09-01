# v5 campaign -- what is running and how to resume it

Started 2026-08-30 19:25. Everything is resumable: work already done is skipped,
so any of these can be restarted at any time without losing or duplicating a
frame.

## The experiment

MD: `tc_tandem_v5_campaign/tandem_v5_1250ps_500fs.dcd`, 2 fs step, 300 K,
saved every 0.5 ps. DCD frames 0-399 are NPT and discarded, 400-899 are NVT
equilibration and discarded, **900-2899 are production**.

| stream | frames | per frame | settings |
|---|---|---|---|
| TDDFT | stride 4 -> **500** (900, 904, ... 2896) | 2 sites | def2-SVP/wB97X-D3, 7 roots, full-system embedding (75,035 charges) |
| STEOM | stride 500 -> **4** (900, 1400, 1900, 2400) | 2 sites | DLPNO-STEOM-CCSD/def2-SVP, same geometry and same field |

Both quantum-chemistry streams read files produced by the *same*
`prepare_production_frame.py` invocation, so the STEOM references sit on
geometries the TDDFT ensemble also contains. That is what makes the
TDDFT -> STEOM shift measurable rather than assumed.

## Units

| unit | what it does |
|---|---|
| `v5camp` | supervisor: 3 TeraChem workers + the STEOM stream, restarts and retries |
| `steomdisk` | resource throttle: caps STEOM at one running frame pair, by SIGSTOP |
| `steomcomp` | re-runs any STEOM job left incomplete, 2 at a time |

    systemctl --user status v5camp steomdisk steomcomp

## Resume after an interruption

    systemctl --user start v5camp        # workers skip completed rows
    systemctl --user start steomdisk steomcomp

Nothing else is needed. The three TDDFT workers hold **disjoint round-robin
slices** of the 500 rows in `manifests/v5camp_{a,b,c}_frames.csv`; they are
separate files because `save_manifest()` rewrites the whole csv after every
frame and concurrent workers would clobber each other. Round-robin rather than
three contiguous blocks, so a partial campaign is still a uniform sample of the
whole nanosecond rather than only its first third.

## Progress and results

    python scripts/analyse_v5_campaign.py

merges the three slices back into `manifests/v5camp_frames.csv`, reports the
coupling and detuning ensembles against v4full (J = 27.43 +/- 0.18, n = 99) and
the rigid-fit value the campaign replaces (32.82), and does the STEOM
calibration once matched points exist. Writes `results/v5camp_summary.json`.

## Measured, against the plan's estimates

- **Three-way GPU concurrency gains ~10%, once the box is not swapping.**
  Measured effective rates: 137 s/frame with one worker and no ORCA; 148 s/frame
  with three workers against four ORCA jobs *while thrashing*; **124 s/frame**
  with three workers once STEOM was capped at one frame pair (15 frames in
  31 min, per-worker times settled at 307-323 s). Projection for 500 frames:
  **~17 h**.

  The middle number is why the throttle matters, and why an early reading of it
  was misleading: three workers looked no better than one, which invited the
  conclusion that the GPU was the wall. It was swap. Do not benchmark
  concurrency on this box without first checking `pswpin` in `/proc/vmstat`.
- The GPU is nonetheless the per-frame cost: 18 Davidson iterations for 7 roots
  against 19 s of SCF, so cutting roots is the only large lever left -- and it is
  the state-tracking that needs those roots, so it is not a free one.
- **STEOM scratch is 27-38 GB per job**, four concurrent, peaking during the
  RI-PNO transformation before ORCA releases intermediates (it settles to
  ~24 GB in the CCSD iterations).
- **Memory, not disk, is the binding constraint.** Twenty ORCA ranks in CCSD
  alongside three TeraChem workers exhausted 7 GB of swap and drove sustained
  thrash -- 4-11 MB/s swap-in, 21% idle -- while disk still had 325 GB free.
  Capping STEOM at one running frame pair took swap-in from ~10,600 to
  ~530 pages/s.

  The throttle caps concurrency; it does not threshold on swap. A threshold
  version flapped: swap goes calm *because* a job was paused, that reads as
  permission to resume it, and it unpaused both pairs within four minutes and
  went straight back into thrash. Load-shedding has to be judged on the load,
  not on the symptom the shedding just relieved.

## Two things checked before trusting any of it

- `PARTNER_CHARGE_CONTROL.md` -- the partner chromophore was at 0.0 e instead of
  -1 e until this campaign; worth 36 cm^-1 rms on the detuning.
- `SELFINTERACTION_VERIFICATION.md` -- `pointcharges_self_interaction no` is a
  constant on the ground state and the better-converged convention.

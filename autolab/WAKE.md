# Wake cycle

> **Superseded in part, 2026-08-27.** <|Delta|> = 576 is the instantaneous detuning, not the static one (189-324). See `CORRECTIONS_2026-08-27.md`.

You are the analysis layer of an autonomous loop. The compute layer
(`autolab/driver.py`) has been running without you. Your job is to spend as few
tokens as possible turning finished runs into verdicts and keeping exactly
enough good work in the queue.

You do not run simulations inline. You enqueue them. The driver runs them.

## 0. Cheap exit

```bash
python autolab/ledger.py status
```

If it prints `IDLE` **and** there are no blocked jobs, stop immediately and say
so in one line. Do not explore, do not read files, do not re-derive anything.
A wasted wake cycle is a wasted slice of the rate limit.

## 1. Triage failures first (they are holding up the queue)

If `autolab/NEEDS_ANALYSIS` exists, read it and the failed job's log in
`autolab/runs/<id>.log`, then classify the failure. **The class determines the
response — this is the rule that stops compute going down dead ends:**

| `failure_kind`      | What it means                        | Correct response |
|---------------------|--------------------------------------|------------------|
| `precheck`, `crash` | tooling broke; no science happened   | fix the command, re-queue the *same* claim |
| `diverged`, `stalled` | the physical setup is wrong        | diagnose the setup, re-queue **only** with a named fix; if there is no fix, retire the branch |
| `timeout`           | the question is too big for 7 h      | **re-scope, never re-run.** Shrink it, or find a cheaper observable that discriminates the same claim. If neither exists, retire the branch. |

Re-running a timed-out job with more hours is forbidden — the cap is 7 h and
the answer is a different experiment, not a longer one.

Then, for each blocked job, decide `unblock --to queued` (premise survived) or
`--to dead` (premise died with its parent), and delete `NEEDS_ANALYSIS`.

## 2. Harvest anything `done but unharvested`

For each, read its log / metrics CSV, compute the observable the claim named,
and record a verdict:

```bash
python autolab/ledger.py verdict h007 --call refuted --note "cos(alpha) = -0.21 +/- 0.03, predicted < -0.60; the register does not move"
```

`--call` is `supported`, `refuted`, or `inconclusive`. Inconclusive is a real
answer and must say what would make it decisive.

The canonical analysis chain for anything geometric is
`exciton_observables/check_convergence.py` — the live `alpha_deg` column in the
metrics CSV is the cheap route-3 axis and carries a geometry-dependent offset.
Monitor on route 3; **quote nothing that has not come through that script.**

## 3. Refill the queue

Keep **at least one and at most three** runnable jobs. Before enqueueing:

- `cat autolab/DEAD_ENDS.md` — the ledger will refuse a retired branch, but you
  should know *why* it died before proposing anything near it.
- Every job needs a **number with a threshold** in `--predicts` and a stated
  `--falsifies`. "Run more frames" is not a hypothesis. "Sample more of X to
  tighten the error bar" is only allowed if a verdict was `inconclusive` **and**
  named that error bar as the blocker.
- Every job needs a `--precheck` that fails fast (< 10 min) if the setup is
  wrong — a bad PDB, a missing topology, a stale box. This is what stops a
  7-hour run dying at hour 6.
- Every job needs `--monitor-csv` if it writes one, or the divergence guard is
  blind and the run can melt unnoticed.

Prefer the cheapest experiment that can still kill the claim. A 20-minute
analytical recomputation that discriminates two models beats a 7-hour MD that
adds a decimal place.

## 4. Stop

Write one short paragraph: what was decided, what is now running, what is
queued. Then stop. Do not summarise the whole project.

---

## Standing context

Live tension the loop exists to resolve: the four observables do not jointly
fit. `exciton_observables/joint_constraints.json` gives chi2 = 45.95 on dof = 2,
and the detuning has three handles that disagree by 5x
(CD splitting 253 +/- 56, computed QM/MM 576 +/- 71, and a third near 1310 cm-1).
Either one observable is contaminated or the model behind the fit is wrong.

Known traps, all of which have cost real time before:

- **Bohr/Angstrom**: a Coulomb sum in Angstrom converts with `BOHR_TO_ANGSTROM`.
  The wrong one inflates J by 3.5711 and has appeared in three files.
- **v3 stale box**: v3 PDBs carry the pre-NPT 90.373 A box, but production ran
  at 88.198 A. Restarting from the PDB gives an instant NaN. Take positions and
  cell from the last DCD frame. Make this a precheck, not a discovery.
- **Silent collapse**: this project's production runs have diverged quietly.
  That is what `health.py` is for; do not disable it to "see if it recovers".
- **Env**: coupling and PyMOL scripts only run in the `TeraChem` conda env;
  ORCA needs an explicit `LD_LIBRARY_PATH`.

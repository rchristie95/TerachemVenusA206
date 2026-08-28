# Dead ends

Branches retired from the queue. `ledger.py add` refuses any `--branch` listed
here as a `## branch:` heading, so a later pass cannot re-spend compute on a
question that has already been answered in the negative.

Retire a branch with:

```bash
python autolab/ledger.py dead <branch> --why "..." --evidence "..."
```

Seeded below from results already established in this repo. Do not re-open one
of these without new evidence that specifically breaks the stated reason.

## branch: linker-tension

- retired: pre-existing result (overnight_linker_release/, commit 67c3461)
- why: releasing linker tension moves the interface angle the WRONG way, toward
  90 deg rather than the 131.3 deg the anisotropy requires. The 36 ns unbiased
  control drifts toward the crystal 106.6 deg, not 131.3, and the chirality
  triple product is negative in every frame.
- evidence: overnight_linker_release/control/control_metrics.csv,
  release/release_metrics.csv, analyze_morning.py

## branch: two-state-degenerate

- retired: pre-existing result (commit 67c3461)
- why: the dimer is detuning-dominated, |Delta| ~ 550-580 cm-1 >> 2|J| = 65.6,
  so the degenerate two-state model that produces alpha = 48.7 deg is the wrong
  branch of the solution. Superradiance forces cos(alpha) < 0.
- evidence: exciton_observables/superradiance_prediction.json

## branch: terachem-magnetic-dipoles

- retired: pre-existing result
- why: TeraChem does print magnetic dipoles, but they are unusable — the gauges
  disagree in sign (208/220) and R(len) moves 100% per 0.79 A of origin shift.
  PySCF reproduces it, so this is TDA plus range-separated exchange, not a code
  bug, and no amount of compute fixes it.

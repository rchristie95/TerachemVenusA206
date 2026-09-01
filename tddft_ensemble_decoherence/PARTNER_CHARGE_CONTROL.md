# The partner chromophore's charge in a full-system embedding

2026-08-30. Five matched frames, def2-SVP, full-system embedding, same
geometries; the only difference is what charge the *other* chromophore carries in
the field.

## The defect

`load_or_build_charges()` deletes CR2 by residue name and leaves 0.0 e behind.
`prepare_production_frame.py --retain-partner-cr2-charges` puts the partner's
AMBER RESP charges back, and every `--full-system-embedding` window is supposed
to pass it, because in those windows the field reaches the other site.

`v5camp` was not in that set in `run_ensemble.py`. Measured on frame 900, both
sites:

| | net field charge | charges pinned at exactly 0.0 |
|---|---:|---:|
| without `--retain-partner-cr2-charges` | +1.849200 e | 29 |
| with it | +0.849200 e | 0 |

The missing e sits ~24 A from the QM region -- the partner chromophore itself.

## Size of the error

New minus old, over frames 900/904/908/912/916:

| quantity | mean | rms | range |
|---|---:|---:|---:|
| site A energy | -10.67 | 50.86 | [-71.35, +69.60] |
| site B energy | -14.29 | 16.61 | [-22.49, -3.69] |
| detuning `Delta` | +3.63 | 46.20 | [-52.57, +73.30] |
| `J` screened | +1.72 | 3.01 | [+0.21, +6.66] |

All cm^-1.

The mean shift is small because the error is nearly symmetric between the two
sites; the rms is not, because the two chromophores are not equivalently placed
with respect to each other's field. **Against an ensemble `<|Delta|>` of ~621
cm^-1 this is a ~7% perturbation on the detuning, not a qualitative change.**
It is worth having right, and it is not what drives the detuning.

The single large `J` shift, +6.66 cm^-1 at frame 916, is a root reassignment,
not electrostatics: frame 916 has `|Delta| = 23 cm^-1`, the smallest of the five,
and it was the one frame flagged `ambiguous_overlap_margin` during calibration.
Excluding it, the residual shifts are **36.40 cm^-1 rms on the detuning and
0.53 cm^-1 rms on `J`**.

## Status

Fixed at `scripts/run_ensemble.py:198`; `v5camp` joins the retain set. The five
pre-fix frames are kept as this control under `frames/v5camp_nopartner/` and
`results/v5camp_nopartner_control.csv`, and the production ensemble was restarted
from zero rather than mixing the two embeddings.

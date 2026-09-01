# The QM/MM boundary treatment dominates the site-energy gap

2026-08-30. `run_boundary_test.sh`, 5 frames, three boundary policies with an
IDENTICAL 44-atom QM region. Baseline is the production window `v4full`.

| frame | baseline (link-only) | 1.8 A excl | shift | 3.0 A excl | shift |
|---|---:|---:|---:|---:|---:|
| 0 | 10.9 | 2029.1 | **+2018.2** | 1030.1 | **+1019.2** |
| 26 | 911.8 | 503.0 | -408.8 | 720.4 | -191.4 |
| 48 | 756.4 | ambiguous | - | 335.4 | -421.0 |
| 1005 | 1078.0 | 1420.8 | +342.8 | 830.5 | -247.5 |
| 1880 | 546.5 | ambiguous | - | ambiguous | - |

rms shift: **1205 cm^-1** (1.8 A), **573 cm^-1** (3.0 A).
For scale: ensemble `<|Delta|>` = 574-680 cm^-1, statistical s.e. +/-54-74 cm^-1.

## The prediction was wrong

Before running this I argued that boundary errors would cancel in `Delta`,
because the boundary is constructed identically at both sites (MM net charge
agrees to 13 decimals). **That reasoning is wrong.** A symmetric selection
*rule* does not give a symmetric *error*: the local environment at each site
differs frame by frame, so deleting the same class of charges perturbs E_A and
E_B by different amounts. The shifts are the same order as `Delta` itself, and
an order of magnitude larger than the statistical error bar we have been
quoting.

Frame 0 is the sharpest case. Its baseline gap of 10.9 cm^-1 -- the value
`EMBEDDING_FIX.md` built its alarm on, and the smallest of all 92 frames --
becomes 1030-2029 cm^-1 under a differently-drawn boundary. That near-degeneracy
is an artefact of the boundary policy, not a property of the frame.

A further symptom: 3 of 10 site calculations came back
`ambiguous_overlap_margin`, i.e. the boundary change perturbed the electronic
structure enough that the density-overlap tracker lost the bright state in
frames where the baseline was clean.

## What this does and does not establish

**Does not** say which policy is correct. None of the three is: all delete
charge without redistributing it, so a larger exclusion digs a larger
unphysical charge hole. `prepare_production_frame.py` has no
`--conserve-boundary-residue-charge` (that flag exists only in
`prepare_candidate2_frame.py` and was never back-ported).

**Does** establish that the spread across defensible boundary policies is a
systematic uncertainty on `Delta` of order 500-1200 cm^-1. The +/-54-74 cm^-1 we
quote is a statistical error on one policy and badly understates the real
uncertainty.

Caveat on the statistics: 5 frames, 3-4 usable per variant. The rms values are
noisy. The magnitude is not -- 1000-2000 cm^-1 shifts cannot be noise.

## What survives, and what does not

- **Survives**: the qualitative claim `|Delta| >> 2|J| = 65.6 cm^-1`. Every
  variant gives hundreds of cm^-1, so "detuning-dominated" is not at risk.
- **Survives**: the impossibility argument, which needs `Delta` only in the
  189-576 range at J = 32.82 (`notes/impossibility_note.tex` Table 1) -- the
  contradiction holds across that whole span.
- **Does NOT survive**: any claim resting on the precise value of `<|Delta|>`,
  including the agreement between the full-system 574 and the 12 A production
  576. That agreement now looks fortuitous rather than confirmatory: two
  policies happening to land close, inside a systematic spread of ~1000 cm^-1.

## Next step, which this reverses

I previously ranked boundary treatment LAST, behind polarisable embedding and QM
region size, on the cancellation argument. That ranking is withdrawn.
Back-porting `--conserve-boundary-residue-charge` from
`prepare_candidate2_frame.py` (~1 h) and repeating this test is now the priority:
charge conservation is the standard fix for exactly this artefact, and if it
collapses the spread then the production numbers can be trusted again.

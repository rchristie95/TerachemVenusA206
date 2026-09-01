# Figure data provenance audit, round 3

2026-08-30. Every figure in `JPCB_tandem_round_3.tex` traced to its generator
and its input data.

## Summary table

| figure(s) | generator | data | method | current? |
|---|---|---|---|---|
| Fig_Tandem_Separation, Fig_Tandem_Histogram (Fig. 5) | `make_tandem_panels.py --data-dir results/coupling_tandem_v5` | v5 campaign, 500 frames | per-frame transition densities, full-system embedding | **current** |
| Fig_Bloch_Grid, Fig_SSE_Site, Fig_SSE_Adiabatic, Fig_Purity_Coupling, Fig_ME_Site, Fig_ME_Adiabatic | `open_quantum_dynamics.py` | `J_OPT = 27.0571` | v5 ensemble mean | **current** |
| frame_deloc/*, frame_loc/* (12 panels) | `make_frame_panels.py` | 466 ps and 714 ps frames, each with its OWN recomputed J | v5 extremes | **current** |
| Fig_Detuning_Histogram | `make_frame_panels.py` | v5 manifests, 500 frames | full-system embedding | **current** |
| Fig_Bridge_Detuning | `bridge_detuning_analysis.py` | v5, 447 frames + `results_detuning_decomposition_v5.json` | per-frame TDDFT, full-system embedding | **current** |
| steom_transition_{dimer,qm}, steom_difference_{dimer,qm} | PyMOL renders | static density isosurfaces | no ensemble dependency | n/a |
| FigS1_qmregion_{a,b,c} | PyMOL renders | QM region illustration | no ensemble dependency | n/a |
| toc_steom | `label_toc.py` | composite of the above renders | n/a | n/a |

**All data-dependent figures were regenerated on 2026-08-31** against the v5
campaign. Three generator bugs were fixed in the process, recorded at the end of
this file.

## The measured size of the problem

The two coupling datasets can be compared directly:

| | rigid fit, static STEOM density, 1000 frames | per-frame recalculated TDDFT, 92 frames |
|---|---:|---:|
| `J` (screened) | **32.82 +/- 0.05** cm^-1 (sd 1.55) | **27.72 +/- 0.14** cm^-1 (sd 1.35) |
| inter-dipole angle | 100.78 deg | 99.64 deg |
| centroid separation | 24.69 A | 24.63 A |

Ratio new/old = **0.845**.

The geometry agrees closely (angle within 1.1 deg, separation within 0.06 A),
so the rigid placement of the density is not itself introducing a geometric
error. Most of the coupling difference is the transition-dipole magnitude:
STEOM gives 9.807 D and TDDFT 8.677 D, and `J ~ |mu|^2` predicts a ratio of
0.783. The residual, 0.845 / 0.783 = **1.079**, is the part attributable to
recalculating the density per frame rather than rigidly placing one.

So the outdated method is worth about **8%** in `J`, on top of a **22%** effect
from the choice of electronic-structure method. Both are smaller than the
factor-level errors this project has previously found, but 8% is not nothing and
the paper does not currently disclose it.

## Caveat on that comparison -- now RESOLVED by the v5 campaign

It was not a clean controlled test. The two datasets differed in three ways at
once: electronic-structure method (STEOM vs TDDFT), trajectory (the 1 ns
production run vs the 8 ps v4 high-cadence run on the corrected 507-residue
construct), and sample size (1000 vs 92). The 8% figure was an estimate, not a
measurement.

**The v5 campaign closes this** (2026-08-31, `tddft_ensemble_decoherence/`,
see `V5_CAMPAIGN_RUNBOOK.md`). 500 frames over a fresh 1 ns trajectory --
matching the rigid-fit dataset's sampling -- with per-frame transition densities
and full-system 75,035-charge embedding, 500/500 complete and zero failures:

| | rigid fit, static density, 1000 frames | **v5, per-frame, 500 frames** |
|---|---:|---:|
| `J` screened | 32.82 +/- 0.05 (sd 1.55) | **27.057 +/- 0.113** (sd 2.518) |
| `J` PDA screened | - | 25.019 +/- 0.107 |
| `<\|Delta\|>` | - | 608.6 +/- 20.7 (sd 462.6) |

Ratio new/old = **0.824**, i.e. the outdated method inflates `J` by **17.6%**,
not the 8% the confounded comparison suggested. The v5 value also lands within
0.4 cm^-1 of the independently recalculated v4full ensemble (27.43 +/- 0.18,
n = 99) computed on a *different* trajectory, which is the cross-check that the
earlier estimate could not offer.

Figure-ready data in the shape `make_tandem_panels.py` reads is emitted to
`tddft_ensemble_decoherence/results/coupling_tandem_v5/` by
`scripts/analyse_v5_campaign.py --emit-coupling-samples`.

## What this means for each figure

**Figure 5 panels and all the open-quantum-system panels** rest on
`J = 32.82 cm^-1` from the rigid-fit route. That value is the paper's stated
method and is documented in the Computational Methods section, so the figures
are not wrong as labelled. They are, however, the older of two available
treatments, and the newer one gives 27.72.

**Fig_Detuning_Histogram** uses the 12 A per-site embedding. That scheme is now
known to carry a boundary systematic of 500-1200 cm^-1 rms
(`tddft_ensemble_decoherence/BOUNDARY_TEST_RESULT.md`). The full-system
ensemble gives `<|Delta|>` = 574 cm^-1 against the 576 behind this histogram, so
the central value is unaffected, but the width and the tails of the plotted
distribution inherit the boundary sensitivity.

**Fig_Bridge_Detuning** is the only data figure built on the full-system,
per-frame-recalculated ensemble.

## Recommendation

Regenerating Figure 5 from `v4full` is not a like-for-like swap: it would change
the trajectory, the construct, the sample size and the electronic-structure
method simultaneously, and 92 frames over 8 ps is a weaker statistical basis for
a coupling distribution than 1000 frames over 1 ns. The defensible options, in
order of cost:

1. **Disclose.** State in the Methods that the coupling ensemble uses a rigidly
   placed STEOM density, that a per-frame recalculated TDDFT comparison on 92
   frames gives 27.72 cm^-1, and that the difference is dominated by the
   transition dipole with roughly 8% attributable to the placement. This is
   accurate, costs nothing, and is what a referee would want.
2. **Recompute the coupling per frame on the production trajectory** at the
   STEOM level, which isolates the placement effect properly. Expensive.
3. **Regenerate Figure 5 from v4full** and re-label it as the high-cadence
   corrected-construct ensemble. Changes what the figure claims.

Option 1 is recommended, with option 2 flagged as the calculation that would
settle it.


## Generator defects found while regenerating (2026-08-31)

**1. `make_frame_dynamics.py` dephased twice as fast as it claimed.** It set
`GAMMA = 2.0 / T2`, but for `L = sqrt(gamma/2) sigma_z` the dissipator is
`(gamma/2)(sz rho sz - rho)` and sigma_z flips the sign of an off-diagonal, so
`d(rho_AB)/dt = -gamma rho_AB` and the rate is `1/T2`. Numerically confirmed:
the constant delivered a 30 fs coherence time while the docstring said 60 fs and
the adjacent panels from `open_quantum_dynamics.py` (`gamma = 1/t2_star`) used
60 fs. Fixed.

**2. The bridge analysis correlated across trajectories by frame number.** Frame
ids are per-trajectory DCD indices. `results_detuning_decomposition.json` holds
v4full ids (0-1999) and v5camp uses 900-2896; they share **8 integers that
correspond to unrelated configurations**, and correlating on them returned
`r = -0.465` for the electrostatic model where the properly matched calculation
gives `+0.203`. The script now refuses any match covering less than half the
frames, and a v5 decomposition was computed.

**3. `decompose_detuning.py` was pinned to the v4 topology and rejected the
partner chromophore.** It loaded the v4 PDB (74,857 atoms) regardless of window,
which silently misaligns the MM index map against v5 (75,079). Its charge guard
also failed outright on v5 because `--retain-partner-cr2-charges` restores CR2
values the cache holds as 0.0; it now verifies that the difference is exactly the
29 CR2 atoms summing to -1.0000 e rather than either failing or ignoring it.

**4. The Fig. S1 caption contradicted itself**, giving both "44-atom QM region"
and "a count of 45 atoms". The geometry files hold 44, and the caption's own
stoichiometry C19H18N3O4 sums to 44. Corrected.

---

## 2026-09-01 — T2 reframing (v6t2 campaign)

**Campaign**: `v6t2a/b/c`, three 100 fs windows at 1 fs on the corrected 511-aa
construct with full-system embedding, seeded on v5 frames 900/1400/1900 to span
the detuning at near-constant J. 300/300 frames, zero failures, 17h15m. Zero new
STEOM: all 8 references reused from the v5 sweep.

**Figures regenerated** at `T2* = 25` fs (was 0.060 ps):

| generator | constant | figures |
|---|---|---|
| `open_quantum_dynamics.py:63` | `T2_STAR = 0.025` | `Fig_Bloch_Grid`, `Fig_ME_*`, `Fig_SSE_*`, `Fig_Purity_Coupling`, `Fig_T2_Sweep` |
| `manuscript/make_frame_dynamics.py:34` | `T2 = 0.025` | `Fig_FrameDynamics` |

25 fs is the mean of the two windows that reach `1/e` (19.7 and 30.7 fs). The
near-degenerate window does not reach it within 100 fs at all, so this is a
representative value for an illustrative figure, **not** an ensemble constant;
the text carries the 20--31 fs range and says so.

The previous 0.060 ps was documented at `open_quantum_dynamics.py:13` as
"borrowed from photosynthetic systems" -- it was never measured for this system,
and the text meanwhile quoted a computed 15.6 fs, so figures and text disagreed
by 4x. Both now come from the same measurement.

**Text changed at `:179` and `:660`.** The claim that the doublet is
spectroscopically unresolved no longer rests on `1/(pi c T2*)` being a
homogeneous linewidth. It rests on the CD couplet width, 343 cm^-1 computed with
no free parameters against 366 cm^-1 measured, ~6x the corrected
`2|J| = 54.1 cm^-1` -- which holds however that width divides between vibronic,
inhomogeneous and homogeneous parts.

**Why the reframing.** The gap modulation is dominated by a ~1350 cm^-1 mode at
`hbar*w/kT = 6.5`, near its quantum ground state, which a classical trajectory
represents as a stochastic bath when it actually produces a Franck--Condon
progression. And every mode thermally populated enough to be classical
(`hbar*w < kT`) has period > 160 fs, quasi-static over the coherence decay. There
is no window in this system where a classical bath yields homogeneous dephasing,
so `1/(pi c T2*)` is no longer called a homogeneous linewidth.

**New supporting numbers, all traced to data files:**

| quantity | value | source |
|---|---|---|
| `T2*` per window | 19.7 / none / 30.7 fs | `results/t2_campaign_summary.json` |
| reorganisation energy, per site | 412 cm^-1 | `results/gap_spectral_density.json` |
| Huang--Rhys, upper bound at 1370 cm^-1 | `S <= 0.301` | same |
| independent `S` from electronic structure | 0.253 | `terachem_site_energy_cd`, unchanged |
| `sigma_site_slow`, same-trajectory | 337 cm^-1 | `results/fast_slow_disorder_split.json` |
| population delocalisation, near-degenerate window | 0.47 -> 0.20 over 100 fs | `results/delocalisation_lifetime.json` |

`sigma_site_slow = 337` is NOT a replacement for the model's 159, and the 500 fs
windows have now settled why. The two numbers are quoted at different fast/slow
boundaries -- ours at 100 fs, the model's at 592 fs -- and a shorter boundary
leaves more variance on the slow side by construction, so they were never
comparable.

Measured at a 396 fs boundary on single trajectories (`v6lfa/b/c`), the
within-window sd reaches 528-540 cm^-1 per site against a 1 ns total of 529-558,
i.e. **~97% of the site-energy disorder is realised within 400 fs**. The slow
remainder, taken directly as the between-window spread of the window means rather
than as a quadrature difference, is 72 cm^-1 (site A) to 184 cm^-1 (site B) --
**bracketing the model's 159 and consistent with it**.

Two cautions on that. With three windows the between-window estimate carries only
2 degrees of freedom, so its 95% interval spans roughly 0.52x-6.3x: this bounds
the slow disorder rather than measuring it. And the quadrature route
(`sigma_slow^2 = sigma_total^2 - sigma_fast^2`) is ill-conditioned once the fast
part captures 97%, differencing two numbers near 295,000 that each carry ~12%
sampling error; by that route the slow variance is not significantly different
from zero.

**The CD forward model therefore does not need refitting.** `SIGMA_SITE_SLOW =
159.0` survives the same-trajectory check; its own caveat that the value is
approximate stands, but it is not wrong, and the computed 343 cm^-1 against the
observed 366 cm^-1 is unaffected.

Compile after these edits: 28 pages, 0 undefined references.

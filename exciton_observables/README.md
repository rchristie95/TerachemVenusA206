# Superradiance as an independent test of the dimer geometry

Round-2 status, 2026-08-12. This directory makes contact with an experimental
observable the manuscript quotes but never inverts against our own geometry:
the **57 ± 4 ps** faster fluorescence decay of dVenus-TD relative to the
single-chromophore dVenus-TDX control, on a 3.026 ns lifetime (Nguyen et al.
2025, Fig. 7a).

## The relation, with the two corrections that matter

For the two-site Hamiltonian the eigenstate dipole strengths are
`|mu_pm|^2 = mu^2 (1 +/- sin(2 theta) cos alpha)` with `sin 2 theta = 2J/Omega`
and `Omega = sqrt(Delta^2 + 4 J^2)`. The manuscript's Eq. (superradiance) uses
this at the ensemble-mean detuning and without a thermal factor. Both matter:

1. **Thermal population.** Emission comes from a thermally populated mixture of
   both exciton states, and the upper one is *sub*radiant, cancelling part of
   the lower state's enhancement. The correct expression is

       Dtau/tau = Phi * tanh(Omega / 2kT) * (2J/Omega) * (-cos alpha)

   At Omega ~ 300-600 cm^-1 against kT = 208.5 cm^-1 that factor is 0.6-0.9.

2. **Ensemble vs mean.** `sin 2 theta` is strongly convex in Delta and the
   computed detuning distribution is broad (SD ~700 cm^-1), so the mean of the
   observable is not the observable at the mean. `predict_superradiance.py`
   averages frame by frame, with J, Delta and cos(alpha) taken from the same
   frame and the same transition densities so their relative signs are physical.

## Result 1 — the sign fixes the branch, and Nguyen picked the wrong one

A *faster* lifetime requires `cos alpha < 0`. The anisotropy analysis is even in
cos(alpha): |cos alpha| = 0.660 admits **48.7 deg or 131.3 deg**. Nguyen quote
48.7 deg. With their own J and detuning that predicts **-106 ps — subradiance, a
slower lifetime**, contradicting the measurement on the same page.

So the experimental angle is **131.3 deg, obtuse** — the same side as our
computed 96-110 deg. The structural model was never qualitatively wrong.

## Result 2 — a hard floor no detuning can rescue

`tanh(Omega/2kT) * (2J/Omega)` is maximised as Delta -> 0 at `tanh(2J/2kT)`
= 0.156 for J = 32.8 cm^-1, so reproducing `(1/Phi)(Dtau/tau) = 0.0330`
requires **|cos alpha| >= 0.212 at any detuning whatsoever**. This is a
geometry constraint, not an electronic-structure one.

## Result 3 — the crystal register lies on the superradiance constraint curve

**Read the caveat under this table before quoting it.** Superradiance is ONE
equation in TWO unknowns, `(|cos alpha|, |Delta|)`. It defines a curve, not a
point, and the crystal register lying on that curve is a consistency check the
crystal *passes* — not a unique determination of the geometry.

Rigid STEOM density placed by Kabsch fit on the 19 shared CR2 heavy atoms:

| geometry | alpha (deg) | cos alpha | J (cm^-1) | predicted Dtau (ps) |
|---|---:|---:|---:|---:|
| release arm, 36.6 ns (perturbed) | 90.2 | -0.003 | 28.70 | 0.7 |
| production ensemble, 1 ns | 100.8 | -0.187 | 32.82 | 18.5 |
| control arm, 36.6 ns | 102.3 | -0.213 | 29.51 | 33.5 |
| **1MYW biological dimer (crystal)** | **110.4** | **-0.349** | **30.82** | **57.0** |
| *measured* | | | | **57 ± 4** |

The crystal register matches at |Delta| = **570 cm^-1**, against an
independently computed QM/MM detuning of **549 / 576 / 581 cm^-1** across three
ensembles.

**Caveat (added after the fact — the first version of this file overstated
it).** Other points on the same curve fit equally well, and one of them is the
experimental angle:

| \|Delta\| (cm^-1) | needs \|cos alpha\| | alpha (deg) | |
|---:|---:|---:|---|
| 0 | 0.212 | 102.2 | unconditional floor |
| 253 | 0.237 | 103.7 | what the CD splitting profile requires |
| 576 | 0.330 | 109.3 | our computed QM/MM detuning; crystal is 0.349 |
| 1310 | 0.663 | 131.5 | where the anisotropy's \|cos alpha\| = 0.660 lands |

So the crystal geometry is *allowed* by superradiance at our computed detuning,
and the MD ensemble (|cos alpha| = 0.213, needing 0.330) is *excluded* by it.
But (131.3 deg, 1310 cm^-1) satisfies superradiance AND the anisotropy
simultaneously, and our detuning calculation contradicts that pair at ~10 sigma.
The tension is real and unresolved; see "What is actually still open" below.

## Result 3b — the escape route is closed: transfer really is ultrafast

The obvious way to rescue our geometry would be for the energy transfer to be
incomplete within the 26.5 ps IRF, which would let a smaller |cos alpha| produce
the observed R0 drop. A generalized-Forster/Marcus rate from our own J and
detuning kills that: with lambda = 355-508 cm^-1 the transfer time is
**0.9-1.7 ps**, 30-70x faster than the IRF, at every detuning considered.
Nguyen's complete-transfer assumption is sound and their |cos alpha| = 0.660
cannot be explained away.

## Result 4 — it is the MD relaxation that breaks the agreement

Both force fields relax the interface 10-14 deg away from the crystal, and that
relaxation is exactly what destroys agreement with experiment. The unbiased
control arm drifts *back* toward the crystal over its own trajectory
(99.3 -> 104.8 deg across 36.6 ns) and had not converged. A long continuation
is running in `../overnight_linker_release/control_long/` to see whether alpha
settles at ~110 deg and the prediction at 57 ps.

## Result 5 — alternative lattice registers are excluded

`notes/lattice_dimer_scan.py` offers two obtuse packing alternatives, and
`build_alt_register.py` materialises them. Both are killed by their couplings:

| register | alpha (deg) | separation (A) | J (cm^-1) | needs \|cos alpha\| |
|---|---:|---:|---:|---:|
| op5, lattice (1,0,-1) | 122.6 | 35.9 | -1.10 | 6.26 |
| op6, lattice (0,2,-1) | 169.6 | 44.0 | 6.34 | 1.09 |

Both require |cos alpha| > 1, i.e. impossible. Despite op5's attractive
26.5 A linker span (vs 54.2 A for the biological dimer), the biological dimer
is the only viable register.

## Files

- `predict_superradiance.py` — per-frame predictor over the QM/MM site-energy
  ensembles; writes `superradiance_prediction.json`
- `build_alt_register.py` — materialises lattice registers from 1MYW
- `alt_register_op5.pdb`, `alt_register_op6.pdb` — the two candidates
- `transforms_{control,release}.npz` — per-frame rigid CR2 fits for the 36.6 ns
  ensembles (from `extract_cr2_transforms.py`)
- `geom_{control,release}.npz` — per-frame cos alpha, angle, J, separation

## What is actually still open

**The detuning is the weakest link, and three handles on it disagree by 5x:**

| source | \|Delta\| (cm^-1) |
|---|---:|
| CD splitting profile inversion | 253 |
| QM/MM CAM-B3LYP, n = 95 | 576 ± 71 |
| anisotropy + superradiance combined | 1310 |

Whichever of these is right decides between alpha ~ 110 deg (crystal) and
alpha ~ 131 deg (experimental). The 5x spread is a *systematic* disagreement,
not a statistical one — the QM/MM standard error is only 71 cm^-1, so more
frames cannot bridge it.

Two arguments make the computed 576 harder to dismiss than it looks: the QM
region already contains Tyr202 (the omission that
`steom-vs-eomccsd-validation` showed breaks site energies), and the detuning is
a *difference* between two chemically identical chromophores, so functional and
basis-set errors largely cancel. That points the finger away from the
electronic structure and toward either the structural model or the
two-state framework being applied to all three observables at once.

## Compute triage (2026-08-12) — what is NOT worth GPU time, and why

Each of these was checked before committing time to it, and none pay:

- **Regenerating the cap-masked density for true TDC couplings.** Across the
  1000-frame production ensemble the TDC/PDA ratio is 1.1872 ± 0.0114 and the
  sign agrees in 100% of frames, so the scaled point-dipole coupling is good to
  ~1%, i.e. ~0.6 ps on a 57 ps observable. (Only 4 of the 7 NTO pair cubes
  survive anyway.) The caveat in the first version of this file was overstated.
- **More MD.** alpha has an integrated autocorrelation time of 2.93 ns, so the
  67.6 ns in hand already gives SEM ~1.2 deg; halving that again costs ~34 h and
  the qualitative verdict (the MD ensemble sits ~8 deg below the crystal) is
  already at ~5 sigma.
- **More QM/MM detuning frames.** SEM is 71 cm^-1 against a 5x systematic
  spread. Statistics are not the problem.

## Remaining caveats

- Phi = 0.57 and tau = 3.026 ns are taken from the manuscript's own reading of
  Nguyen; the prediction scales linearly in Phi.
- All three observables are interpreted in the same two-state
  (one-exciton, two-site) framework. Nguyen themselves hedge that a model with
  "only two excited states might be overly simplistic".

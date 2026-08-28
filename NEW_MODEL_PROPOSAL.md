# A four-observable model with no fitted parameters

2026-08-27. Constraint: every parameter comes from simulation. Nothing is
fitted to the spectra it is meant to predict.

## What is actually broken

Not the numbers -- three specific things the current two-state Holstein model
cannot do:

1. **The lobe asymmetry.** We compute 1.30, observed 2.34. A conservative
   couplet is a *theorem* of the two-state model, so the asymmetry has to come
   from a term the model does not contain. Vibronic borrowing with the computed
   S = 0.253 gets to 1.30-1.51 and no further.
2. **The transition dipole.** We use 9.81 D; the length and velocity gauges
   disagree by 1.5x and diverge with better basis, and three experimental lines
   want 7.5-7.9 D.
3. **The anisotropy/superradiance impossibility.** `cos^2 x <= |cos x|`, the
   anisotropy needs `<cos^2> = 0.436` and superradiance allows `<|cos|> <=
   0.287`. This is distribution-free: **no choice of geometry or parameters
   escapes it.** That is a statement that the physics is wrong, not the numbers.

(3) is the one that forces a new model rather than better parameters.

## The crack in (3)

The two experiments are treated as constraining the same angle. They do not.

- **Superradiance** uses the angle between the two chromophores' *transition*
  dipoles -- absorption and emission share the operator.
- **Depolarisation** compares the dipole that *absorbed* on one chromophore
  with the dipole that *emitted* on the other.

If a chromophore's emission dipole is rotated by `beta` from its absorption
dipole through excited-state relaxation, the anisotropy angle is `a + beta` and
the inequality no longer binds `a`. `beta` is a single-chromophore property and
is computable.

| alpha | beta needed so a+beta = 131.3 | superradiance satisfied? |
|---|---|---|
| 100.8 (1 ns MD) | 30.5 | no, alpha too small |
| **110.4 (crystal)** | **20.9** | yes |
| **111.2 (converged MD)** | **20.1** | yes |
| 119-120 (AlphaFold3) | 11.3-12.3 | yes |
| 131.3 (anisotropy branch) | 0 | yes, but over-predicts 10 sigma |

So `beta ~ 20 deg` reconciles the crystal / converged-MD geometry with **both**
experiments simultaneously. Running now: CAM-B3LYP/6-31G* S1 optimisation of
the 41-atom chromophore, same level as the S and omega already in use
(`emission_dipole/`).

**This is a prediction, not a fix.** If beta comes back at 5 deg the
impossibility stands and something else is wrong -- most likely the anisotropy's
complete-transfer assumption, or a state outside the two-level manifold.

## Rejected on the way (recorded so it is not re-tried)

The 57 ps lifetime drop as a **refractive-index effect**. The arithmetic works
-- `dn = 0.0088` reproduces it, against `dn ~ 0.01` already published for eGFP
tandems. But the control kills it: TDX deletes Tyr66, *not the barrel*, so TD
and TDX have identical protein mass and the same local n. `dn ~ 0` between
them. The chromophore's own polarizability gives `dn ~ 5e-4`, 17x too small.

## The model

One Hamiltonian, one ensemble, four observables, nothing fitted.

| parameter | source | status |
|---|---|---|
| `E_A(t), E_B(t)` | QM/MM TDDFT, bright-moment estimator | **have** |
| spectral density `C(w)` | full gap ACF, not one effective HWHM | **have data, needs analysis** |
| `J` | transition-density coupling | have (dipole caveat) |
| `S, omega_vib` | computed Hessian + ES gradient | **have** (0.253, 1369.7) |
| `alpha` distribution | marginalised over register, not one value | **new** |
| **`beta`** | **S1 relaxation, abs vs emission dipole** | **running** |
| monomer rotational strength | QM/MM, chiral protein field | **new, hardest** |

Two changes of principle:

**Marginalise the register instead of choosing one.** Today's trunk probe
showed alpha spans 55 deg across AlphaFold variants (110 crystal / 119 AF3 /
128 AF2-thinned / 165 AF2-full) with each method internally reproducible to a
few degrees. No method determines it. A model that picks one value is claiming
knowledge nobody has; the honest object is a distribution, and every observable
should be reported as a band over it.

**Use the full spectral density, not one effective linewidth.** The current
fast/slow split at 592 fs is a defensible cut but still a cut, and the CD width
result depends on where it falls. `C(w)` from the gap ACF removes the choice:
the cumulant handles every timescale correctly by construction.

## What each fix buys

| discrepancy | addressed by | expected |
|---|---|---|
| lobe asymmetry 1.30 vs 2.34 | monomer rotational strength | the missing non-conservative term; the difference spectrum retains ONE monomer's intrinsic CD, which is not conservative |
| dipole 9.81 vs 7.5-7.9 D | converged transition moment | gauge test says current value is not converged |
| anisotropy/superradiance | `beta` | breaks a distribution-free contradiction |
| couplet width | already closed | 342.7 vs 366.4, no free parameters |

## Order of work

1. **`beta`** -- running, cheapest, and decides whether (3) is solvable at all.
2. **Spectral density** -- data in hand, removes an arbitrary choice from the
   one result that currently works.
3. **Dipole** -- gates J, the CD amplitude and the lifetime; needs the TZVPD
   run that died, or a decision to quote TDDFT.
4. **Monomer rotational strength** -- hardest; the magnetic-dipole gauge
   problem is documented and unresolved, so this may not be computable at the
   accuracy needed.

## The discipline

No observable gets its own parameter. Absorption, CD, anisotropy and radiative
rate all come from the same Hamiltonian and the same ensemble. Where the model
misses, that is reported as a residual -- not absorbed by tuning a parameter
that only that observable constrains.

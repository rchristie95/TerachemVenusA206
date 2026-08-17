# Screening reconciliation: this work vs Cusick et al. 2026

All couplings in cm^-1. `epsilon = 1.77` is the optical dielectric used
throughout this repo; Cusick apply no screening at all.

## The comparison in both conventions

| Quantity | Structure | R (A) | Screened (eps=1.77) | Unscreened |
|---|---|---|---|---|
| This work, TDC | tandem MD, n=1000 | 24.69 +/- 0.32 | **32.82 +/- 1.55** | 58.09 |
| This work, PDA | tandem MD, n=1000 | 24.69 +/- 0.32 | 27.64 +/- 1.28 | 48.93 |
| Cusick, dVenus-vdW (crystal contact) | point-dipole from structure | 25.4 | 15.8-21.5 | **28.0-38.0** |
| Cusick, dVenus-TD (tandem, linkered) | from the Omega two-photon analysis | 26.0 | 18.1-22.6 | **32.0-40.0** |
| Cusick Table S2 (1MYW) | -- | -- | 15.3-24.3 | 27.0-43.0 |

Bold marks the number each group actually reports. Reading across the bold
cells compares a screened value against an unscreened one -- that is the
entire content of the apparent agreement.

## Like for like

- Screened throughout: ours 32.8 vs theirs 18.1-22.6 -- we are HIGH by ~1.61x.
- Unscreened throughout: ours 58.1 vs theirs 32.0-40.0 -- we are HIGH by ~1.61x.
- The discrepancy is the same in either convention, as it must be: the
  screening factor cancels from the ratio. What does not survive is the
  claim of a ~1% coincidence, which existed only because the two numbers
  were quoted in different conventions.

## The point-dipole partner

- Genuine partner of 32.82: **27.64 cm^-1** (same 1000 frames, same density, same convention), ratio **1.19**.
- NOT 13.31, which pairs with the superseded 74.38 on venus_dimer.pdb (crystal), a TDDFT S2 cube, 1 frame. reproduce_paper.py:115 records that TDC value as corrected to 20.83.
- So the near-field enhancement is 1.19x, not 5.6x. At R ~ 24.7 A the point-dipole approximation is accurate to ~19%, which makes cheap point-dipole
  estimates adequate for fluorescent-protein dimers at biological separations.

## What survives untouched

- sigma/<J> = 4.7% over 1000 frames at 293 K. This is a
  ratio, so it is invariant under both the units correction and the choice
  of epsilon. The coupling is remarkably rigid against thermal motion.

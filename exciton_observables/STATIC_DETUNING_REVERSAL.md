# Superradiance at the static detuning: the crystal match breaks, and J = 21 rescues it

2026-08-27. `predict_superradiance.py --delta-sd-cm`, run frame by frame on
`ens_v2_all` (n = 95). Follows from the fast/slow decomposition in
`tddft_ensemble_decoherence/audit/REPRODUCIBILITY_AND_TIMESCALE.md`.

## Why the detuning input had to change

The site-energy gap has tau_c = **8.5 fs** against hbar/2J = **81 fs**. The
fluctuation is ~10x faster than the exciton can respond, so the mixing angle is
set by the 81 fs-averaged detuning (sd 225-406) and not the instantaneous one
(sd 690). This is the same argument that fixed the CD couplet width
(`terachem_site_energy_cd/results/CD_WIDTH_RESOLVED.md`); it cannot be applied
to one observable and withheld from the other.

## Result 1 -- the published crystal match does not survive it

Predicted Dtau (ps), measured 57 +/- 4, frame by frame with fixed (alpha, J):

| geometry | alpha | sd 690 (instant) | sd 406 | sd 300 | sd 225 |
|---|---|---|---|---|---|
| production MD, 1 ns | 100.8 | 34.2 (-5.7s) | 41.2 (-3.9s) | 44.3 (-3.2s) | 46.5 (-2.6s) |
| control arm, 36.6 ns | 102.3 | 35.0 (-5.5s) | 42.2 (-3.7s) | 45.4 (-2.9s) | 47.6 (-2.3s) |
| 1MYW crystal dimer | 110.4 | **59.8 (+0.7s)** | 72.0 (+3.8s) | 77.5 (+5.1s) | 81.3 (+6.1s) |

The crystal register's headline agreement is specific to the instantaneous
detuning. At the static value it overshoots by 3.8-6.1 sigma.

## Result 2 -- the MD IS converged, and it converges ONTO the crystal

`check_convergence.py` on the 67.7 ns `control_long` continuation: alpha
settles at **111.19 deg** over the last half (J 29.8-30.3, separation 25.2 A),
against the crystal's 110.44. The question left open in README Result 4 is
answered -- **the 1 ns production ensemble at 100.8 deg was simply not
converged**, and the force field does reproduce the crystal interface.

(Careful: `control_long_metrics.csv` reports alpha = 107.6 +/- 4.2 deg over the
same window. That column uses the MD's own CR2 long axes; the superradiance
table uses the rigid STEOM density Kabsch-fitted on the 19 shared CR2 heavy
atoms, which is 111.19. Do not mix the two.)

So geometry is no longer the free parameter. Both plausible structures sit at
110-111 deg, and both overshoot at the static detuning.

## Result 3 -- the coupling is the free parameter, and it wants J = 21

At the converged geometry (alpha 111.19, cos alpha -0.3615):

| J (cm^-1) | source | sd 690 | sd 406 | sd 300 | sd 225 |
|---|---|---|---|---|---|
| 32.82 | STEOM dipole, 9.80 D | 66.0 (+2.2s) | 79.5 (+5.6s) | 85.5 (+7.1s) | 89.7 (+8.2s) |
| 30.00 | control_long TDC | 60.4 (+0.8s) | 72.7 (+3.9s) | 78.2 (+5.3s) | 82.1 (+6.3s) |
| 26.00 | intermediate | 52.4 (-1.2s) | 63.1 (+1.5s) | 67.9 (+2.7s) | 71.3 (+3.6s) |
| **21.00** | **TDDFT dipole, 7.84 D** | 42.4 (-3.7s) | **51.1 (-1.5s)** | **54.9 (-0.5s)** | **57.7 (+0.2s)** |

J scales as `|mu|^2`, so 32.82 * (7.84/9.80)^2 = 21.0.

**With the smaller dipole, the converged geometry reproduces the measurement
across the whole static range: -1.5 to +0.2 sigma.** With the STEOM dipole it
is excluded at 5.6-8.2 sigma.

## Three independent observables now agree the STEOM dipole is too large

| observable | says |
|---|---|
| Cusick 2026's own numbers | `|mu|` = 7.6-7.9 D, not 9.80 |
| absolute CD scale (`M ~ |mu|^4`) | epsilon implies ~7.5 D against our 9.80 |
| **superradiance at the static detuning** | **J = 21-26, i.e. `|mu|` = 7.8-8.7 D** |

A plausible mechanism already exists: the EOM-CC left/right eigenvector
asymmetry behind the TrEsp spectroscopic normalisation
([[steom-coupling-normalization-and-thermal]]). Note the STEOM dipole was shown
basis-converged (-0.54% for diffuse functions, against TDDFT's -10%), so this
is a normalisation question, not a basis-set one.

## What this does and does not settle

- Settles: the MD converges to the crystal interface; geometry is not the
  problem, and "the MD relaxation breaks the agreement" (README Result 4) is
  **withdrawn** -- it was an artifact of a 1 ns unconverged ensemble.
- Settles: the instantaneous/crystal pairing that reproduced 57 ps was two
  compensating errors -- a detuning ~3x too large against a coupling ~1.6x too
  large.
- Does NOT settle: which dipole is right. Three observables point one way and
  the basis-convergence argument points the other. **This is now the single
  highest-value open question in the project.**
- Unchanged: the anisotropy remains irreconcilable. At the static detuning
  superradiance needs `|cos alpha|` ~ 0.287, while the anisotropy needs
  `<cos^2 a>` = 0.436 and `cos^2 a <= |cos a|`. The gap widens from
  0.330-vs-0.436 to 0.287-vs-0.436, so
  [[anisotropy-superradiance-impossibility]] holds a fortiori.

## Reproduce

```bash
python3 exciton_observables/predict_superradiance.py \
  --ensembles terachem_site_energy_cd/results/ensembles/ens_v2_all.npz \
  --delta-sd-cm 300 --geometry-table
python3 exciton_observables/check_convergence.py --blocks 6
```

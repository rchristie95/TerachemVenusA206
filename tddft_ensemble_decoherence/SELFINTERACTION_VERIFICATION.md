# `pointcharges_self_interaction no` -- verification

2026-08-30. The v5 campaign runs with TeraChem's MM-MM Coulomb self-energy
switched off, which saves an O(N^2) sum over 75,035 charges every job. The term
contains no electron density, so it cannot enter the Fock matrix and cannot move
an excitation energy. This is the check that TeraChem agrees.

Frame 900 site A, def2-SVP/wB97X-D3, 7 roots. Within each pair the point-charge
file is byte-identical; only the flag differs.

## The term is exactly a constant

| field | total energy ON - OFF | spread of that shift across the 7 roots |
|---|---:|---:|
| 2,500 charges | -73.03559018 au | 2e-08 au |
| 75,035 charges | -7591.49681373 au | 1e-08 au |

The MM-MM sum lands identically on every state, as it must.

## But the excitation energies are not identical

Max |dE| between ON and OFF: **0.02 cm^-1** at 2,500 charges, **8.96 cm^-1** at
75,035.

A control settles what that is. Rerunning the *identical* flag-off input a
second time reproduces all seven roots to **0.000e+00 cm^-1** -- TeraChem is
bit-reproducible here. So the 8.96 cm^-1 is caused by the flag, and since the
flag cannot reach the Fock matrix, it reaches the answer through the **SCF
convergence path**: the convergence test is applied to a total energy 7.3x
larger when the constant is included, so the two runs stop at slightly
different densities, and CIS amplifies that into the roots.

## Consequence for the campaign

The choice is a convention, applied identically to every frame and both sites,
and the resulting dataset is bit-reproducible. It is also the *better-converged*
of the two conventions: with the constant excluded the total energy is 7.3x
smaller, so the same absolute convergence threshold is a tighter relative one.

The cross-convention difference, <= 9 cm^-1, is ~1.5% of the ensemble
`<|Delta|>` of ~621 cm^-1 and about a quarter of the 36 cm^-1 rms that
restoring the partner chromophore's charge was worth
(`PARTNER_CHARGE_CONTROL.md`).

Raw outputs under `selfint_check2/`.

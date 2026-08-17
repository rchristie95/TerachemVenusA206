# How much of the delta gap is the axis definition?

Cusick estimate delta from a **structural proxy** — the vector joining OH and
CB of the tyrosine chromophore precursor — because AlphaFold3 "does not
consider any post-translational modifications" and so never builds the mature
chromophore (their SI, Note S4). We estimate it from the STEOM **transition
dipole**. Evaluated on the same structure, the two definitions differ:

## `venus_dimer.pdb`

crystal / vdW dimer, built from 1myw symmetry

| axis used | theta_ab (deg) | delta (deg) | offset vs dipole |
|---|---|---|---|
| **STEOM transition dipole** | — | **38.16** | — |
| OH->CB2 (structural proxy) | 125.19 | **27.40** | −10.75° |
| OH->CG2 (structural proxy) | 122.76 | **28.62** | −9.54° |
| CZ->CB2 (structural proxy) | 126.05 | **26.97** | −11.18° |
| *Cusick published* | — | *31* | — |

Separations on this same structure, by convention:
- `R_CB2_CB2_A` = 25.38 Å
- `R_CG2_CG2_A` = 22.65 Å
- `R_CZ_CZ_A` = 17.64 Å

## `tandem_dimer_production_cr2.pdb`

tandem dimer, MD starting structure

| axis used | theta_ab (deg) | delta (deg) | offset vs dipole |
|---|---|---|---|
| **STEOM transition dipole** | — | **39.61** | — |
| OH->CB2 (structural proxy) | 129.34 | **25.33** | −14.28° |
| OH->CG2 (structural proxy) | 127.72 | **26.14** | −13.47° |
| CZ->CB2 (structural proxy) | 129.06 | **25.47** | −14.14° |

Separations on this same structure, by convention:
- `R_CB2_CB2_A` = 25.53 Å
- `R_CG2_CG2_A` = 22.85 Å
- `R_CZ_CZ_A` = 17.73 Å

## Verdict

On the crystal dimer — the one structure where we and Cusick are
unambiguously looking at the same coordinates — switching from the
transition dipole to their structural proxy moves delta from
**38.16° to 27.40°**, i.e. by 10.8°, landing near their
published 31°. The separation agrees
to the digit: our CB2–CB2 = 25.38 Å against
their Table S2 value of 25.4 Å.

So most of the apparent geometric disagreement is a difference in what the
angle is *between*, not where the atoms are. Comparing our dipole-derived
delta directly against their proxy-derived delta is not like for like.

## Two further cautions about their AlphaFold numbers

- Table S2 lists five AlphaFold3 structures with delta = 9, 9, 8, 15, 9°. Four of the five are 8–9°; only
  structure #4 gives 15°. The main text quotes **15°** — the single value
  that overlaps their spectroscopic 14–20° range. The modal prediction, 9°,
  does not overlap it at all.
- Those same five structures give couplings of 35, 31, 43, 29, 30 cm⁻¹ — a
  spread of 1.48×
  from one prediction run. That is the signature of a poorly determined
  inter-domain orientation, which is expected for two barrels joined by a
  flexible 33-residue linker.

## What this does not resolve

- Even under their proxy our tandem sits at ~25°, still above their
  spectroscopic 14–20° and well above the modal AlphaFold 9°. The definition
  accounts for most of the gap, not all of it.
- The proxy offset is itself structure-dependent; it is measured here on two
  static structures, not averaged over the MD ensemble.
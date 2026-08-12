# Deep audit of the tandem-dimer geometry

The inter-chromophore transition-dipole angle computed here (~100 deg,
|cos alpha| = 0.19) disagrees with what Nguyen's limiting-anisotropy drop
requires (|cos alpha| = 0.660). Because the absorption first-moment shift is
exactly `J cos(alpha)`, and the exciton CD rotational strength carries
`R_AB . (mu_A x mu_B)`, that angle propagates into two of the three
observables the round-2 argument rests on. This audit works through every step
that sets it.

## Ruled out

**1. Rigid-density placement.** Three mutually independent routes agree:

| route | angle | |cos alpha| |
|---|---|---|
| rigid STEOM density, Kabsch on 19 CR2 heavy atoms (production, n=1000) | 100.78 +/- 2.72 | 0.187 |
| per-frame QM/MM TDDFT dipoles, computed fresh by TeraChem (n=88) | 96.84 +/- 3.16 | 0.120 |
| CR2 long axes alone, no dipoles at all (n=103) | 98.60 +/- 3.02 | 0.149 |

The angle is a property of the structure, not of how a density was placed.

**2. The wrong crystal dimer.** `venus_dimer.pdb` carries the header
"GENERATED CLOSEST DIMER WITH LATTICE SHIFT", which sounds like a
packing-contact selection. It is not: its chain B centroid is
`[17.82, 79.12, 1.53]`, identical to the image of 1MYW chain A under the
author-determined `REMARK 350 BIOMT2`. It is the deposited biological unit.
(The file does have a real defect: a literal `\n` glues atom 1 onto the HEADER
line, so strict parsers silently drop `MET0 N`.)

**3. The transition-dipole direction inside the chromophore.** In the shared
monomer frame the STEOM-CCSD and CAM-B3LYP transition dipoles agree to
**3.37 deg**, and each lies within 3.5 deg of the phenolate-O -> imidazolinone
long axis. A ~12 deg intramolecular tilt would have reconciled everything; it
is not there, at two very different levels of theory.

## The geometry, exactly

The 1MYW biological dimer is a clean **C2** dimer (BIOMT2 is a proper 180 deg
rotation, det +1, trace -1, axis `[0.866, -0.5, 0]`). For a C2 dimer

    cos(alpha) = 2 cos^2(theta) - 1

with theta the angle between a chromophore dipole and the C2 axis. This holds
to machine precision in the crystal: theta = 53.28 deg gives cos(alpha) =
-0.2850, matching the directly computed long-axis angle of 106.56 deg.

`alpha` is therefore controlled by a single interface parameter, and

| theta | cos(alpha) | alpha |
|---|---|---|
| 47.6 (v3 ff19SB/OPC) | -0.10 | 95.8 |
| 50.2 (production ff14SB) | -0.19 | 100.8 |
| **53.3 (crystal)** | **-0.285** | **106.6** |
| 65.6 (what the anisotropy needs) | -0.659 | 131.2 |

## Two findings

**(a) The MD moves the interface the wrong way, in both force fields.**
Crystal 106.6 -> 100.8 (ff14SB) -> 95.8 deg (ff19SB/OPC). Both trajectories are
flat over the 1 ns production window (net drift -0.4 and +5.1 deg), so the
offset is established during equilibration, not accumulated later. Whatever the
truth is, the simulation is not sitting on the crystal interface.

**(b) The crystal register puts the linker under tension.** The tandem's
33-residue linker must span C-term(barrel 1) -> N-term(barrel 2) =
**54.2 A** in the crystal register, against ~38 A RMS end-to-end for a
33-residue disordered chain. Across all 1000 v3 production frames the linker
sits at **52.3 +/- 0.7 A** — pinned, with essentially no slack, for a
nominally flexible tether. The tandem was *built* in the crystal register and
1 ns cannot re-dock two beta-barrels, so an alternative register has never
been sampled. This is the leading candidate for the discrepancy and the one
worth testing.

## What a 12 deg correction would and would not fix

Putting the dipole on a cone of half-angle theta about the C2 axis and scanning
the azimuth: the minimal change reaching cos(alpha) = -0.66 is a **12.3 deg**
tilt (theta 53.3 -> 65.6 at fixed azimuth). That would reconcile the anisotropy
and the red shift. It leaves the chirality triple product **negative**
(-5.05 -> -3.97 A, a 21% reduction). Reaching a positive triple product needs a
40-76 deg tilt, which is not a plausible structural error.

So the angle discrepancy and the CD-handedness discrepancy are **not the same
problem**, and should stop being treated as one.

## The handedness is probably a subtraction-order sign error

Nguyen define (main text, CD data analysis): "Difference CD spectra were
obtained by subtracting the TD from TDX", i.e.

    Delta[theta] = [theta]_TDX - [theta]_TD

The excitonic couplet exists only in TD, so the couplet seen in `Delta[theta]`
is **minus** the couplet of the dimer itself. In `Delta[theta]` for dVenus-TD
the bands are negative at 511.8 nm and positive at 521.7 nm (their Fig. 4d /
Table S3: A1 = +3.59 on the low-wavenumber component). Inverting, the TD
dimer's own CD is **positive to the blue, negative to the red** — a
**negative** exciton chirality, which is the sign we compute in every frame of
both ensembles.

Before attributing the handedness to the geometry, check whether the
comparison was made against `Delta[theta]` rather than against the TD couplet.
See [[venus-absolute-cd-handedness-disagreement]].

## Both open questions above were settled on 2026-08-12

**The linker-tension hypothesis is refuted.** A two-arm 36.6 ns experiment
(`overnight_linker_release/`) steered the linker end-to-end from 52.7 to 41 A,
held, then released. Alpha moved the *wrong* way (down to ~86 deg, away from
the experimental requirement); after release the linker kept ~4 A of slack, so
the stored tension is real, but alpha only recovered toward the control. The
unbiased control meanwhile drifted *up*, 99.3 -> 104.8 deg, toward the crystal.

**And finding (a) above inverts: the crystal is right and the MD is wrong.**
The superradiant lifetime drop (57 +/- 4 ps) is an independent handle on
alpha that neither paper inverts. Thermally corrected,

    Dtau/tau = Phi * tanh(Omega/2kT) * (2J/Omega) * (-cos alpha)

A faster lifetime demands **cos alpha < 0**, which resolves the anisotropy's
cos^2 branch degeneracy in favour of **131.3 deg, not the 48.7 deg Nguyen
quote** — so both experiments agree on an obtuse angle, the same side as this
audit's 96-110 deg. It also imposes a floor |cos alpha| >= 0.212 that no
detuning can rescue. The 1MYW biological dimer (alpha 110.44, cos -0.349,
J 30.82) then reproduces the measured 57 ps *exactly* at |Delta| = 570 cm^-1,
against a computed detuning of 549-581 cm^-1, with no free parameters, while
every MD ensemble under-predicts in proportion to how far it has relaxed away
from the crystal interface. See `exciton_observables/README.md`.

The 12.3 deg tilt discussed above is therefore not needed: it was an artifact
of comparing against the wrong branch (48.7 deg) of the anisotropy fit.

## Reproduce

```bash
python3 notes/interdipole_angle_probe.py tandem_dimer.pdb venus_dimer.pdb
```

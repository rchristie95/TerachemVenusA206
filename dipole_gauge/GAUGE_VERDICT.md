# The STEOM transition moment fails its own consistency check

2026-08-27. Same geometry (`geom_cthrp.xyz`), same point-charge embedding
(`field.pc`), ORCA 6.1.1. Length and velocity gauges must agree for an exact
wavefunction in a complete basis, so their ratio is an internal check that needs
no experimental input.

| method | basis | E (eV) | f (length) | f (velocity) | vel/len | \|mu\| length |
|---|---|---|---|---|---|---|
| STEOM-CCSD | def2-SVP | 2.4724 | 0.9163 | 1.2273 | 1.34 | 9.886 D |
| STEOM-CCSD | def2-SVPD | 2.3272 | 0.8488 | 1.2833 | **1.51** | 9.807 D |
| TDDFT CAM-B3LYP | def2-SVPD | 3.0681 | 0.8795 | 0.8568 | 0.97 | 8.694 D |
| TDDFT CAM-B3LYP | def2-TZVPD | 3.0574 | 0.8729 | 0.8690 | **1.00** | 8.677 D |

**The basis is exonerated.** TDDFT reaches gauge agreement (0.97 -> 1.00,
improving toward exact) in the same bases where STEOM disagrees by 1.34 -> 1.51
and *worsens*. The earlier reading -- that the STEOM dipole was
"basis-converged to -0.54%" -- was measuring the stability of the length-gauge
number while the underlying transition moment was not converged at all.

So the defect is in the STEOM transition moment, not the basis and not (as
tested earlier) the left/right normalisation, which accounts for only 1.3%.

## What |mu| to use

TDDFT gives **8.68 D**, gauge-consistent and basis-converged (8.694 -> 8.677,
-0.2%). Propagating through `J ~ |mu|^2`:

    J = 32.82 x (8.677/9.807)^2 = 25.7 cm^-1

against 32.82 from the STEOM dipole. Superradiance at the corrected static
detuning wanted J = 21-26, so 25.7 lands at the top of that range -- a large
improvement on 32.82, which was excluded at 5.6-8.2 sigma.

## The honest limit

Gauge agreement is necessary, not sufficient. TDDFT places the excitation at
3.06 eV (405 nm) against STEOM's 2.33 eV (532 nm) and experiment's ~515 nm: it
gets the DIPOLE self-consistently and the ENERGY badly wrong, while STEOM does
the reverse. A method that is wrong by 0.7 eV on the energy can be consistently
wrong on the moment too.

**Neither |mu| is trustworthy on its own.** What the test does establish:

1. the STEOM value fails an internal check that worsens with basis, so it should
   not be quoted as the converged number it was taken for;
2. the true value is very unlikely to be as high as 9.81 D -- three independent
   experimental lines (7.5-7.9 D) and the one gauge-consistent calculation
   (8.68 D) all sit below it;
3. J should be quoted as a range, **21-26 cm^-1**, not as 32.82.

Resolving it properly needs a method that gets both the energy and the moment
right -- EOM-CCSD with explicit left/right moments, or CC2 -- not a bigger basis.

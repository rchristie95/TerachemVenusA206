# ORCA Tyr-link-cap audit

This control compares the manuscript's retained ORCA monomer TDDFT calculation
against an otherwise identical calculation in which only the Tyr203 CG link
hydrogen is reflected through its parent CG atom.

- `orca_correct_cap_tddft`: byte-identical input geometry and MM field from
  `neo_model/orca_dft`.
- `orca_reversed_tyr_cap_tddft`: identical except atom 44 is moved from
  `(41.89052226, 40.94235807, 40.60983043)` Angstrom to
  `(41.79147774, 42.34764193, 38.94616957)` Angstrom.
- Tyr203 CG is atom 30 at `(41.841, 41.645, 39.778)` Angstrom. Both cap
  positions are 1.09 Angstrom from CG, on opposite sides.
- Method: ORCA 6.1.1, TDA-wB97X-D3/6-311G**, RIJCOSX, TightSCF, five roots,
  charge -1, singlet, and the unchanged 2,350-charge embedding field.

The retained STEOM and ORCA TDDFT geometry has the correct cap: the unit vector
from CG to the link H has dot product 1.000 with the unit vector from CG to the
omitted CB atom. The reversed-cap control represents the bug found only in the
new trajectory-frame extractor.

## Results

Both ORCA 6.1.1 calculations used eight MPI processes and terminated normally.
The correct-cap run took 3 minutes 6 seconds.

| Quantity | Correct retained cap | Deliberately reversed cap | Reversed - correct |
|---|---:|---:|---:|
| Brightest root | 2 | 2 | |
| Brightest-root energy (eV) | 3.435623 | 3.419290 | -0.016333 |
| Brightest-root wavelength (nm) | 360.878 | 362.602 | +1.724 |
| Brightest-root oscillator strength | 1.032983 | 0.728895 | -0.304088 |
| Root-3 wavelength (nm) | 358.602 | 360.648 | +2.046 |
| Root-3 oscillator strength | 0.253484 | 0.572858 | +0.319374 |
| Roots 2+3 strength-weighted energy (eV) | 3.439920 | 3.427444 | -0.012476 |
| Roots 2+3 strength-weighted wavelength (nm) | 360.428 | 361.740 | +1.312 |
| Roots 2+3 summed oscillator strength | 1.286467 | 1.301753 | +0.015286 |

The manuscript reports 360.42 nm and oscillator strength 1.03. The retained
correct-cap ORCA rerun gives 360.878 nm for the single brightest root
(`f=1.032983`), while the strength-weighted centroid of the two mixed bright
roots is 360.428 nm. The latter agrees with 360.42 nm to 0.008 nm.

The reversed cap does not simply translate one isolated state. It mixes roots
2 and 3, reducing root 2 from `f=1.0330` to `f=0.7289` and increasing root 3
from `f=0.2535` to `f=0.5729`. The two-root bright-band centroid redshifts by
1.312 nm (0.012476 eV or 100.63 cm-1).

Machine-readable results are in `summary.json`; complete ORCA outputs and exit
codes are retained in each control directory.

# v1 ensemble — trajectory with known MD setup defects

**Status: superseded. Retained for comparison only. Do not use v1 numbers for
publication without the caveats below.**

Frozen 2026-08-10 at **87 usable frames** of a planned 200-frame 5 ps grid. The
QM campaign was stopped deliberately once the underlying trajectory was found
to have setup defects; the remaining ~113 frames were not run because they
would have measured the same flawed system.

## What v1 is

| item | value |
|---|---|
| trajectory | `tc_tandem_nvt_rerun_20260722_retry2/tandem_nvt_1000.dcd` |
| topology | `tc_tandem_nvt_rerun_20260722_retry2/solvated_protonated.pdb` |
| system | 69,456 atoms, 491-residue VenusA206 tandem dimer, 20,524 waters, 67 Na + 56 Cl |
| QM | CAM-B3LYP/6-311+G\*\*, TDA, 7 roots, charge −1, TeraChem 1.97B |
| embedding | full-system point charges, link-only exclusion (both introduced during this work) |
| QM directories | `results/hi_camb3lyp_frame_*` (5 ps grid), `results/ct_camb3lyp_frame_*` (5 ps correlation test) |

## Defects found in the trajectory (audit, 2026-08-10)

**1. System is not electrically neutral: net −2.000000 e.**
`run_nvt.py` calls `addSolvent(neutralize=True, ionicStrength=0.15 M)` against a
force field in which the generic fallback template assigns **charge 0.0 to every
CR2 atom** (`run_nvt.py:453`). The real −1 e per chromophore is transplanted
afterwards directly onto the `NonbondedForce` (`run_nvt.py:1858`), with no
re-neutralisation. The 123 ions present were placed only to reach 0.15 M ionic
strength, never to neutralise the chromophores.

Consequence: **the two anionic chromophores have no counter-ion atmosphere**, and
PME compensates the net charge with a uniform background rather than local
screening — at exactly the two sites whose energies are being measured.

**2. Box is 5.2% under-dense; no NPT stage exists.**
Density 0.9530 g/cm³ against 1.00–1.02 for a solvated protein; water number
density 27.46 nm⁻³ against a bulk 33.4, i.e. **water is ~18% under-dense**.
There is no `MonteCarloBarostat` anywhere in `run_nvt.py`; the protocol is
minimise → NVT → minimise with the box frozen at the `addSolvent` padding.
Under-dense water weakens dielectric screening, which is a plausible mechanism
for inflated site-energy fluctuations.

**3. No heating ramp and no velocity seeding.**
No `setVelocitiesToTemperature` call exists; on a fresh run production begins
from zero velocities and the thermostat heats the system inside the analysed
trajectory. At 1 ps⁻¹ friction this affects only the first few ps of 1 ns, so it
is minor here, but nothing trims it.

**4. Linker cysteine not mutated.**
Nguyen et al. mutated the linker Cys→Gly to prevent spontaneous disulfide
formation. Our construct retains **Cys246**. It sits ~60 residues from either
chromophore so it should not shift site energies directly, but our linker is not
their linker.

**5. Minor.** `addSolvent(model="tip3p")` places water geometry while the system
is parameterised with `amber14/tip3pfb.xml`. Switching function, Ewald error
tolerance and dispersion correction are never set explicitly. Histidine
tautomers are left to PDBFixer defaults, uninspected and unlogged — relevant
because His147/His409 each contribute ~+670 cm⁻¹ to the site energy.

Not a defect, but noted: the trajectory is not fully relaxed — Cα RMSD still
rising at 0.75 Å/ns, chromophore separation drifting −0.83 Å/ns. Tested and found
NOT to propagate into the site energies: the detuning spread is statistically
indistinguishable between the first and second halves of the trajectory
(sd 659.7 vs 734.1 cm⁻¹, overlapping bootstrap intervals).

## v1 results, and how much to trust them

| quantity | v1 value |
|---|---|
| mean \|Δ\| | 549.3 cm⁻¹ |
| sd Δ | 691.1 cm⁻¹ |
| 2\|J\| | 65.63 cm⁻¹ |
| ⟨mixing⟩ | 0.18 |
| chromophore separation | 25.31 ± 0.43 Å (crystal: 25.39) |
| \|T\| | 76.82 ± 15.55 Å·au² (crystal: 94.27) |
| detuning decorrelation | ≤ 5 ps |

**Robust to the defects:** |Δ| ≫ 2|J|. The ratio is 8.4×, and none of the defects
above plausibly manufactures an order of magnitude. The interface geometry is
also validated independently against the crystal dimer.

**Not robust:** the *magnitude* of the disorder, and anything derived from it —
the apparent-coupling value, the CD couplet width, and the comparison of
absolute CD intensity to experiment. Defects 1 and 2 both act on the local
electrostatics at the chromophores.

## What v2 changes

Rebuild with: neutralisation performed *after* the CR2 charges are correct; an
NPT stage to relax the density; explicit velocity seeding; and the linker
Cys246→Gly mutation. v2 directories will be tagged `v2_` and the two ensembles
compared frame-for-frame on the same 5 ps grid.

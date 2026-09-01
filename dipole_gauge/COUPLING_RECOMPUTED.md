# The coupling, recomputed rather than rescaled (item 2)

2026-08-27. `build_tzvpd_coupling.py` (new), `build_tddft_coupling.py` (existing).
All three densities are placed on the **identical** dimer geometry
(`venus_dimer.pdb`, OLD_MONOMER frame) by the same pipeline, and all three QM
runs use the byte-identical 44-atom geometry md5 `6d55e271...`.

| density | \|mu\| D | gauge v/l | J_TDC | J_PDA | TDC/PDA |
|---|---:|---:|---:|---:|---:|
| STEOM-DLPNO-CCSD / def2-SVPD | 9.81 | 1.512 | 30.45 | 24.81 | 1.227 |
| wB97X-D3 / 6-311G** (TDA) | 8.90 | -- | 22.06 | 19.53 | 1.130 |
| **CAM-B3LYP / def2-TZVPD (TDA false)** | **8.68** | **0.995** | **21.50** | 18.35 | **1.172** |

> **2026-08-30, J_PDA corrected.** The two TDDFT rows previously read 18.12 /
> 1.217 and 16.48 / 1.305. `build_tddft_coupling.py` and
> `build_tzvpd_coupling.py` placed the point-dipole site centroids at the
> *unweighted* mean of the threshold-retained voxels, while the STEOM row (and
> `compute_steom_coupling.py`) used the |q|-weighted centroid for the same
> quantity. The unweighted mean moves 0.13 A over a 100x threshold change; the
> weighted centroid moves <0.001 A. All three rows are now |q|-weighted and
> directly comparable. **J_TDC is unchanged in every row** -- the Coulomb
> double sum never used a centroid -- so the headline coupling is unaffected.

The third row is the one to quote: it is computed from a density whose
length/velocity gauges agree to 0.5%, with **no |mu|^2 rescaling anywhere**.

    single-geometry 21.50  x  (32.82/30.454 ensemble factor)  =  23.2 cm^-1

so the ensemble coupling from a gauge-verified density is **~23 cm^-1**, sitting
inside the quoted 21-26 range and near its lower half.

## Rescaling is NOT equivalent to recomputing

STEOM couples 1.383x more strongly than the 6-311G** TDDFT density while the
dipole ratio squared is only 1.214 -- a residual **1.139x that is density
SHAPE**, not magnitude. So the two routes do not coincide:

| route | renormalised to 7.84-8.68 D |
|---|---|
| recomputed from the TDDFT density | 17.1 - 21.0 cm^-1 |
| rescaled from the STEOM ensemble | 21.0 - 25.7 cm^-1 |

They meet only at ~21. The directly computed gauge-verified value (23.2
ensemble-equivalent) lands between them, which is the honest place to sit.

## The reviewer's point is untouched by any of this

| source | TDC/PDA |
|---|---|
| Reviewer 2, Q-Chem gas phase | 1.127 |
| STEOM, this geometry | 1.227 |
| wB97X-D3/6-311G** | 1.217 |
| CAM-B3LYP/def2-TZVPD | 1.305 |
| production ensemble (n = 1000) | 1.187 +/- 0.011 |

Three different densities, two different functionals, one coupled-cluster
method, and an independent code all put the near-field correction at 1.13-1.31.
Nothing here is near 4.

## A trap worth recording

`orca_plot` state numbering depends on TDA. For a `TDA false` run the states must
be entered as **2*state - 1** (the X+Y vectors); `2*state` gives X-Y. Entering
the plain root number 3 for the def2-TZVPD run returned the X+Y vector of state
**2** -- a dark state -- which showed up as a cube dipole of 0.302 au against the
3.414 target and **cos(theta) = -0.598** against the ORCA transition dipole. The
resulting J of 1.85 cm^-1 with TDC/PDA = 0.41 was pure artefact. The two
diagnostics printed by the build script (raw-vs-target |mu|, and cos(theta))
catch this immediately -- always check cos(theta) ~ +1 before believing a
coupling. The earlier wB97X-D3 run was TDA, where the numbering is direct, so
that script was correct as written.

# Harvested site-energy ensembles

Compact extracts of the per-frame QM/MM TDDFT trees in `../`. These arrays are
what every round-2 number in the manuscript is computed from — the per-frame
directories are ~10 GB of TeraChem scratch, these are 15 KB.

**Rescued into the repository 2026-08-11.** They had been written to a session
scratch directory under `/tmp` and existed nowhere else. A `/tmp` sweep would
have destroyed the round-2 result and forced a full recomputation of the QM
ensemble.

Regenerate with `../../scripts/parse_frame_ensemble.py --output <file>.npz`.

## Contents

Every `ens_*` / `ensemble_*` file carries the same ten keys, one row per frame:

| key | meaning |
|---|---|
| `frame` | frame index in the source DCD |
| `e_a_cm`, `e_b_cm` | bright-state site energies of chromophores A and B, cm^-1 |
| `mu_a_au`, `mu_b_au` | transition dipole vectors, atomic units, phase-aligned along the CR2 long axis |
| `r_a_ang`, `r_b_ang` | chromophore centroids, Angstrom |
| `triple_product` | R_AB . (mu_A x mu_B), the chirality factor T |
| `separation_ang` | centroid separation |
| `margin` | bright-state assignment margin (state-tracking confidence) |

| file | n | source trajectory | method |
|---|---|---|---|
| `ens_v2.npz` | 88 | v2 (`tc_tandem_nvt_v2`, corrected MD) | CAM-B3LYP/6-311+G**, TDA |
| `ens_final.npz` | 87 | v1 (`tc_tandem_nvt_rerun_20260722_retry2`, defective MD) | CAM-B3LYP/6-311+G**, TDA |
| `ensemble_camb3lyp.npz` | 23 | v1, early partial | CAM-B3LYP/6-311+G**, TDA |
| `ct.npz` | — | v1, 5 ps correlation-time test | CAM-B3LYP |
| `pol_v2.npz` | — | v2, perturbative QM/MMpol polarization pass | see `../../polarizable_embedding.py` |
| `surrogate_fit.npz` | — | classical detuning-surrogate fit (rejected, R^2 = 0.063) | n/a |

`ens_v2` vs `ens_final` is the v2-vs-v1 comparison that showed the MD rebuild
left the disorder statistically unchanged (Mann-Whitney p = 0.45). Note that
comparison is underpowered — see `../../PROVENANCE_v1_DEFECTIVE_MD.md`.

`ens_live.npz` / `ens_now.npz` were intermediate snapshots of `ens_v2` taken
while the campaign was still running and were deliberately not rescued.

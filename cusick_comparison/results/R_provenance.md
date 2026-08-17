# Provenance of the 27.6 Å centroid separation

Read-only investigation. **No manuscript was edited.** This records what could and
could not be established about a number in a live submission.

## The claim under investigation

`manuscript/Submit-JPCL-21April2026.tex` states a chromophore centroid separation
of 27.6 Å in three places — the abstract (:59), the methods (:192) and the
discussion (:217) — always alongside the `J = 74.38` / `J_PDA = 13.31` / "5.6×"
near-field claim. The methods line also reports an inter-dipole angle of 92.85°.

## What the tree actually contains

| Source | R (Å) | Inter-dipole angle | Provenance |
|---|---|---|---|
| Production tandem ensemble, n=1000 | **24.687 ± 0.320** | 100.778 ± 2.721° | `coupling_geometry.npz`, verified by closure test |
| Static crystal / vdW dimer | **25.209** | 103.686° | `coupling_paper_steom_thermal/dipole_geometry.json` |
| `venus_dimer.pdb` CR2 heavy-atom centroids | 25.388 | — | measured directly |
| **Submitted JPCL manuscript** | **27.6** | **92.85°** | **matches nothing in the tree** |

No result file, JSON, CSV or npz anywhere in the repository contains 27.6 as a
separation. Neither does any file contain 92.85 as an inter-dipole angle — every
grep hit for "92.85" is a coincidental substring of an unrelated float.

## Hypothesis tested and REJECTED: different density → different centroid

R is defined as the separation of the two density point-cloud centroids
(`coupling_dcd_steom.py:182`, `density_origin = points_full.mean(axis=0)`), so a
different transition density genuinely gives a different R. The legacy JPCL
numbers used a TDDFT S2 cube; the production numbers use a cap-masked STEOM
density. That could in principle explain the gap.

It does not. Measured directly, in the shared "oldframe":

```
tddft_transdens_specnorm_oldframe.npz   n=285975  centroid=[23.951, 89.275, -23.315]
steom_transdens_specnorm_oldframe.npz   n=157474  centroid=[24.267, 88.921, -23.160]
|centroid difference| = 0.499 Å
```

A 0.50 Å shift in the monomer-frame centroid can move the dimer separation by at
most ~1.0 Å, and less in practice because the same shift applies to both
monomers under the dimer's near-C2 symmetry. Closing the 2.4 Å gap between 25.21
and 27.6 is out of reach. **The density choice does not explain 27.6 Å.**

## What 27.6 *is* elsewhere in the repo

The numeral 27.6 does occur, twice, as a **coupling in cm⁻¹** — never as a length:

- `J_pda_cm` mean over the production ensemble = **27.642 cm⁻¹**
- `manuscript/JPCB_tandem.tex` quotes "PDA 27.6" (recorded at `CLEANUP_MANIFEST.md:22`)

So the same numeral appears as a cm⁻¹ coupling in one manuscript and an Ångström
separation in another. That is suggestive of a units/quantity transcription, but
it is **not proof**, and the chronology does not obviously support a simple
copy: the JPCL submission predates the JPCB draft.

## Verdict

Both JPCL geometry numbers — R = 27.6 Å and the 92.85° inter-dipole angle — are
**not reproducible from anything currently in the tree**, and both travel with
the superseded 74.38 / 13.31 pair that `reproduce_paper.py:115` already records
as corrected (74.38 → 20.83). The most defensible reading is that the entire v1
geometry block belongs to the same lost or defective calculation as the coupling
numbers it accompanies, consistent with the separate v1 units bug (v1's 74.38
came from applying no conversion at all, ×1.8897, a *different* error from the
×3.5711 that affected production).

Whether it is a transcription slip or a geometry that no longer exists, the
conclusion for the arXiv v2 correction is the same:

> **Do not carry 27.6 Å or 92.85° into v2.** Quote the reproducible values:
> R = 24.69 ± 0.32 Å and inter-dipole angle 100.78 ± 2.72° for the production
> tandem ensemble, or R = 25.21 Å / 103.69° for the static crystal dimer —
> stating which construct each belongs to.

## What would settle it

Only an artefact outside this repository: the original working notebook, scratch
directory, or PyMOL session from April 2026. If none survives, the number is
simply unverifiable and should be replaced rather than defended.

# arXiv 2605.00027 v1 → v2: what changes and why

v1 is live and citable. Its headline — a coupling of 74.38 cm⁻¹, "5.6 times
stronger" than a point-dipole estimate, attributed to near-field multipolar
effects — is an artefact of a unit-conversion error. This document catalogues
every claim that changes, the corrected value, and the artefact in the repository
that establishes it.

> **v2 source: `Submit-arXiv-v2.tex` is `manuscript/JPCB_tandem_round_2.tex`
> plus a correction note in the abstract.** An earlier attempt patched the v1
> JPCL letter instead; that was abandoned because the round-2 revision already
> contained the corrected analysis in fuller form — 18 of the 18 claims below
> against 3, at 24 pages against 10 — and because the v1 letter is organised
> around a headline that no longer stands. Same title, same author list, so the
> replacement is the same paper revised rather than a different one. Keeping a
> separately patched arXiv text would have recreated exactly the two-divergent-
> drafts problem that produced the contradiction this correction resolves.

**Every number below is reproducible from committed scripts.** The production
ensemble is `coupling_nvt_production_cr2_1000_20260721/` (1000 frames); the
derived-parameter comparison is `cusick_comparison/` (six scripts).

## The error

The Coulomb double sum is evaluated with coordinates in Ångström, so the
conversion to atomic units carries `BOHR_TO_ANGSTROM`, not its reciprocal.
v1 applied *no* conversion at all, inflating J by 1.8897. (Production code
separately carried a different form of the same mistake, ×3.5711 — the two are
not the same factor, so **v1 numbers must not be rescaled by 3.5711.**)
Corrected on v1's own static geometry, its calculation gives 39.36 cm⁻¹ with a
TDC/PDA ratio of 2.96 — still not 74.38, and still not 5.6.

But v1's geometry is also not the right comparison: it paired a TDDFT-cube
transition density on a single static crystal frame against a point-dipole value
from a different structure. The correct like-for-like pair, over the same 1000
frames of the same trajectory with the same density, is below.

## Claim-by-claim

| v1 claim | location | corrected | source |
|---|---|---|---|
| J = 74.38 cm⁻¹ | abstract, §Results | **32.82 ± 1.55 cm⁻¹** (screened, ε = 1.77) | `coupling_distribution.json` |
| PDA = 13.31 cm⁻¹ | abstract | **27.64 ± 1.28 cm⁻¹** | `coupling_samples.csv`, column `J_pda_cm` |
| "5.6 times stronger" | abstract, §Results | **1.187 ± 0.011** | `cusick_comparison/angle_inference_bias.py` |
| "near-field multipolar effects … dominant" | abstract, intro, §Results | **PDA is accurate to ~19% at this separation** | as above |
| centroid separation 27.6 Å | abstract, §Results | **24.69 ± 0.32 Å** | `coupling_samples.csv`, `separation_A` |
| Davydov splitting 2\|J\| ≈ 149 cm⁻¹ | §Results | **2\|J\| = 65.6 cm⁻¹** | derived from J above |
| "better agreement with Nguyen 131–186 cm⁻¹" | §Results | that range is itself an artefact — see below | `cusick_comparison/cd_splitting_resolution.py` |
| "delocalised exciton superposition" | abstract | **detuning-dominated**: ⟨\|Δ\|⟩ = 576 cm⁻¹ ≫ 2\|J\| = 65.6, eigenstates localised, mean mixing 0.18 | `ens_v2_all.npz` (n = 95) |

### The 27.6 Å is probably a transcription of a coupling

No result file contains 27.6 Å. The production separation is 24.69 ± 0.32 Å.
Meanwhile `27.64` does exist in the data — as `J_pda_cm` in cm⁻¹. It also
appears in Cusick's Table S2 as a 2D coupling in cm⁻¹. A cm⁻¹ value appears to
have been transcribed into an Ångström slot.

### The CD comparison does not support the coupling either way

v1 argued that 2|J| ≈ 149 cm⁻¹ agreed better with an apparent 131–186 cm⁻¹ read
from CD. That reasoning does not survive:

- Reading the couplet **separation** as 2J requires |μ| ≈ 29 D or R ≈ 10 Å.
- Reading it as the exciton gap Ω = √(Δ² + 4J²) fails too: feeding Ω = 553.9
  through the actual lineshape produces extrema separated by **652 cm⁻¹**, not
  the 548 measured.
- The separation is a **non-monotonic, two-valued, lineshape-dominated** function
  of the gap. It falls from the zero-gap derivative limit, passes a minimum near
  ~190 cm⁻¹, then rises — so a measured separation maps to two candidate gaps
  and is not a usable estimator in either direction.

Cusick et al. hit the identical artefact in their two-photon polarisation
spectrum — features separated by 750 cm⁻¹ against an H–J splitting of 60–90
cm⁻¹ — identified it as the first derivative of a Gaussian with extrema at
FWHM/√(2 ln 2), and fitted **amplitudes** instead. Only the couplet amplitude,
suppressed by 2|J|/Ω, carries J.

## What replaces the headline

Three results that are stronger than the claim they displace:

1. **The point-dipole approximation is adequate for FP dimers at biological
   separations** — accurate to ~19% at 24.7 Å. Unglamorous, but it is a usable
   methods result, and it is the opposite of what v1 reported.

2. **The dimer is detuning-dominated, not delocalised.** ⟨|Δ|⟩ = 576 cm⁻¹
   against 2|J| = 65.6 cm⁻¹. The eigenstates are localised (mean mixing 0.18),
   which explains the spectroscopy without invoking a protected delocalised
   exciton.

3. **The dielectric screening convention is now the dominant systematic** —
   larger than the near-field correction v1 was reporting. This is a real
   methodological result and it is defensible; see the next section.

## Screening: the comparison that must be stated in both conventions

Our J is screened by ε = 1.77. Cusick state explicitly, after their eq 13:
*"We disregard any dielectric attenuation by the protein matrix because of the
nanometer scale of the system."*

| convention | this work | Cusick 2026 |
|---|---|---|
| screened (ε = 1.77) | 32.82 ± 1.55 | 18.1–22.6 |
| unscreened | 58.09 | 32–40 |

**The factor cancels from the ratio: we are 1.61× high in either convention.**
Any claim of agreement that compares a screened number to an unscreened one is
not a real agreement, and must not be made.

## The transition dipole is the dominant remaining discrepancy

Inverting Cusick's eq 13 for |μ| from their own (J, R, δ) gives **7.6–7.9 D** —
the experimental extinction/Strickler-Berg value. Our spec-normalised STEOM
density gives **9.81 D**. Since J ∝ |μ|², that single input accounts for most of
the residual 1.61× gap.

Independently, a basis-convergence sweep gives TDDFT |μ| converging to **7.84 D**
at def2-TZVPD, matching experiment. A STEOM/def2-TZVPD calculation testing
whether 9.81 D is basis-limited is **in progress**; v2 should not be posted with
a strong claim about |μ| until it lands.

## Recommended framing for v2

Post a short **"Corrected"** note in the abstract rather than waiting for the
journal revision to conclude. The coupling of the two is what keeps a wrong
number in circulation, and anyone who has already built on v1 needs the
correction to be findable.

The honest and more interesting statement is: *at this geometry the three
methods agree to within the screening convention, and the choice of ε is now the
dominant systematic.* Not: *we reproduce the measured coupling.*

## Open items to resolve or flag before posting

- STEOM/def2-TZVPD |μ| (running) — decides whether 9.81 D is a basis artefact.
- γ₀ = 38.7–44.4° against Cusick's directly measured 22°. Unexplained. Our value
  rests on *unrelaxed* state dipoles; no relaxed block exists in the tree.
- Residual δ gap: ~25° (ours, under their structural axis convention) against
  their 14–20° spectroscopic and 9° modal AlphaFold3. ~11–14° of the original
  apparent gap was purely definitional; the remainder is real and points at the
  tandem register.

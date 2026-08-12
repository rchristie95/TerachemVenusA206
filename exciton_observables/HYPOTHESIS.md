# Working hypothesis: J ≈ 30 cm⁻¹, and how every Nguyen observable follows

2026-08-12. Standing assumption (fixed by decision, and well supported): the
tandem coupling is **J ≈ 30 cm⁻¹**. The TDC production ensemble gives
32.82 ± 1.55, the candidate-2 ensemble 32.25 ± 1.80, and the two registers are
indistinguishable in J — so the coupling is settled and cannot be the thing
that reconciles anything. What remains free is the geometry (α) and the
site-energy detuning (Δ).

## The hypothesis

**H: the tandem adopts a candidate-2-like obtuse register whose true
inter-dipole angle is α ≈ 131°, and the site-energy detuning in that register
is large and broad, ⟨|Δ|⟩ ≈ 1000–1300 cm⁻¹ with comparable spread and ps-scale
fluctuations.**

Two subsidiary claims carry the load:

1. **The ~10° force-field relaxation is an artifact, not information.** ff19SB
   shaves ~10° off both interfaces it has been handed (crystal 110.4 → 102.3;
   candidate-2 131.5 → 121.1, each a stable basin). The real register holds
   the as-built ~131°. This is the hypothesis's weakest link and is listed as
   test D.
2. **The detuning is register-specific.** Δ is driven by the differing
   electrostatic environments the two barrels present to each other's
   chromophore. Every computed detuning to date (576 ± 71 cm⁻¹) came from
   crystal-register MD. If the structure is candidate-2-like, that number is
   the right answer to the wrong question. The hypothesis predicts the
   candidate-2 ensemble gives ⟨|Δ|⟩ well above 576, toward ~1000–1300.

## What it reproduces, observable by observable (numbers verified at J = 30)

| # | Nguyen observable (dVenus) | mechanism under H | check |
|---|---|---|---|
| 1 | Absorption red shift, −19.7 ± 4 cm⁻¹ after electrostatic correction | J·cos α = −19.8 cm⁻¹; exactly Δ-independent (first-moment identity) | ✓ |
| 2 | Limiting anisotropy R0 0.52 → 0.30 | transfer is ultrafast (Marcus 1–4 ps ≪ 26.5 ps IRF even at Δ = 1310); R0 = 0.30 requires ⟨cos²α⟩ = 0.436, i.e. α = 131.3° exactly; held-131.5° basin with ±0.09 fluctuations gives R0 = 0.303 | ✓ |
| 3 | Superradiance, Δτ = 57 ± 4 ps on 3.026 ns | Φ·tanh(Ω/2kT)·(2J/Ω)·\|cos α\| with \|cos α\| = 0.66 reproduces 57 ps at Ω ≈ 1190 cm⁻¹ — inside the predicted detuning range | ✓ |
| 4 | CD couplet, apparent splitting 261.58 cm⁻¹ (constrained fit) | at typical Δ ≈ 1200 the couplet amplitude is suppressed 20× (2J/Ω = 0.05); the visible couplet is dominated by the near-degenerate tail of the broad Δ distribution, and Nguyen's constrained 3-Lorentzian fit of that unresolved, tail-weighted feature returns an apparent splitting that is not the exciton gap. Their own hedges ("apparent", "merger of multiple unresolved split components") anticipate this | forward-model (test B) |
| 5 | Third CD band at 481 nm | vibronic, S ≈ 0.35–0.5 at ω ≈ 1450 cm⁻¹; unchanged from the existing treatment | ✓ (existing) |
| 6 | TDX controls: τ = 3.026 ns, R0 = 0.52, no couplet | TDX is a single chromophore; 0.52 fixes the 14° intramolecular abs/em angle, which enters observable 2 through the r_TDX baseline | ✓ trivially |

The 5× detuning "trilemma" (253 / 576 / 1310 cm⁻¹) dissolves rather than being
won by one side: **1310 becomes the prediction** for the true register, **576
is structure-conditional** (right calculation, wrong structure), and **253 is a
fit-protocol artifact** to be demonstrated, not assumed.

## What would kill it, and the tests in order of leverage

**A. Recompute the QM/MM site-energy detuning on the candidate-2 ensemble.**
The 20.5 ns trajectory exists (`tc_candidate2_ff19sb_opc/candidate2.dcd`), the
per-frame machinery exists (`terachem_site_energy_cd/`). Prediction:
⟨|Δ|⟩ ≫ 576, toward ~1000. **If it comes out ≈576 again, the hypothesis is
dead** — there would then be no consistent (α, Δ) pair and the failure moves to
the two-state framework itself. This is the decisive, GPU-worthy computation.

**B. Forward-model the CD end-to-end.** Generate ensemble spectra from
per-frame (Δ, J, geometry), apply Nguyen's *exact* constrained fit protocol
(couplet centres pinned symmetrically about the TDX absorption peak, δ free),
and check that the apparent δ ≈ 130 cm⁻¹ and a sensible amplitude scale emerge
from a broad-Δ ensemble whose true 2J is 60. This converts claim 4 from
plausible to demonstrated (or refutes it).

**C. Compute the two-photon absorption bias q.** Excitation at 20,000 cm⁻¹ is
blue of both site origins, so the higher-energy site absorbs preferentially and
q (the transferred fraction) exceeds ½ — which pushes the anisotropy-required
angle slightly *more* obtuse. Refinement of observable 2, not a rescue.

**D. Attack the ~10° force-field systematic.** Options: a third independent
starting register, interface restraints calibrated on the crystal dimer's
B-factors, or a polarizable force field on a short window. Until then the
131° ← 121° attribution rests on two observations of the same drift.

## Status

- [x] arithmetic consistency at J = 30 (this file)
- [ ] test A — detuning on candidate-2 ensemble
- [ ] test B — CD forward model through Nguyen's fit protocol
- [ ] test C — 2P absorption bias
- [ ] test D — force-field systematic

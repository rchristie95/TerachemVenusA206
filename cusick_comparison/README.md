# Cusick et al. 2026 comparison — derived parameters, no 2PA required

We never ran two-photon absorption. Omega is Cusick's instrument, not their
result: their analysis outputs ordinary molecular parameters, and this pipeline
produces all of them natively. The comparison therefore happens at the level of
derived parameters. Derivations and conventions are in
[`notes/cusick_derived_comparison.md`](../notes/cusick_derived_comparison.md).

Everything here reads existing artefacts. No GPU, no new QM, no new MD.

## Results

All Cusick values below are verified against the published paper
(JPC A [10.1021/acs.jpca.6c02663](https://doi.org/10.1021/acs.jpca.6c02663)) and
its SI, in `notes/Cuisick2026/`.

| Quantity | This work | Cusick | Verdict |
|---|---|---|---|
| delta, tandem | **39.61 ± 1.36°** | 14–20° (spectroscopic), 15° (AlphaFold3) | **structural disagreement** |
| delta, crystal/vdW | **38.16°** | 31° (1myw crystal) | 7° geometry gap |
| kappa vs their 2D form −(1+cos²δ) | ratio 0.954 ± 0.015 | — | agrees to 4.7% |
| gamma_0 (monomer) | **38.7–44.4°** | 22° (from Ω = 0.709, eq 7) | disagrees, ~2× |
| 1PA shift, excitonic | **−6.16 ± 1.68 cm⁻¹** | ≈ −23 cm⁻¹ | ours matches our own TD−TDX (−6.1) |
| J (screened, ε=1.77) | 32.82 ± 1.55 cm⁻¹ | 18.1–22.6 screened | high by 1.61× |
| J (unscreened) | 58.09 cm⁻¹ | 32–40 unscreened | high by 1.61× |
| implied \|mu\| | 9.81 D (STEOM) | **7.6 / 7.9 D** | theirs matches experiment, ours does not |

**Most of the delta gap is a difference in what the angle is measured *between*.**
Cusick estimate δ from a structural proxy — the OH→CB axis of the tyrosine
chromophore precursor — because AlphaFold3 "does not consider any
post-translational modifications" and never builds the mature chromophore (their
SI Note S4). We use the STEOM transition dipole. On the *identical* crystal
structure the two definitions give **38.16° vs 27.40°**, an 11° offset that lands
next to their published 31°; on our tandem the proxy gives 25.3° rather than
39.6°. The separation cross-checks exactly: our CB2–CB2 = 25.38 Å against their
Table S2 value of 25.4 Å. See
[`results/delta_definition_offset.md`](results/delta_definition_offset.md).

Two earlier candidate explanations were tested and are **not** the cause.
Point-dipole inference does bias their fitted angle (removing our f = 1.1872 ±
0.0114 moves 14–20° to 37–40°), but AlphaFold3 reaches ~15° with no spectroscopy
at all, so that cannot be it —
[`results/angle_inference_bias.md`](results/angle_inference_bias.md). And their
own Table S2 shows the five AlphaFold3 structures giving δ = 9, 9, 8, **15**, 9°
and couplings spanning 29–43 cm⁻¹; the main text quotes the single value that
overlaps their spectroscopic range, and the modal prediction of 9° does not
overlap it.

**What does survive is the transition dipole.** Inverting their eq 13 for |μ|
using their own J, R and δ gives **7.6–7.9 D** — the experimental value, and
independent corroboration of the weekend basis-convergence result that our STEOM
9.8 D is the outlier. Since J ∝ |μ|², that single input accounts for most of the
residual 1.61× coupling gap, and it is unaffected by the δ question.

**They also avoid the trap our own draft falls into.** Cusick observe their Ω
features separated by 750 cm⁻¹ against an H–J splitting of only 60–90 cm⁻¹,
identify the shape as the first derivative of a Gaussian with extrema at
FWHM/√(2 ln 2) = 2σ, and deliberately fit amplitudes rather than reading a
splitting off the separation — exactly the error the JPCB round-2 CD argument
makes.

## Three corrections to the record

1. **The screening bullseye is not real.** Our 32.82 is screened by ε = 1.77;
   Cusick explicitly do not screen. Put both in a single convention and we are
   high by **1.61×** either way — the factor cancels from the ratio. The
   apparent 1% agreement existed only because 1.61 happens to sit near 1.77.
   ε is now the dominant systematic, larger than the near-field correction
   originally reported.

2. **13.31 cm⁻¹ is not the point-dipole partner of 32.82.** The genuine partner
   is on disk: `J_pda_cm` over the same 1000 frames, mean **27.64 cm⁻¹**, ratio
   **1.19×**, not 5.6×. The 13.31 value belongs to a different structure
   (crystal `venus_dimer.pdb`), a different density (TDDFT cube) and a single
   static frame. At R ≈ 25 Å the point-dipole approximation is good to ~19%, so
   cheap estimates are adequate for FP dimers at biological separations.

3. **R = 27.6 Å in the submitted JPCL manuscript is not reproducible.** Neither
   is its companion 92.85° inter-dipole angle. The density-centroid explanation
   was tested and rejected (centroids differ by only 0.50 Å). See
   [`results/R_provenance.md`](results/R_provenance.md). No manuscript was edited.

## The CD couplet

Two drafts in this repo commit to incompatible readings of the couplet
separation. Neither survives. Feeding the round-2 Omega = 553.9 cm⁻¹ through the
actual lineshape gives a separation of **652 cm⁻¹**, not the 548 claimed — the
"parameter-free match" equates two quantities instead of computing one from the
other. But the separation is not pure linewidth either: it is a **non-monotonic,
two-valued, lineshape-dominated** function of the gap, so it is not a usable
estimator in either direction. The couplet *amplitude* (2|J|/Omega) is where the
coupling actually lives. See
[`results/cd_splitting_resolution.md`](results/cd_splitting_resolution.md).

## Scripts

| Script | Produces |
|---|---|
| `derived_parameters.py` | `results/derived_parameters.{json,csv}` — delta, kappa, 1PA shift over 1000 frames + static |
| `parse_orca_state_dipoles.py` | `results/gamma0.json` — gamma_0 from ORCA state dipoles |
| `screening_table.py` | `results/screening_reconciliation.{md,json}` |
| `cd_splitting_resolution.py` | `results/cd_splitting_resolution.{md,json}` |
| `angle_inference_bias.py` | `results/angle_inference_bias.{md,json}` — multipolar bias in their inferred delta (tested, insufficient) |
| `delta_definition_offset.py` | `results/delta_definition_offset.{md,json}` — transition dipole vs their structural axis proxy |

Run with any numpy-capable interpreter, e.g.

```bash
/home/robson/anaconda3/envs/adcc_env/bin/python cusick_comparison/derived_parameters.py
```

## Verification built in

- `derived_parameters.py` re-derives the stored `J_pda_cm`, `angle_deg` and
  `separation_A` columns from the raw npz vectors and **aborts** if they do not
  match (currently 2×10⁻⁷, 1.4×10⁻¹⁴, 0). This gates the au/Ångström convention
  and the A→B direction of `r_hat`, which the npz does not record.
- `cd_splitting_resolution.py` reproduces this repo's own reconstruction of the
  published Table-S3 fit (500.2 against 500.0 cm⁻¹) before drawing conclusions,
  and aborts otherwise.
- `parse_orca_state_dipoles.py` guards the root↔IROOT mapping by energy rather
  than trusting index order, and reports gamma_0 under all three EOM-CC
  transition-moment conventions so the spread is visible.

## Known limitations

- gamma_0 uses **unrelaxed** state dipoles; no relaxed state-dipole block exists
  anywhere in this tree. delta_mu is a difference of two large, nearly
  cancelling vectors, so this caveat is load-bearing.
- The CD analysis uses Lorentzian and Gaussian profiles; the real lineshape is
  Voigt-like and was not fitted.
- Cusick's delta ranges are read from a summary of their Table S2, not from the
  paper directly. The multipolar argument depends specifically on their *tandem*
  delta being spectroscopically inferred rather than structural — **verify that
  against the source before publishing**, because the explanation collapses if
  their tandem angle is also structural.

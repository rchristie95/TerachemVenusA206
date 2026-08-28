# Superseded results, 2026-08-27

One authoritative list. Files across the repo still assert the left-hand column
as current; each carries a dated pointer to this file. Nothing here is a
correction to arithmetic -- every entry is a change of interpretation or a
failed consistency check found on this date.

| superseded | replacement | established in |
|---|---|---|
| `<\|Delta\|>` = 576 cm^-1 treated as STATIC detuning | 576 is the INSTANTANEOUS value; tau_c = 8.5 fs against hbar/2J = 81 fs, so the exciton sees **189-324** | `tddft_ensemble_decoherence/audit/REPRODUCIBILITY_AND_TIMESCALE.md` |
| `J = 32.82 +/- 1.55` cm^-1 as settled | **21-26 cm^-1**, a range not a point | `dipole_gauge/GAUGE_VERDICT.md` |
| `\|mu\| = 9.81 D` from the STEOM density | STEOM's length/velocity gauges disagree 1.51x and **worsen** with basis; TDDFT is gauge-consistent at **8.68 D** | `dipole_gauge/GAUGE_VERDICT.md` |
| STEOM dipole "basis-converged to -0.54%" | that measured the stability of the length-gauge number, not the moment | same |
| homogeneous linewidth unresolved (300 / 611 / 41.5) | **HWHM 346 cm^-1** from the measured C(t), no fast/slow cut | `terachem_site_energy_cd/results/spectral_density_lineshape.json` |
| "the couplet width requires `<\|Delta\|>` <~ 115" | **dead** -- it fed the OBSERVED width in as the homogeneous input (circular). Width predicts 342.7 vs observed 366.4, parameter-free | `terachem_site_energy_cd/results/CD_WIDTH_RESOLVED.md` |
| couplet width as evidence for the detuning mechanism | the width measures the LINEWIDTH, not the disorder; it is **silent** on the mechanism in both directions | same |
| crystal register reproduces 57 ps "exactly, no free parameters" | that used the instantaneous detuning. At the static value it **overshoots at +5 sigma** | `exciton_observables/STATIC_DETUNING_REVERSAL.md` |
| "it is the MD relaxation that breaks the agreement" (`exciton_observables/README.md`, Result 4) | **WITHDRAWN.** The 1 ns ensemble was not converged; 67.7 ns settles at alpha = 111.19 deg, on the crystal interface | same |
| converged MD as independent support for the crystal register | **it is not.** The tandem was BUILT in that register and cannot re-dock in 67.7 ns, so it was never able to disagree | `tandem-construct-is-truncated` (memory) |
| J ~= 30 cm^-1 working hypothesis | dead; killer test gave `<\|Delta\|>` = 463.5 against the >=1000 required | `exciton_observables/HYPOTHESIS.md` header |
| AlphaFold3 delta = 15 deg as structural evidence | **no AlphaFold variant determines this register**; alpha spans 55 deg across methods, each internally reproducible | `alphafold_tandem/ALPHAFOLD_REGISTER.md` |
| AlphaFold3 weights are gated behind a request | freely downloadable since before 2026-08; AF3 also accepts precomputed MSAs, so it needs no genetic database | same |
| our 495-aa tandem = the experimental construct | the deposit (Addgene #250938) is **511 aa**; ours drops LGMDELYK from both barrels and carries a wrong linker Cys (fixed in `tandem_dimer_v2.pdb`) | `alphafold_tandem/pdVenus-TD_Addgene250938_511aa.fasta` |
| tandem MD boxes are neutral | **net -2 e**: CR2 is correctly anionic (-1.0 each) but neutralisation supplied only +11 against a true -13 | `cr2_amber_parameter_transplant.json` |

## Two things that did NOT change

- The **site energy** (527-532 nm) and the STEOM excitation energy are unaffected;
  only the transition MOMENT fails its consistency check.
- The dimer remains **detuning-dominated and localised** -- by 3-5x rather than
  9x, but the eigenstates are ~98% localised either way.

## Open, and now the largest discrepancies

1. **Lobe asymmetry 1.30 computed vs 2.34 observed.** Needs a single-signed
   component at 40% of the couplet amplitude. A conservative couplet is a
   theorem of the two-state model, so this requires a term the model lacks --
   most likely the monomer's intrinsic induced CD, which enters
   `Delta[theta] = TDX - TD` at full weight.
2. **The anisotropy/superradiance impossibility.** `cos^2 x <= |cos x|` makes
   the two experiments contradictory for ANY geometry. Candidate resolution: they
   do not measure the same angle -- superradiance uses the inter-transition-dipole
   angle, depolarisation uses absorption-on-one vs emission-on-the-other, which
   differ by the intra-chromophore relaxation angle `beta`. See
   `NEW_MODEL_PROPOSAL.md`.

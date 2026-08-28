# Can AlphaFold determine the dVenus-TD barrel register?

2026-08-27. Supersedes `AF3_REGISTER_VERDICT.md`, `SEED_SWEEP_RESULT.md` and
`TRUNK_PROBE_RESULT.md`, which were three passes at one question.
Raw output in `af3_server_run/` and `af3_seed_sweep/`; local runs in
`colabfold/`.

**Answer: no. And neither can anyone else's AlphaFold run, including the one in
Cusick et al.**

## 1. The barrels are excellent; the register is not

AlphaFold Server, 511-aa construct (Addgene #250938):

| quantity | value |
|---|---|
| pLDDT per barrel / fold RMSD vs crystal | 90.2 / **0.70 A** |
| linker 240-272 pLDDT | **30.7, 100% below 70** |
| PAE within a barrel | 2.7-2.9 A |
| **PAE between barrels** | **18.7-21.8 A, 0.0% of pairs below 5 A** |
| pTM | 0.57 |

PAE saturates near 31.75 A. AF3 folds each barrel to 0.70 A and then reports
that it does not know how they sit relative to each other.

Superposing barrel 1 on crystal chain A and carrying barrel 2 through the same
transform gives **11.7 A** against 0.70 A for the fitted barrel -- not the
crystal register, but by its own PAE not any *particular* register.

Contact counts do NOT discriminate: the crystal dimer interface is itself small
(35 CA pairs within 8 A, centroids 29.6 A) and AF3's is comparable (18-44,
30.0-30.5 A). Venus is a weak dimer; do not cite interface size as evidence.

## 2. We DO reproduce Cusick's delta -- their formula is not the inter-dipole angle

Note S4 defines delta as the angle between `r_a` and `r_a - r_b` (r = OH->CB of
the tyrosine precursor), which for unit dipole directions is exactly
**delta = 90 - theta/2**.

| source | delta | R |
|---|---|---|
| our AF3, 14 seeds | **13.5-16.1** | 25.6 A |
| Cusick's AF3 | **8, 9, 9, 9, 15** | 24.3-27.5 A |
| crystal, ours / theirs | 27.4 / 31 | 25.4 A |

Both delta and R reproduce, across different constructs (their tyrosines are
116/388, ours 67/331). An earlier version of this file claimed we got 27 against
their 15; that compared our inter-dipole angle to their delta.

## 3. The SEED cannot probe the register -- do not repeat this experiment

14 distinct seeds (7 x 495-aa, 7 x 511-aa, 70 structures). One-way ANOVA on
delta grouped by seed: **F = 2.29 p = 0.064** (495) and **F = 0.30 p = 0.930**
(511). Seeds are statistically indistinguishable; the between-seed spread of
means is no larger than the 5 diffusion samples alone give.

The register is fixed upstream of the seed (sequence -> MSA -> trunk); the seed
samples only diffusion noise. To perturb the register you must perturb the
trunk.

Also established: **the construct truncation is irrelevant to the register.**
495-aa vs 511-aa moves delta 0.9 deg, alpha 1.2 deg, register RMSD not at all.

## 4. Perturbing the trunk: every method disagrees

| method | alpha | register RMSD vs crystal |
|---|---|---|
| **crystal 1MYW** | **110.4** | **0** |
| AF3, 14 seeds | 119-120 +/- 3 | 11.7 |
| AF2 tandem, full MSA | **165.1 +/- 4.5** | **35.0** |
| AF2 tandem, MSA 64:128 | 127.8 +/- 7.7 | 19.2 |
| AF2 tandem, MSA 32:64 | 129.1 +/- 9.4 | 19.8 |
| AF2 tandem, MSA 16:32 | 128.7 +/- 20.3 | 24.5 |
| AF2-multimer, free dimer | 134.3 +/- 18.2 | 21.4 (ipTM 0.21-0.46, failed) |

**alpha spans 55 degrees** while each method reproduces itself to a few degrees.
Internal reproducibility is not accuracy. All are genuinely DOCKED (centroid
separations 31-32 A vs the crystal's 29.6), so AF2 is not failing to fold -- it
is confidently docking elsewhere. Thinning the MSA moves alpha 37 deg, so the
full-MSA answer is one basin a deep alignment collapses onto.

## 5. Consequence

**No AlphaFold result is evidence about this register, in either direction** --
including Cusick's delta = 15. Their AF3-vs-spectroscopy agreement is not the
independent structural confirmation it reads as.

It does **not** support the crystal register either. It removes structure
prediction from the argument and leaves the spectroscopy to decide.

Caveat: AF2 monomer mode on a two-domain fusion is not its strength and
AF2-multimer scored ipTM < 0.5, an outright failure. Some of the spread is
method weakness. But we cannot tell which method is the weak one, which is
exactly why none of them settles it.

## Reproduce

```bash
python3 alphafold_tandem/analyse_af.py                     # server runs
bash    alphafold_tandem/colabfold/run_probe.sh            # local AF2 probe
```

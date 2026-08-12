# Explaining the Nguyen/Kim experiments without a coherent exciton

2026-08-12. Literature deep dive plus direct calculation. Working assumption
J ~ 30 cm^-1; computed detuning |Delta| ~ 455-576 cm^-1 on two independent
registers, so Delta/2J ~ 8, and the dimer is ~99% localized.

## The headline: the "Davydov splitting" is Omega, not 2J

The exciton gap of a **detuned** dimer is `Omega = sqrt(Delta^2 + 4J^2)`, which
for this system is dominated by the detuning:

| \|Delta\| (cm^-1) | Omega | at 516 nm | couplet amplitude 2J/Omega |
|---:|---:|---:|---:|
| 0 | 65.6 | 1.75 nm | 1.000 |
| 455 | 459.7 | 12.24 nm | 0.143 |
| **550** | **553.9** | **14.75 nm** | **0.119** |
| 576 | 579.7 | 15.44 nm | 0.113 |

**Kim et al. measured 14.6 +/- 0.3 nm. Our computed J and Delta give 14.75 nm,
inside their error bar, with no free parameters.** Reading the same number as
2J gives J = 274 cm^-1, which at 25 A separation needs a transition dipole of
~29 D, or a separation of ~10 A (bacteriochlorophyll special-pair distance) —
impossible for two beta-barrels of diameter >= 22 A.

Corollary: the couplet **amplitude** (suppressed 8x by 2J/Omega = 0.119) is the
observable carrying J; the couplet **width** is a linewidth observable that
saturates at the bandwidth. An independent simulation reproduces a 14.7 nm
couplet from J = 33 cm^-1, and J = 274 gives 19 nm — a *worse* fit to Kim.

## Observable by observable

| observable | Nguyen/Kim reading | what it actually measures | status |
|---|---|---|---|
| CD "Davydov splitting" 14.6 nm | 2J -> J = 274 | Omega = sqrt(Delta^2+4J^2) = 14.75 nm | **explained** |
| "apparent coupling" 102-357 cm^-1 | J | half a couplet separation, i.e. a linewidth; our vibronic ensemble gives an apparent delta of 408 cm^-1 from a true J of 30 (13.6x) | **explained** |
| absorption red shift 34 cm^-1 | excitonic | J cos(alpha), tens of cm^-1 — consistent with J ~ 30, not with 274 | **explained** |
| lifetime -57 +/- 4 ps | superradiance | see below — at least three non-excitonic mechanisms of exactly this size | **explained** |
| anisotropy 0.52 -> 0.30 | alpha = 48.7 deg | genuine geometry, but with ~+/-8 deg systematics | **partly open** |
| antibunching at 20 nm | coupling | J(200 A) ~ 0.06 cm^-1; cannot be excitonic | **not coupling** |

## The lifetime shortening does not need coherence

- **Already published as a refractive-index effect for the same architecture.**
  Teijeiro-Gonzalez 2021 (Biophys J 120:254): eGFP tandem dimer lifetime
  "slightly and consistently lower" than the monomer (~0.8%), attributed to the
  neighbouring barrel raising the local refractive index, explicitly not
  excitons. van Manen 2008 (Biophys J 94:L67) gives -1.8% lifetime per 0.01
  refractive-index unit, so the observed 1.88% needs only **Delta n ~ 0.010**.
- **The same consortium published the null result.** Sanchez-Pedreno Jimenez,
  Puhl, Vogel & Kim 2023 (PCCP 25:19532): dEGFP-TD 2.61 ns vs mEGFP 2.64 ns
  (1.1%), concluding dimerisation "does not affect ... the average time the
  fluorophore remains in its excited state."
- **A 1.9% dark-chromophore fraction reproduces 57 ps exactly**, and sits inside
  Nguyen's own absorption bound (<=2.7%) and FCS brightness bound (2.0 +/- 0.1).
- **Superradiance has never been shown in an FP.** Huff 2021 (JPCB 125:10240):
  Cy5 dimers at J = 387-589 cm^-1 show lifetime reductions *"inconsistent with
  ... superradiance"*; nonradiative decay dominates.

And where it *is* read as superradiance, it favours our coupling. Predicted
Dtau at our Delta:

| \|cos alpha\| | J=30 | J=33 | J=131 | J=186 |
|---:|---:|---:|---:|---:|
| 0.213 (MD) | 36 | 39 | 146 | 194 |
| 0.349 (crystal) | 58 | 64 | 239 | 318 |

**Measured 57 +/- 4.** Their own CD-derived couplings overpredict by 3-6x.

## The anisotropy is the one real structural observable — and it is softer than claimed

It is independently confirmed: **Kinoshita 2026 (PCCP 28:14488) measure 49 deg
for the eYFP dimer by ONE-photon anisotropy**, agreeing with Nguyen's 48.7 deg.
So the two-photon tensor is not the problem — and for EGFP the tensor is
measured to be single-element and collinear with the S0->S1 dipole
(Masters 2018, JCP 148:134311).

But the central value carries systematics Nguyen do not quote:

| variant | alpha | obtuse branch |
|---|---:|---:|
| as published (R_TDX = 0.52) | 48.7 | 131.3 |
| using the dTomato reference 0.571 | 52.7 | 127.3 |
| 90% paired/mature | 52.3 | 127.7 |
| 80% paired/mature | 57.1 | 122.9 |
| **Kim 2019, same construct, 950 nm** | **56.7** | **123.3** |

Between 2019 and 2025 the monomer reference moved 0.42 -> 0.52 and the tandem
0.20 -> 0.30 — both by exactly +0.10, across a change of excitation wavelength
and IRF. That is a calibration shift, and since the *ratio* sets the angle it
moves the answer from 56.7 to 48.7 deg. A defensible value is **49 +/- 8 deg,
or its 131 +/- 8 deg mirror**; the superradiance sign selects the obtuse
branch. Our computed geometry is 93-110 deg — still outside, but the gap is
~15-20 deg, not the 30 deg originally supposed.

Incomplete pairing is the largest single lever and is the same group's own
finding: Jimenez 2023 report the dEGFP tandem limiting anisotropy rising from
0.30 to 0.32-0.37 with glycerol and attribute the amplitude to "the fraction of
dEGFP molecules that are dimerised in tandem dimers."

## The self-contradiction at the centre of the claim

Nguyen frame the result as exciton **delocalization** but analyse the
anisotropy with the **incoherent-hopping** formula
`r_TD = (r_TDX/2)[1 + d2(alpha)]`. In the delocalized limit the two exciton
dipoles are orthogonal, emission from the same exciton preserves r0 entirely,
and the formula does not apply. **The observed anisotropy drop is therefore
evidence for localization plus hopping** — which is exactly what
Delta/2J ~ 8 implies, and what the measured 4.4-20 ps homo-FRET times show
(Jung 2005 ChemPhysChem; Kinoshita 2026).

## Concentration: the two pillars are different species

**eYFP Kd = 36 +/- 4 uM** (Kinoshita 2026, AUC; 20 uM by another route).
Nguyen ran CD at **25 uM** — below the Kd, so a large intermolecular-dimer
fraction — while themselves recommending "~0.75 uM or lower", i.e. 33x over
their own limit. Kim's couplet appears only at high concentration (at 0.98 uM,
A206 and A206K are indistinguishable). Antibunching was done at 200-500 nM,
where nothing is dimerised. **CD and antibunching are measured four to six
orders of magnitude apart in concentration and are not the same species.**

## The consortium has already moved

**Cusick, ..., Kim, Vogel & Drobizhev 2026, JPCA 130:5471** — "...Detects
**Weak** Excitonic Coupling in Fluorescent Protein Dimers." Same authors.
Paywalled; **getting this full text is the highest-value follow-up.**
There is no independent replication of any of this work, and no published
critique, comment or correction of either paper.

## What would settle the remainder

1. **Relative quantum yield of TD vs TDX to +/-1%.** Superradiance predicts
   k_r +3.4% and Phi +1.4%; quenching predicts k_r unchanged and Phi -1.9%.
   Nobody has measured it. This is the decisive experiment.
2. **CD at <=1 uM**, below the dimerisation Kd.
3. **Masters' three-observable protocol on Venus** (linear + circular 2PA plus
   R_L(0), R_C(0)) to pin the tensor and the emission angle independently.

# Comparing with a two-photon experiment without doing two-photon spectroscopy

Supporting note for the Cusick et al. 2026 comparison. Derivations and
conventions; the numbers live in `cusick_comparison/results/`.

## Why this is legitimate

Omega — the two-photon polarization ratio — is Cusick's *instrument*, not their
result. What emerges from their analysis is a set of ordinary molecular
parameters: an interchromophore angle, an orientation factor, a monomer angle
between two dipoles, a coupling, a separation. Every one of those is something a
transition-density pipeline produces natively. Comparing at the level of derived
parameters is standard practice between a quantum-chemistry paper and a
spectroscopy paper, and it requires no 2PA measurement on our side.

It has a second advantage. Three of the four quantities below are dimensionless
angles. They are therefore immune to the dielectric-screening convention that
makes the raw coupling comparison ambiguous (see
`results/screening_reconciliation.md`), which makes them the *cleanest* tests
available rather than a fallback.

## 1. delta from the interdipole angle

Cusick define a dimer frame by

    x_hat  ||  mu_a - mu_b ,    y_hat  ||  mu_a + mu_b

with delta the angle between each monomer transition dipole and x_hat. Because
the two monomers are related by the dimer's near-C2 symmetry, mu_a and mu_b sit
symmetrically about x_hat, each at angle delta, so the angle between them is

    theta_ab = 180 - 2 delta

Equivalently, in the form that is numerically safer because it avoids
constructing the frame at all,

    cos theta_ab = -cos 2 delta        =>        delta = 90 - theta_ab / 2

theta_ab is the plain angle between the two transition-dipole vectors, which
this repo already computes identically in six places (canonical form at
`coupling_dcd_steom.py:226-227`). No frame construction, no alignment, no new
convention.

**Phase.** The overall sign of a transition dipole is arbitrary: mu_a -> -mu_a
is the same physical state. That flips cos theta_ab and hence maps delta to
90 - delta. We fold delta into [0, 90] and say so. The folding is explicit in
`derived_parameters.py`, not a silent `abs()`.

## 2. kappa, and the test of their quasi-2D assumption

Their eq 13,

    nu_bar = mu^2 (1 + cos^2 delta) / (h c R^3)

does not follow from the geometry alone. It requires R_hat || x_hat — the
centre-to-centre vector lying along mu_a - mu_b. That is an assumption of their
two-dimensional model, not a measurement, and they acknowledge only that their
2D and 3D estimates "almost coincide" for the vdW dimer.

We have the real three-dimensional structure, so we can check it directly. The
general Kasha orientation factor is

    kappa = mu_a_hat . mu_b_hat - 3 (mu_a_hat . R_hat)(mu_b_hat . R_hat)

and their assumption is the special case kappa = -(1 + cos^2 delta). Evaluating
both *at the same delta* is a self-consistency test of their model that only a
3D structure can perform. If kappa departs materially from the 2D form, their
nu_bar shifts with it, and because their final answer comes from intersecting
eq 15 with the 1PA-shift constraint, that error would propagate into their
quoted coupling.

Note that `jdd` as computed in this repo uses unnormalised dipoles, so the
normalised factor is `kappa = jdd / (|mu_a| |mu_b|)`. The reusable batched form
is `point_dipole_coupling_cm` at
`terachem_site_energy_cd/scripts/absolute_cd_vs_experiment.py:111`.

## 3. gamma_0, the monomer angle

gamma_0 is the angle between the transition dipole mu and the difference dipole

    delta_mu = mu_excited - mu_ground

This needs no dimer, no MD and no coupling — one excited state of one monomer
settles it, which makes it the cheapest possible validation of the electronic
structure underpinning the whole transition density.

Two practical traps, both live in our ORCA outputs and both handled in
`parse_orca_state_dipoles.py`:

- A STEOM output contains **four** blocks headed `ABSORPTION SPECTRUM VIA
  TRANSITION ELECTRIC DIPOLE MOMENTS`. The first is the CIS/TDDFT guess at a
  completely different energy (331.8 nm against 523.9 nm). The remaining three
  are the RIGHT, LEFT and LEFT-RIGHT conventions of the non-Hermitian
  coupled-cluster response. Every block must be keyed off the
  `SPECTRUM FOR <X> TRANSITION MOMENTS` header that precedes it.
- ORCA's `FINAL STEOM-CCSD ABSORPTION SPECTRUM` block is internally mixed: its
  `fosc`/`D2` are the LEFT-RIGHT geometric mean while its printed `DX/DY/DZ` are
  the RIGHT vector. That is why |mu| from the components (3.8810 au) disagrees
  with sqrt(D2) (3.9070 au) by 0.7%. This is the L/R asymmetry of EOM-CC
  transition moments, not rounding.

**Caveat that matters.** ORCA prints *unrelaxed* excited-state dipoles, and no
relaxed state-dipole block exists anywhere in this tree. delta_mu is a
difference of two large, nearly cancelling vectors, so orbital relaxation
affects it disproportionately. A disagreement with Cusick's 22 deg is a finding
about this level of theory; it is not a settled refutation of either side.

## 4. The 1PA shift

Kasha's shift for the dimer relative to the monomer,

    delta_nu = J (tan^2 delta - 1) / (tan^2 delta + 1)

simplifies, since cos 2d = (1 - tan^2 d)/(1 + tan^2 d), to

    delta_nu = -J cos 2 delta = J cos theta_ab

The second form is worth noting: the red shift is just J cos(interdipole angle),
which is the same relation recorded independently elsewhere in this project. It
is also an ordinary one-photon absorption observable, measured by Cusick at
-35.3 +/- 0.4 cm^-1 with at most 12 cm^-1 of Stark contribution, leaving about
-23 cm^-1 excitonic. So it tests J and delta jointly against a measurement that
involves no two-photon physics at all.

## What the comparison is not

It is not a validation exercise looking for agreement. Three of the four
comparisons disagree with Cusick (delta, gamma_0, and the 1PA shift), and one
agrees well (kappa against their 2D assumption). That pattern is informative:
their *model* is sound where we can test it, and the disagreement is
concentrated in the *geometry*. Reporting it the other way round — or reporting
only the agreement — would be the failure mode to avoid.

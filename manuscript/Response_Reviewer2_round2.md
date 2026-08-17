# Response to Reviewer 2

We thank the reviewer, and in particular for performing an independent Q-Chem
coupling calculation. That calculation is correct, our published ratio was not,
and tracking down why has substantially improved the paper.

---

## Major point: the dipole-dipole approximation does not underestimate the coupling by a factor of four

**The reviewer is right, and we withdraw the claim.** The factor of $\sim$5.6
reported in the original manuscript was an error in unit conversion, not
physics.

The Coulomb double sum was evaluated with grid coordinates in Angstrom. The
reciprocal distance must then be converted with
$1\,\mathrm{\AA}^{-1} = 0.529177\,a_0^{-1}$. The original driver applied no
conversion at all, which inflates every transition-density coupling by
$1/0.529177 = 1.8897$ while leaving the point-dipole values — evaluated in
bohr throughout — untouched. The published ratio was inflated by exactly that
factor. (A later production script applied $1.8897$ *instead of* $0.529177$,
inflating by $1.8897^2 = 3.5711$; the two errors are different and we note it
so that no reader rescales the older numbers by the wrong constant.)

**Our corrected numbers reproduce the reviewer's calculation.** We had already
run a method-matched gas-phase control on the crystal geometry in ORCA 6.1.1.
Against the reviewer's Q-Chem values:

| | Reviewer (Q-Chem) | This work (ORCA) | agreement |
|---|---:|---:|---:|
| TDC | 6.02 meV (48.7 cm$^{-1}$) | 6.087 meV (49.10 cm$^{-1}$) | 1.1% |
| PDA | 5.36 meV (43.2 cm$^{-1}$) | 5.361 meV (43.24 cm$^{-1}$) | 0.01% |
| TDC/PDA | 1.123 | 1.135 | — |

The point-dipole terms agree to three significant figures and the distributed
densities to 1%. There is no remaining disagreement.

Across the 1000-snapshot production ensemble the corrected ratio is
**TDC/PDA = 1.187 $\pm$ 0.011**, i.e. the near-field correction is a ~19%
effect, exactly as the reviewer's chemical intuition for well-separated
chromophores requires. A multipole decomposition of the same transition density
accounts for it term by term: dipole–dipole 24.81, dipole–quadrupole 4.03,
dipole–octupole 0.73, quadrupole–quadrupole 0.16 cm$^{-1}$, summing to 29.74
against the full 30.45 cm$^{-1}$.

The corrected ensemble coupling is $J = 32.82 \pm 1.55$ cm$^{-1}$. We have added
a regression test that checks a closed-form two-charge reference and greps the
source for the wrong constant, so this class of error cannot recur.

**Consequences, which we now make the centrepiece of the paper.** The corrected
coupling is 4–8$\times$ smaller than the values inferred from CD couplet
splittings. The revised Results show that this is not a conflict, but for a
different reason than we gave in the previous version, and we are grateful to
have been pushed to check it numerically rather than argue it analytically.

We had proposed that the couplet separation reports the detuned exciton gap
$\Omega = \sqrt{\Delta^2 + 4J^2}$, and that the computed $J$ and $\Delta$
reproduce the measured $14.6 \pm 0.3$ nm with no adjustable parameter. **That
argument does not survive simulation and has been withdrawn.** $\Omega$ is a
Hamiltonian gap; the measurement is the peak-to-peak separation of a bisignate
lineshape, and for bands this broad the two are not equal. Propagating
$\Omega = 553.9$ cm$^{-1}$ through the lineshape of the published constrained
fit gives a couplet separated by **652 cm$^{-1}$, not the observed 548**.

What replaces it is stronger, and does not rely on a numerical coincidence. On
that same lineshape the peak-to-peak separation is a **non-monotonic,
two-valued and linewidth-dominated** function of the gap: the published latent
gap of 261.58 cm$^{-1}$ yields extrema 500 cm$^{-1}$ apart (a 1.93$\times$
inflation), a gap of *zero* still yields 690 cm$^{-1}$, and the function passes
through a minimum near $\Omega \approx 190$ cm$^{-1}$. The couplet separation is
therefore not an estimator of the exciton gap in either direction, and the
factor of four never required a physical explanation. The detuning remains
real and consequential — it localises the eigenstates and suppresses the
couplet *amplitude* by $2|J|/\Omega = 0.12$ — but it is the amplitude, not the
separation, that carries $J$. The same artefact is visible in the two-photon
data of Cusick et al., whose $\Omega$ spectrum shows extrema 750 cm$^{-1}$
apart against an H–J splitting of only 60–90 cm$^{-1}$, and who correctly
analyse it by fitting amplitudes rather than separations.

We also note that Cusick et al. (*J. Phys. Chem. A* **2026**, 130, 5471), from
the laboratories that reported the larger couplings, now measure 32–40
cm$^{-1}$ for this dimer. **We have corrected our comparison to this work.**
Their values carry no dielectric screening, which they omit explicitly, whereas
ours are screened by $\varepsilon_{\rm opt} = 1.77$; the previous version
compared the two directly and reported agreement that the conventions do not
support. Like for like we obtain 58 cm$^{-1}$ unscreened against their 32–40,
or 32.8 against 18–23 screened — a ratio of 1.6 in either convention, since the
screening factor cancels from it. Inverting their point-dipole expression for
the transition dipole implied by their own parameters returns 7.6–7.9 D against
the 9.8 D of our spectroscopically normalised density; as $J \propto |\mu|^2$
this accounts for a factor of 1.5 and is now identified in the Discussion as
the dominant remaining systematic in the calculation.

---

## Minor point 1: transition atomic charges vs transition density cubes

Correct, and the wording was wrong. We use the transition-density-cube approach
of Krueger, Scholes and Fleming; $q_i$ and $q_j$ are transition densities
integrated over individual **grid voxels** of the cap-masked cube, not
atom-centred transition charges. No atomic-charge fitting is performed at any
stage. The Methods now state this explicitly and cite Krueger et al., and
distinguish it from the TrEsp/Renger atomic-charge scheme.

## Minor point 2: the 44-atom QM region

Apologies — the count was underspecified. The region is:

- the complete anionic CR2 chromophore residue, **29 atoms**;
- the $\pi$-stacked Tyr203 **phenol group only**, **12 atoms** (CG, CD1, CD2,
  CE1, CE2, CZ, OH + HD1, HD2, HE1, HE2, HH), cut at the C$_\gamma$–C$_\beta$
  bond, with the C$_\beta$ and backbone excluded;
- **3 link hydrogens** capping the severed boundary bonds.

Total 44, stoichiometry C$_{19}$H$_{18}$N$_3$O$_4$, net charge $-1$. A count of
45 arises if the tyrosine C$_\beta$ is retained rather than replaced by a link
hydrogen. The Methods now give this breakdown, and a fully labelled structure
diagram has been added as Fig. S1.

## Minor point 3: gas-phase controls in Table 1

Added, as a new table. Both methods were re-run on the identical frozen 44-atom
anion with every point charge removed:

| Method | Gas phase | Embedded | Shift |
|---|---|---|---|
| TDA-wB97X-D3/6-311G** | 27 324.5 cm$^{-1}$ (366.0 nm) | 27 710.2 cm$^{-1}$ (360.9 nm) | +385.7 cm$^{-1}$ |
| DLPNO-STEOM-CCSD/def2-SVPD | *[in progress]* | 19 088.2 cm$^{-1}$ (523.9 nm) | *[in progress]* |

The bright state is the same state in both environments (gas-phase and embedded
transition dipoles collinear to $\cos\theta = 0.9992$), so the shift is a
genuine electrostatic effect on a fixed state rather than a reordering. We note
that with the geometry frozen this isolates electrostatics only, and does not
include the geometric relaxation the protein also imposes.

## Minor point 4: sensitivity to the time-dependent screening

A fair question, and the answer is that it is immaterial to every spectroscopic
conclusion. Over the sub-picosecond window, $1/\varepsilon$ relaxing from 0.56
to 0.50:

| quantity | $1/\varepsilon = 0.56$ | $1/\varepsilon = 0.50$ | change |
|---|---:|---:|---:|
| $J$ | 32.53 cm$^{-1}$ | 29.04 cm$^{-1}$ | $-10.7\%$ |
| Rabi period $h/2\lvert J\rvert$ | 0.513 ps | 0.574 ps | $+12\%$ |
| $\Omega$ | 553.83 cm$^{-1}$ | 553.06 cm$^{-1}$ | $-0.14\%$ |
| minor-site weight | 0.69% | 0.55% | — |

The coupling itself moves by 11%, but the exciton gap by 0.14%, because
$|\Delta| \approx 550 \gg 2\lvert J\rvert \approx 65$ cm$^{-1}$: the detuning
sets the gap, the mixing and the exciton composition, essentially independently
of $J$. A 10% drift in $J$ during the first picosecond cannot turn a 99%
localised system into a delocalised one. This is now stated in Limitations.

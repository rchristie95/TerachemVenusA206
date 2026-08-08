# Prompt for the Linux GPU agent: frame-resolved TeraChem site energies and CD

Work in the `TerachemVenusA206` repository on branch
`codex/round-2-coupling-correction`. Pull that branch from `origin` and record the
starting commit. Use the licensed TeraChem installation and the available NVIDIA
GPUs. Do not edit `manuscript/JPCB_tandem.tex`. Do not edit
`manuscript/JPCB_tandem_round_2.tex` or make scientific claims in any manuscript
until the calculations and validation below are complete. Do not mention private
editorial correspondence in reader-facing text or generated scientific reports.

## Scientific question

Determine whether instantaneous inequivalence of the two Venus chromophores can
help explain the experimentally observed apparent circular-dichroism (CD) peak
separation. For every molecular-dynamics frame, use the nondegenerate two-site
Hamiltonian

\[
H_i=\begin{pmatrix}E_{1,i}&J_i\\J_i&E_{2,i}\end{pmatrix},\qquad
\Omega_i=\sqrt{(E_{1,i}-E_{2,i})^2+4J_i^2}.
\]

The audited trajectory result is
`J = 32.8165 +/- 1.5540 cm^-1`, so the resonant splitting is only
`2J = 65.6330 cm^-1`. The two published fits imply apparent separations of about
`262` and `372 cm^-1`; at the mean calculated coupling these would require
`|E1-E2|` of approximately `254` and `366 cm^-1`, respectively. Large detuning
also reduces excitonic mixing by
`2|J|/sqrt((E1-E2)^2+4J^2)`, so success requires reproducing the CD intensity and
line shape, not just a peak separation.

The experimental source is Nguyen et al., *Biophysical Journal* 124 (2025)
4293-4309, DOI `10.1016/j.bpj.2025.10.022`. Treat its fitted coupling as
"apparent": its spectral components may contain unresolved electronic or
vibronic contributions.

## First inspect and preserve the existing provenance

Read these files before changing code:

- `ORCA_ONLY_PIPELINE.md`
- `README.md`
- `qmmm_tddft_pipeline.py`
- `neo_model/build_qm_phenol_relaxed.py`
- `tc_simple_anionic/qm_selection.json`
- `tc_qmmm_opt_constrained/opt.in`
- `coupling_nvt_production_cr2_1000_20260721/README.md`
- `coupling_nvt_production_cr2_1000_20260721/coupling_samples.csv`
- `coupling_nvt_production_cr2_1000_20260721/coupling_geometry.npz`
- `absorption_cd_spectra.py`
- `reference/SHA256SUMS`

Locate the external 1 ns NVT trajectory and its matching topology on the Linux
machine. Validate their frame count, atom ordering, timestep, and SHA-256 values
against the repository records. The 1000 saved frames are 1 ps apart. They can
describe a static inhomogeneous ensemble, but they are too coarsely sampled to
derive a femtosecond energy-gap autocorrelation function or a homogeneous
linewidth.

Record `git rev-parse HEAD`, `terachem` path/version, GPU model(s), driver,
TeraChem license availability, Python environment, trajectory/topology paths,
checksums, and every numerical setting in a machine-readable run manifest. Do
not commit licensed binaries, raw trajectories, scratch files, orbitals, or
large cube files.

## Efficient calculation strategy

Do not begin with all 1000 frames.

1. Run a smoke test on one central frame for both labelled chromophores.
2. Run a deterministic stratified pilot of 20 frames spanning the trajectory
   and the observed low/median/high `J` and separation values.
3. If atom mapping, charges, state identity, and restart behaviour pass, expand
   to 50-100 stratified frames. Estimate convergence with bootstrap confidence
   intervals.
4. Only propose the full 1000-frame calculation if the pilot changes the
   scientific conclusion or the spectra remain unconverged.

Use one independent TeraChem process per GPU, assign it explicitly with
`CUDA_VISIBLE_DEVICES`, and write each `(frame, site)` job into an isolated
directory. Avoid GPU oversubscription. Extract all requested geometries and
point-charge files once, then run TDDFT jobs from those immutable inputs. A
failed job must be restartable without overwriting completed results. If a
validated TeraChem orbital restart substantially reduces SCF time, reuse a
nearby frame from the same site, but retain periodic cold-start checks to show
that restarting does not change the tracked state.

## QM/MM model

Perform vertical single-point calculations on the NVT coordinates; do not
geometry-optimize individual frames.

Build two strictly equivalent site calculations per frame. In the site-1 job,
chromophore 1 is QM and chromophore 2 remains part of the MM environment; reverse
the roles in the site-2 job. Use identical atom-selection, capping, charge,
basis, embedding, and state-selection rules for the two sites.

For the fast production pilot, first construct the existing 44-atom
CR2-plus-Tyr203 capped model independently at each site, following
`neo_model/build_qm_phenol_relaxed.py` and the audited atom order. Confirm rather
than assume its net charge and multiplicity; the retained correlated model is
charge `-1`, multiplicity `1`, whereas the older 274-atom TeraChem setup has a
different total charge because it contains additional residues. Never mix the
charge or atom count of these two QM definitions.

Validate the compact model on at least 5 representative frames against the
larger retained TeraChem QM region, if computationally feasible. Report whether
the compact model preserves the distribution of the site-energy difference
`Delta = E1-E2`; a common absolute energy shift is less important than a biased
or compressed detuning distribution.

Use electrostatic embedding that preserves the instantaneous, asymmetric
protein/solvent environment. Keep the MM atom identities and charge convention
consistent across frames so that moving cutoff membership does not create
artificial energy jumps. Prefer a residue-complete fixed selection or a tested
smooth alternative. Explicitly test the sensitivity to embedding radius and to
the treatment of solvent. Do not silently combine explicit solvent point
charges with an equilibrium dielectric in a way that double-counts solvent
response. Include at least one clearly defined point-charge/no-PCM calculation;
an optical dielectric (`epsilon` about 1.77) may be a labelled sensitivity test
if supported by the installed TeraChem version. Do not use `epsilon = 78.39` as
the sole vertical-excitation model without demonstrating why it is appropriate.

## TDDFT/TDA calculations and state tracking

Begin with the repository method, `wB97X-D3`, using TDA and at least 5 roots.
Use `6-31G*` for a timed smoke/pilot if necessary, then check a representative
subset with `6-311G**`. Preserve the comparison of detunings even if a common
absolute shift is applied later. Document every deviation from the retained
input settings.

Do not select the brightest root independently in every job and assume it is
the same electronic state. Root ordering can change. Save all requested
excitation energies, oscillator strengths, electric transition dipoles, and
available excited-state descriptors. Track the same local bright pi-pi-star
state across frames and sites using, in descending order of preference:

1. NTO or transition-density overlap in a common atom/grid representation;
2. transition-dipole direction and magnitude continuity;
3. energy/oscillator-strength continuity plus inspection of ambiguous cases.

Flag every root switch or ambiguous assignment. Determine from the installed
TeraChem manual and a small validation calculation whether this version can
produce magnetic transition dipoles or rotational strengths. If it cannot,
state that clearly and reconstruct only the exciton-chirality contribution to
CD; do not label a normalized interaction-induced signal as absolute molar CD.

For every accepted `(frame, site)` result, retain a compact parsed record with
at least:

- frame index and time;
- site label and exact QM atom mapping;
- QM charge, multiplicity, method, basis, and embedding metadata;
- tracked root, excitation energy in eV and cm^-1, oscillator strength, and
  electric transition-dipole vector in a documented coordinate frame;
- state-tracking score and ambiguity flag;
- SCF/excited-state convergence status, wall time, GPU, input/output hashes.

## Coupling and spectral reconstruction

For the first pilot, pair the calculated site energies and transition dipoles
with the existing corrected frame-specific `J_i` values and geometries from
`coupling_nvt_production_cr2_1000_20260721`. Verify frame identity rather than
joining tables by row position alone. This is the quickest clean test of the
site-energy hypothesis.

As a second-stage sensitivity test, use the two frame-specific TeraChem
transition densities to calculate `J_i` consistently for each pair. Reuse the
audited Coulomb/TDC infrastructure and its corrected reciprocal-distance unit
conversion. Do not substitute an identical rigid transition-density template
at both sites when claiming a fully frame-specific TDDFT coupling. Compare this
TDDFT `J_i` distribution with the retained STEOM-density distribution.

Diagonalize each frame's nondegenerate Hamiltonian, retain its eigenvectors,
and transform the two site transition dipoles into exciton transition dipoles.
Implement a tested Rosenfeld/exciton-chirality expression for rotational
strength with explicit unit and sign conventions. Include intrinsic monomer
rotational strengths only if TeraChem actually supplies them. The generalized
implementation must pass these limits:

- `E1 = E2` recovers `E+/- = E0 +/- J` and the current degenerate model;
- `J = 0` gives localized site states and no interaction-induced exciton
  couplet;
- swapping site labels leaves absorption and CD spectra invariant;
- reversing the dimer handedness reverses the exciton-chirality CD sign;
- large `|E1-E2|` reduces excitonic mixing as `2|J|/Omega`.

Generate ensemble-averaged absorption and CD spectra with clearly separated
inhomogeneous sampling and assumed homogeneous/vibronic broadening. A common
calibration shift may align the mean TDDFT excitation energy to the retained
STEOM or experimental band origin, but it must not alter `Delta_i`. Compare
simultaneously with absorption and with the experimental `TD - TDX` or
`TDX - TD` trace using the documented subtraction/sign convention. Do not tune
independent arbitrary site shifts or linewidths merely to force agreement.

Report peak positions, apparent separation, integrated positive/negative CD,
absorption line shape, excitonic mixing, localization weights, and bootstrap
uncertainties. Analyse the joint distribution and correlations of
`E1`, `E2`, `Delta`, `J`, geometry, dipoles, `Omega`, and rotational strength;
averages of those quantities are not interchangeable.

## Required deliverables

Create a self-contained, restartable workflow in a new clearly named directory
or scripts; do not overwrite archived production data. Deliver:

1. extraction/input-generation code with deterministic frame/site mapping;
2. GPU batch launcher with bounded parallelism and resume support;
3. robust TeraChem parser and state-tracking diagnostics;
4. per-frame CSV/Parquet plus a JSON run manifest and hashes;
5. generalized nondegenerate absorption/CD code and unit tests;
6. convergence, correlation, spectrum, and ambiguous-root figures;
7. a concise `TERACHEM_SITE_ENERGY_REPORT.md` stating what was run, what remains,
   and whether the pilot supports or rejects site-energy detuning as an
   explanation while preserving the observed CD strength;
8. exact commands and an estimated GPU-hour cost for scaling to 100 and 1000
   frames.

Run existing relevant tests before and after the changes, add focused tests for
all new scientific invariants, and inspect all numerical outliers rather than
discarding them silently. Commit only source, compact results, and documentation
that are safe to publish. At the end, report the commit hash, tests, completed
frame/site counts, failures, GPU time, and the strongest conclusion justified by
the data. Do not claim that a matching peak separation proves the mechanism if
the calculated CD amplitude, sign, absorption, or state identity disagrees.

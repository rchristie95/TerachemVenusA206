# Tandem-Venus solvation and decoherence numerical test

> **Coupling correction (8 August 2026):** The shared-scale spectral sensitivity
> values in this report have been regenerated with the corrected production
> J = 32.82 +/- 1.55 cm-1 documented in `reference/orca_validation.json`.
> The independently calculated electrostatic-gap traces and 22--24 fs
> dephasing estimates are unaffected by the TDC unit correction.

## Scope

This is an equilibrium QM/MM electrostatic-gap test, not a final quantum-bath
calculation. Two contiguous 4 ps NVT segments were propagated from the retained
1 ns production checkpoint with a 2 fs timestep and the complete 69,471-atom
protein/water/ion system saved every 4 fs (2,000 frames; 8 ps total). The exact
retained AMBER CR2 charges and bonded parameters were transplanted at both
chromophore sites.

The S1 excited-minus-ground charge probe was reconstructed from all seven
printed STEOM NTO pairs (total represented occupation 0.99045993). Each orbital
was normalized on its cube grid before forming

\[
\Delta\rho = \sum_k n_k (|p_k|^2-|h_k|^2).
\]

The neutral difference density was partitioned onto the 41 physical QM atoms.
For each trajectory frame and each chromophore site, its electrostatic
interaction with the explicit MM charges gave an environmental energy-gap
shift. The probed site's own CR2 and stacked Tyr were excluded. Site and
differential energy-gap correlations were converted to a classical
second-cumulant coherence envelope.

## Main 8 ps result

| Environment | sigma(A-B) (cm-1) | rho(A,B), zero lag | A 1/e (fs) | B 1/e (fs) | differential 1/e (fs) | positive differential integral (fs) | classical T2*, 1/e (fs) |
|---|---:|---:|---:|---:|---:|---:|---:|
| All MM charges | 385.0 | 0.065 | 78.5 | 71.6 | 61.3 | 67.8 | **21.80** |
| Protein + water (ions removed from estimator) | 387.2 | 0.069 | 78.8 | 78.9 | 62.1 | 72.4 | **21.62** |
| Protein only | 368.8 | 0.124 | 97.0 | 67.6 | 60.8 | 73.1 | **22.94** |
| Water only | 115.8 | 0.212 | 197.9 | 263.6 | 169.0 | 148.3 | **73.91** |

The local protein electrostatics dominate the differential variance and the
short coherence lifetime in this fixed-charge calculation. Water is slower and
more common-mode, so water alone dephases the inter-site coherence much less
strongly.

The effective differential reorganization scale from the all-MM variance,
sigma^2/(2 kBT), is 355 cm-1 at 300 K. This number is only a classical
variance-based diagnostic, not a quantum reorganization energy.

## Robustness tests

- Independent contiguous segments: T2* = 22.40 and 21.81 fs.
- Eight separate 1 ps blocks: median 24.52 fs, standard deviation 3.73 fs,
  range 19.67--29.99 fs.
- Sampling every 4, 8, 16, and 32 fs: T2* = 21.80, 22.24, 24.21, and 22.52 fs.
- Removing explicit ions from the estimator: 21.80 -> 21.62 fs.
- Dominant NTO pair only: 21.95 fs.
- Seven-pair atomic probe corrected to reproduce the cube dipole exactly:
  22.01 fs.
- Full OpenMM PME cross-check at 8 fs cadence, using a neutral CR2-only probe
  on the actual 29 CR2 atom positions: T2* = 24.01 fs. The PME and matched
  minimum-image differential traces correlate at 0.974; the matched
  minimum-image result is 23.00 fs.
- Halving the PME finite-difference scale changes the first 100-frame PME
  energy traces by 0.00055 cm-1 RMS (maximum 0.0016 cm-1).

The PME site-solvation 1/e times are 71.1 and 62.0 fs. Its differential
correlation first crosses 1/e at 21.5 fs but subsequently rebounds; its positive
integral time is 74.3 fs. This non-monotonic behaviour is why a single fitted
"solvation time" is not an adequate description of the bath.

## Physical interpretation

The calculation supports a shared microscopic origin but not equality of the
timescales. Solvation is the relaxation/memory of the environmental
energy-gap fluctuation. Decoherence is the accumulated random phase generated
by that fluctuation. Its cumulant contains a double time integral of the
differential correlation, so both fluctuation amplitude and memory matter.

Here the site-solvation memory is about 60--80 fs in the PME/local-environment
tests, while the 1/e coherence lifetime is about 22--24 fs. The apparently
similar 60 fs number is therefore a bath-correlation time, not T2*. The large
differential fluctuation, roughly 370--385 cm-1, makes coherence decay before
the bath itself has forgotten its initial configuration.

This also shows that the 8.3 ps bulk-water Debye time is not the relevant sole
clock for the chromophore. The local protein/water response is strongly
multicomponent and contains sub-100-fs motion as well as slower solvent
relaxation.

## Consequence for the corrected-coupling CD sensitivity

The archived reconstruction's 60 fs assumption gives a Lorentzian HWHM of
88.48 cm-1. Using the PME value 24.01 fs gives HWHM = 221.07 cm-1 and
FWHM = 442.13 cm-1, larger than the corrected mean Davydov splitting
2J = 65.63 cm-1. The nominal two exciton centres do not move (522.98 and
524.79 nm), because T2 changes the broadening rather than
J. In the rerun reconstruction:

- normalized CD sign order remains negative-short/positive-long;
- the extrema move from 522.30/525.47 nm to 520.32/527.49 nm;
- extrema separation changes from 3.17 to 7.18 nm;
- the raw model peak magnitude falls to 17.3% of the 60 fs result (independent
  normalization of each curve would hide this loss).

Thus the CD couplet topology survives, but its width, extrema and raw amplitude
are materially affected. A 22--24 fs Lorentzian lifetime would not support a
resolved absorption doublet under the simple FWHM < 2J criterion; that threshold
is T2* = 161.7 fs for J = 32.82 cm-1.

## Important limitations

1. This is a classical high-temperature cumulant from a fixed STEOM difference
   density. It is not a rigorous quantum T2 and lacks detailed balance.
2. It omits geometry-dependent QM excitation energies, intrachromophore
   vibrations, non-Condon effects, and population relaxation.
3. AMBER/TIP3P is non-polarizable, so electronic solvent/protein polarization
   and mutual QM/MM induction are absent.
4. Eight picoseconds is adequate for the observed ultrafast decay and block
   checks, but not for converging slow protein or solvent modes.
5. The pre-solvated production box is net -2e after the two exact -1 CR2 charge
   sets are transplanted. The neutral difference probe makes the PME background
   contribution small, and ion exclusion changes T2 by only 0.2 fs, but a
   publication calculation should rebuild/neutralize the box explicitly.
6. The biexponential tail fits are unstable over this short record. First
   crossings, positive integrals, segment replication, and direct coherence
   envelopes are more defensible diagnostics.

## Route to a publication-quality QM/MM result

Use the present calculation as a screening layer. Next, compute
geometry-dependent vertical gaps on a statistically selected high-cadence
subset with electrostatic-embedding QM/MM (a tractable range-separated TDDFT or
sTDA level calibrated against retained STEOM points). Construct separate A and
B gap traces, their cross-correlation, and a quantum-corrected spectral density.
Then propagate a non-Markovian lineshape or reduced-density-matrix model. In
parallel, launch nonequilibrium ground- and excited-state QM/MM trajectories
from independent equilibrium snapshots to obtain the Stokes-shift solvation
response directly. That separates solvation relaxation from equilibrium pure
dephasing rather than identifying one fitted time with the other.

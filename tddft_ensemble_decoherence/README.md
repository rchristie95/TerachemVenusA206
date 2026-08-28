# Geometry-dependent TDDFT/TDC ensemble decoherence

This directory is the restartable work package for the tandem-Venus
geometry-dependent excitation and decoherence calculation. It was initialized
from repository commit `34775aa82847f654d77f8720182ddc36782095c7` on 22 July
2026. No molecular dynamics has been rerun, and no trajectory, checkpoint, or
failed electronic-structure result has been replaced or interpolated.

## Current input status

The validated compact records identify the requested inputs as:

- short segment 1: `full_highres_0to4ps_4fs.dcd`, 1,000 frames, 0--3,996 fs;
- short segment 2: `full_highres_4to8ps_4fs.dcd`, 1,000 frames, 4,000--7,996 fs;
- long production: the July 21 exact-AMBER-CR2 trajectory, 1,000 frames at
  `(frame + 1) ps`.

Those three raw coordinate trajectories are external to Git and are not
present on the local filesystem. The short compact gap traces and long compact
STEOM-density coupling/geometry tables are present, but neither contains the
instantaneous full-system coordinates required to construct the 44-atom
QM geometries and MM point-charge fields. Therefore they cannot be converted
into TDDFT inputs.

The older local `tc_tandem_nvt_1000/tandem_nvt_1000.pdb` is deliberately not
substituted: it has 69,418 atoms and 31 atoms per CR2, whereas the short audit
has 69,471 atoms and the July 21 production protocol uses the retained 29-atom
AMBER CR2 topology. See `audit/trajectory_audit.json`.

After the requested scope was reduced, that old file was used only for an
explicitly labelled feasibility preflight. Its CR2 geometries failed badly, so
no malformed chromophore was sent to TeraChem. One rigid-reference,
old-environment snapshot was run to exercise state tracking and signed TDC;
see `surrogate_pilot/README.md`. It is not production ensemble data.

## Validated electronic-structure convention

The retained same-geometry benchmark fixes the calculation at:

- QM region: all 29 anionic CR2 atoms plus the 12-atom Tyr203 phenol group and
  3 hydrogen link atoms (44 atoms total; 41 physical atoms);
- charge/multiplicity: `-1/1`;
- TeraChem: `wB97X-D3/6-311G**`, electrostatic point-charge embedding,
  `scf diis+a`, `threall 1e-13`;
- target: the bright anionic chromophore pi-pi* state, root 2 at the reference
  geometry, tracked rather than assumed to remain root 2;
- coupling screening: report vacuum TDC and divide by `epsilon_opt = 1.77` for
  the manuscript convention;
- link atoms are retained in the QM wavefunction but excluded from physical
  transition-density maps and TDC integrals.

`input_templates/` contains exact restart inputs derived from
`tc_tddft_44/tddft.in` and `tc_tddft_44/td.in`. The pilot uses five roots,
because the retained reference shows the physical bright root and its nearest
competitor within the first three roots; root count must be rechecked on the
beginning/middle/end trajectory pilot before production.

## Calibration

The retained TeraChem reference bright state is 3.43441587 eV and the final
ORCA STEOM electric-dipole state is 2.366633 eV. The initial constant correction
is therefore `Delta = -1.06778287 eV`. It changes absolute energies and
wavelengths but cancels from `E_A-E_B`, so it cannot change the differential-gap
correlation or `T2*`. It is not applied to transition densities or couplings.
See `calibration/calibration.json`.

## Commands

Regenerate the planned manifests without overwriting any existing classified
rows:

```bash
python tddft_ensemble_decoherence/scripts/initialize_manifests.py
```

Run the same-geometry reference pilot (skips an already completed output):

```bash
bash tddft_ensemble_decoherence/scripts/run_reference_pilot.sh
```

Once the exact external trajectories are restored, their paths and SHA-256
hashes must be recorded in `provenance.json` before the stratified trajectory
pilot or production jobs are enabled. A coordinate checksum is required for
every manifest row; missing coordinates are never inferred from compact energy
or coupling tables.

## Environment

- GPU: NVIDIA GeForce RTX 4080, 16,376 MiB
- driver: 580.126.09; reported CUDA capability: 13.0
- TeraChem: 1.97B-beta-251105, Git version
  `3fcf25291f21a40754a20c9c6450f71e006a2729`
- ORCA retained outputs: 6.1.1
- OpenMM: 8.4
- TeraChem Conda Python: 3.12.12

Raw trajectories, DCDs, cubes, scratch directories, and ORCA temporary files
remain external/untracked. Only compact numerical products, scripts, figures,
and reports are intended for a final commit after all checks pass.

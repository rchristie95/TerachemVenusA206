# Original ORCA-only tandem-paper pipeline

This manifest records the original scientific workflow behind
`manuscript/JPCB_tandem.tex`. It is an inventory and execution map, not an
installation guide or a Windows port. Paths are repository-relative. Several
recovered scripts retain the author's absolute Linux paths so that the original
provenance is visible; replace those paths in a new implementation.

## Provenance and versions

- Electronic structure: ORCA 6.1.1 RELEASE.
- Correlated method: ORCA's MDCI implementation of DLPNO-STEOM-CCSD. No external
  AutoCI program, input language, or binary was used. The `autoci_*` executables
  found in the local licensed ORCA distribution are not project files and are
  intentionally excluded.
- TDDFT: ORCA 6.1.1, wB97X-D3/6-311G**, RIJCOSX/def2-J, TightSCF, TDA, five roots.
- STEOM: def2-SVPD, AutoAux, RIJCOSX, TightSCF, TightPNO, five root-wise roots,
  `AddL2Term`, `TCutPNOSingles 1e-12`, `NDav 400`, `MaxIter 200`.
- Classical preparation: OpenMM/PDBFixer with AMBER14 XML force fields. The
  environment snapshot is in `environment.yml`; exact production settings are
  below.
- Global random seed: `20260618`. Optical dielectric: `1.77` in the released
  production data (`1.78` when rounded in manuscript prose).

## Directory and filename contract

```text
1MYW.pdb                              crystallographic source monomer
venus_dimer.pdb                       A206 interface geometry
tandem_build/                         retained linker-build intermediates
tandem_dimer.pdb                      FP1-linker-FP2 starting construct
tc_tandem_nvt_1000/                   production MD work directory
tandem_nvt_1000_clean.pdb             external, coupling-ready 1000-frame trajectory
anionic_build/                        deterministic monomer topology/coordinates
tc_simple_anionic/                    anionic monomer frame and frozen QM selection
tc_simple_old/                        original dimer-chain coordinate frame
tc_qmmm_opt_constrained/              constrained QM/MM geometry source
neo_model/orca_steom/                 STEOM templates, fields and geometries
neo_model/orca_dft/                   like-for-like ORCA TDDFT input
coupling_tandem_1000/                 Figure 4 numerical values
lineshape_tandem/                     Figure 5 spectrum values
multipole_out_correct/                multipole table/analysis values
reference/                            checksums and compact validation values
```

## Exact original order

1. Build the crystallographic interface with `python build_dimer.py`, producing
   `venus_dimer.pdb` from `1MYW.pdb`.
2. Build the published tandem with `python build_tandem_dimer.py --dimer
   venus_dimer.pdb --out tandem_dimer.pdb --workdir tandem_build`. The fixed
   linker is `SGLRSENLYFQGPREFCRYPAQWRPLESRPRTT` (33 residues, TEV motif
   `ENLYFQG`). The hand-retained intermediates in `tandem_build/` fix atom order
   and permit validation without reproducing PyMOL/OpenMM placement exactly.
3. Run the production unrestrained NVT calculation:

   ```text
   python run_nvt.py --pdb tandem_dimer.pdb --workdir tc_tandem_nvt_1000 \
     --nvt-steps 500000 --openmm-trajectory-interval 500 \
     --openmm-trajectory-file tandem_nvt_1000.pdb --seed 20260618 --no-video
   ```

   Defaults used by that command are pH 7.0, 10 A solvent padding, 0.15 M ionic
   strength, 300 K, Langevin friction 1/ps, 2 fs step, PME with 1.0 nm cutoff,
   and CUDA. No `--restrain-interface` flag was used. Thus 500,000 x 2 fs = 1 ns
   and reporting every 500 steps gives 1000 production snapshots.
4. Convert the trajectory to the coupling convention:

   ```text
   python tandem_unwrap.py --in tc_tandem_nvt_1000/tandem_nvt_1000.pdb \
     --out tandem_nvt_1000_whole.pdb --n-fp1 229 --n-link 33
   python tandem_to_coupling.py --in tc_tandem_nvt_1000/tandem_nvt_1000.pdb \
     --out tandem_nvt_1000_clean.pdb --rcut 27
   ```

   The first file is for whole-molecule visualisation. The second drops solvent
   and linker and assigns each barrel to chain A or B by its nearest CR2. Every
   frame 0 through 999, in trajectory order, was selected; there was no random
   subsampling. `reference/production_frames.txt` is the selection record.
5. Prepare the original monomer/QM inputs. `anionic_qm_setup.py` and
   `qmmm_tddft_pipeline.py` create the anionic monomer, frozen residue selection,
   AMBER topology, link atoms and constrained QM/MM input. The exact retained
   sources are `anionic_build/monomer_min.pdb`, `monomer_solv.prmtop`,
   `monomer_solv.inpcrd`, `tc_simple_anionic/qm_selection.json`, and
   `tc_qmmm_opt_constrained/qm_opt.xyz`.
6. Extract CR2 plus the Tyr203 phenol model with
   `neo_model/build_qm_phenol_relaxed.py`. The model contains 44 atoms in the
   exact order stored in `neo_model/orca_steom/geom_cthrp.xyz`, is hydrogen
   capped, has charge -1 and multiplicity 1. The selection is CR2 plus Tyr203
   residue id 202 atoms `CG CD1 CD2 CE1 CE2 CZ OH HD1 HD2 HE1 HE2 HH` and its
   link H. `neo_model/build_screen.py` records the tested region variants.
7. Construct the point-charge environment from the retained AMBER topology and
   coordinates. The build scripts use residue-complete atoms within 12 A of the
   QM region, exclude QM atoms and charges within 1.8 A of the boundary, then
   preserve the original MM net charge by distributing the removed charge. The
   exact 2350-charge production field is
   `neo_model/orca_steom/field.pc`; the other `field_*.pc` files are the original
   region-screening variants.
8. Run the like-for-like ORCA TDDFT input
   `neo_model/orca_dft/tddft_wb97xd3.inp`. `reproduce_paper.py --engine orca`
   copies the same `geom_cthrp.xyz` and `field.pc` into this directory and parses
   the brightest oscillator-strength root.
9. Run the production DLPNO-STEOM-CCSD input
   `neo_model/orca_steom/steom_phenol_svpd.inp`. The recovered `go_par.sh`,
   `go_serial.sh`, `go_ckpt.sh`, and `run_orca.sh` show the original execution
   conventions. The valid spectrum is printed before the optional
   `DoSTEOMNatTransOrb` post-step terminates in MDCI; validation must test for the
   parsed spectrum, not only the final process return code.
10. Generate NTO orbital cubes from the new ORCA GBW with ORCA's `orca_plot`
    facility for S1 orbitals 89--96. These generated cubes are deliberately not
    distributed. `rebuild_steom_multipair.py` combines the four retained NTO
    pairs `(92,93,0.94764284)`, `(91,94,0.02484332)`,
    `(90,95,0.00866944)`, and `(89,96,0.00518175)`, choosing pair signs by the
    ORCA transition dipole. `make_steom_specnorm_density.py` scales the result to
    the spectroscopic dipole. `build_steom_density.py` is the earlier dominant-
    pair reconstruction and `build_steom_cubes.py` creates plot cubes; both are
    retained because they were part of the original analysis trail.
11. Align the density with `align_steom_density.py`, fitting shared CR2 atoms
    from `tc_simple_anionic/monomer_relaxed.pdb` onto
    `tc_simple_old/classical_relaxed.pdb`. This produces
    `steom_transdens_specnorm_oldframe.npz` without changing the dipole norm.
12. Compute static TDC with `coupling_ensemble.py` on `venus_dimer.pdb` and
    trajectory TDC on `tandem_nvt_1000_clean.pdb`, using the aligned density,
    rigid placement, epsilon 1.77, and all 1000 frames. The exact per-frame
    published results are in `coupling_tandem_1000/coupling_samples.csv`.
13. Run `multipole_analysis.py` with the aligned density, old monomer frame and
    `venus_dimer.pdb`. The corrected values are in
    `multipole_out_correct/multipole_analysis.csv`.
14. Run `absorption_cd_spectra.py` using the tandem coupling distribution,
    per-frame geometry, `T2*=60 fs`, and the experimental splitting window.
    `lineshape_tandem/lineshape_data.csv` contains the numerical spectrum.
15. Run `open_quantum_dynamics.py --all` for Figure 1 and the sensitivity
    results. `make_paper_figures.py`, the `make_steom_*`/`make_qm_plots.py`
    helpers, `manuscript/make_tandem_panels.py`, and
    `manuscript/make_spectra_panels.py` assemble the remaining panels.

## Manuscript result map

| Manuscript result | Producing source/data |
|---|---|
| Figure 1, stochastic/ME dynamics and purity/coupling | `open_quantum_dynamics.py` |
| Figures 2-3, transition/difference densities | `build_steom_cubes.py`, `make_steom_*`, PyMOL helpers |
| Figure 4a-c, tandem separation/coupling/histogram | `coupling_tandem_1000/coupling_samples.csv`, `manuscript/make_tandem_panels.py` |
| Figure 5a-c, resolvability/absorption/CD | `absorption_cd_spectra.py`, `lineshape_tandem/lineshape_data.csv` |
| ORCA TDDFT reference | `neo_model/orca_dft/tddft_wb97xd3.inp`, `reference/orca_validation.json` |
| DLPNO-STEOM site energy/dipole | `neo_model/orca_steom/steom_phenol_svpd.inp`, `reference/orca_validation.json` |
| Static TDC and PDA | `coupling_paper_steom_static/` |
| 1 ns ensemble, J = 110.56 +/- 2.79 cm-1 | `coupling_tandem_1000/` |
| Multipole enhancement, TDC/PDA = 4.3624 | `multipole_out_correct/multipole_analysis.csv` |
| Algorithm 1 workflow | `reproduce_paper.py` plus steps 1-15 above |

The optional Q-Chem ladder discussed in the paper is method validation, not part
of the requested ORCA-only execution path. It is skipped by `--engine orca`.

## Environment variables and external files

Use `config/original.env.example` as the variable contract. The original runner
scripts also contain literal Linux paths; those document the actual run and are
not portable configuration.

The following cannot or should not be published:

- ORCA executables, libraries, licence material, and bundled `autoci_*`
  executables: obtain a licensed ORCA 6.1.1 distribution externally.
- `tc_tandem_nvt_1000/tandem_nvt_1000.pdb` (5.62 GB) and
  `tandem_nvt_1000_clean.pdb` (584 MB): generated production trajectories. Their
  expected paths and SHA-256 values are in `reference/SHA256SUMS`.
- ORCA GBW, NTO, cube, densities, scratch and MPI files: generated by steps 8-10.
  They are optional caches, not scientific source inputs.
- `steom_transdens_*.npz`: generated deterministically from the ORCA NTO cubes;
  scalar shape/norm/charge checks are in `reference/orca_validation.json`.

Everything needed to implement the scientific logic is present: fixed structure
and linker definitions, atom/residue selections, original coordinate frames,
point charges and exact ORCA inputs, NTO pair/state metadata, frame selection,
coupling/density algorithms, plotting/analysis code, and compact numerical
validation data. A new implementation does not need to infer a missing scientific
step from the manuscript.

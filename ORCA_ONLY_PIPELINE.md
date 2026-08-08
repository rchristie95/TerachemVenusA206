# Original ORCA-only tandem-paper pipeline

This manifest records the original scientific workflow behind
`manuscript/JPCB_tandem_round_2.tex`. It is an inventory and execution map, not an
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
  `FollowCIS`, root homing and per-root checks, `DTol 1e-5`, `AddL2Term`,
  `TCutPNOSingles 1e-12`, `NDav 400`, `MaxIter 400`.
- Classical preparation: OpenMM/PDBFixer with AMBER14 XML force fields. The
  environment snapshot is in `environment.yml`; exact production settings are
  below.
- Global random seed: `20260618`. Optical dielectric: `1.77` in both the
  released production data and the manuscript.

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
coupling_nvt_production_cr2_1000_20260721/  definitive main-text tandem-panel values
coupling_tandem_1000/                 legacy clean-PDB coupling analysis
lineshape_tandem/                     archived spectrum-panel values
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
5. Prepare and relax the original monomer/QM inputs. `anionic_qm_setup.py` and
   `qmmm_tddft_pipeline.py` create the anionic monomer, frozen residue selection,
   AMBER topology, link atoms and constrained QM/MM input. The retained
   TeraChem minimisation (`tc_qmmm_opt_constrained/opt.in`) uses
   wB97X-D3/6-31G*, the MM point-charge field, and freezes every atom in the
   274-atom QM region except the 29 physical CR2 atoms. Its gradient-converged
   coordinates are `tc_qmmm_opt_constrained/qm_opt.xyz`. The other exact
   retained sources are `anionic_build/monomer_min.pdb`,
   `monomer_solv.prmtop`, `monomer_solv.inpcrd`, and
   `tc_simple_anionic/qm_selection.json`.
6. Build the frozen composite electronic-structure geometry with
   `neo_model/build_qm_phenol_relaxed.py`: it substitutes the constrained
   QM/MM-relaxed CR2 coordinates into the classically minimised monomer, keeps
   the Tyr203 phenol coordinates from that monomer, and rebuilds the three link
   hydrogens. The resulting model contains 44 atoms in the exact order stored
   in `neo_model/orca_steom/geom_cthrp.xyz`, is hydrogen capped, has charge -1
   and multiplicity 1. The selection is CR2 plus Tyr203 residue id 202 atoms
   `CG CD1 CD2 CE1 CE2 CZ OH HD1 HD2 HE1 HE2 HH` and its link H.
   `neo_model/build_screen.py` records the tested region variants. The complete
   44-atom composite was not separately geometry-optimised.
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
   `neo_model/orca_steom/steom_phenol_svpd_robust2.inp`. This deterministic
   active-space variant enables `DoDbFilter` and sets `NRootsCISNAT` above the
   requested root count. The recovered `go_par.sh`,
   `go_serial.sh`, `go_ckpt.sh`, and `run_orca.sh` show the original execution
   conventions. The valid spectrum is printed before the optional
   `DoSTEOMNatTransOrb` post-step terminates in MDCI; validation must test for the
   parsed spectrum, not only the final process return code.
10. Generate NTO orbital cubes from the production ORCA GBW with ORCA's
    `orca_plot` facility for S1 orbitals 86--99. These generated cubes are
    deliberately not distributed. `build_capmasked_steom_density.py` combines
    all seven printed pairs (99.046% total NTO occupation), chooses pair signs
    from the ORCA right transition dipole, removes voxels assigned by a
    nearest-QM-atom Voronoi mask to link atoms 42--44, and least-squares scales
    the masked density to the ORCA right transition moment before the relative
    cutoff. The 259,277-point retained grid preserves the source-frame moment
    norm to within 0.02%; PDA and multipole moments are evaluated about the grid
    centroid so the small cutoff residual charge cannot create translation
    dependence. `rebuild_steom_multipair.py`,
    `make_steom_specnorm_density.py`, and `build_steom_density.py` retain the
    earlier four-pair/dominant-pair reconstruction trail but are not the source
    of the definitive coupling archive; `build_steom_cubes.py` creates plot
    cubes.
11. Align the cap-masked density with `align_steom_density.py`, fitting shared CR2 atoms
    from `tc_simple_anionic/monomer_relaxed.pdb` onto
    `tc_simple_old/classical_relaxed.pdb`. This produces
    `steom_transdens_capmasked_oldframe.npz` without changing the stored dipole
    norm.
12. Compute static TDC with `coupling_ensemble.py` on `venus_dimer.pdb`. For the
    definitive trajectory calculation, extract the two CR2 rigid transforms
    from the PBC-whole production DCD with `extract_cr2_transforms.py`, then run
    `coupling_dcd_steom.py` (normally through
    `run_coupling_production_cr2_1000.ps1`) with the complete cap-masked density,
    epsilon 1.77, no spatial binning, and all 1000 frames. The exact per-frame
    manuscript results are in
    `coupling_nvt_production_cr2_1000_20260721/coupling_samples.csv`;
    `coupling_tandem_1000/` is retained only as a legacy clean-PDB analysis.
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
| Figures 3-4, transition/difference densities | `build_steom_cubes.py`, `make_steom_*`, PyMOL helpers |
| Figure 5, tandem separation/coupling/histogram panels | `coupling_nvt_production_cr2_1000_20260721/coupling_samples.csv`, `manuscript/make_tandem_panels.py` |
| Archived resolvability/absorption/CD panels (not in the current main text) | `absorption_cd_spectra.py`, `lineshape_tandem/lineshape_data.csv` |
| ORCA TDDFT reference | `neo_model/orca_dft/tddft_wb97xd3.inp`, `reference/orca_validation.json` |
| DLPNO-STEOM site energy/dipole | `neo_model/orca_steom/steom_phenol_svpd_robust2.inp`, `reference/orca_validation.json` |
| Static TDC and PDA | `coupling_paper_steom_static/` |
| Definitive 1 ns ensemble, J = 32.82 +/- 1.55 cm-1 | `coupling_nvt_production_cr2_1000_20260721/` |
| Static multipole correction, TDC/PDA = 1.2273 | `multipole_out_correct/multipole_analysis.csv` |
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

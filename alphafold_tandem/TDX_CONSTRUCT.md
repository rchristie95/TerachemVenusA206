# dVenus-TDX: sequence, structure, and what it takes to simulate it

2026-09-01. Companion to `ALPHAFOLD_REGISTER.md`, which is about the TD
register and applies here unchanged.

TDX is the single-chromophore control of Nguyen et al., *Biophys. J.* **124**,
4293-4309 (2025): the same tandem dimer as TD, with one barrel's chromophore
extinguished. Their published signal is the difference
`dVenus-TDX - dVenus-TD`, so TDX is half of every number this project compares
against, and until now it existed here only as an electrostatic model
(`terachem_site_energy_cd/scripts/prepare_tdx_partner_variant.py` neutralises
the partner CR2's 29 MM charges). This note adds the construct itself.

## 1. The sequence

**One residue.** TDX is TD with the chromophore tyrosine of one barrel mutated
to cysteine -- Amber, the Vogel lab's non-fluorescent Venus (Koushik, Chen,
Thaler, Puhl & Vogel, *Biophys. J.* **91**, L99, 2006; Addgene #27798). The
tyrosine is the phenol that becomes the chromophore's anionic arm, so without
it the Gly-Tyr-Gly triad never cyclises: the barrel folds, docks, and holds no
chromophore and no -1 e.

| | |
|---|---|
| parent | Addgene #250938, pdVenus-TD, 511 aa (`pdVenus-TD_Addgene250938_511aa.fasta`) |
| chromophore triads | GYG at 66-68 and 338-340; tyrosines at **67** and **339** |
| mutation, dark barrel 2 | **Y339C** (barrel-local Y67C) |
| mutation, dark barrel 1 | **Y67C** |
| length | 511 aa, unchanged |

```bash
python3 alphafold_tandem/make_tdx_sequence.py                  # dark barrel 2
python3 alphafold_tandem/make_tdx_sequence.py --dead-barrel 1  # dark barrel 1
```

writes `pdVenus-TDX_derived_511aa_dark{1,2}.fasta` and the matching
`af3_tdx_job_request_dark{1,2}.json`. Nothing is transcribed by hand: the
parent FASTA is read, the GYG triads are located, one character changes.

Barrel 2 is the default dark barrel so the surviving chromophore keeps the
residue numbering it has in every TD artefact (CR2 at PDB residue 66, sequence
Tyr67).

### What is derived and what is verified

`pdVenus-TD_Addgene250938_511aa.fasta` is the deposited TD sequence and is
verified. The TDX sequence here is **derived from it**, not read off Addgene:
this container's network policy blocks `addgene.org` (and everything else
outside GitHub and the package registries), so neither a browser nor a fetch
can reach the deposit page from inside the session. The two facts that need
five seconds of eyeballing against **https://www.addgene.org/250939/** when you
are at a machine that can open it:

1. the insert is 511 aa and identical to #250938 except at the chromophore;
2. the mutated residue is the chromophore tyrosine and it becomes **C**
   (if the deposit instead deletes the tyrosine, or uses a different
   substitution, rerun the generator -- it is one character);

and one that decides nothing physically but should be recorded: **which**
barrel is dark, N-terminal or C-terminal. Both variants are generated.

## 2. The structure

There is no crystal structure of the tandem construct, TD or TDX, so the fold
has to be predicted -- but *predicted* is not the same as *taken wholesale*.
`ALPHAFOLD_REGISTER.md` is the standing result: AF3 folds each barrel to
0.70 A and then reports 18.7-21.8 A PAE **between** them, and its inter-barrel
register sits 11.7 A from crystal. AF settles the barrels and says nothing
about how they sit relative to each other, and the coupling depends entirely on
how they sit relative to each other.

So TDX is built the way the v4 TD structure was built: crystal register from
`tandem_dimer_v4_af3tails.pdb`, local conformation only from AlphaFold.

```bash
python3 tandem_tdx_build/build_tdx_structure.py --from-af3 fold_..._model_0.cif
```

writes `tandem_tdx_build/tandem_dimer_tdx_dark2.pdb`: 509 residues, because the
one CR2 becomes the three standard residues it formed from. The triad is
superposed on the 15 flanking barrel residues either side, backbone atoms only,
and the graft aborts above 1.5 A -- the same procedure
`tandem_v4_correct_seq/graft_af3_tails.py` used for the C-terminal tails, and
the reason an AlphaFold model is needed at all. The surviving chromophore moves
0.000 A and the live-CA2-to-dead-Cys-CA distance is 26.334 A, the TD CR2-CR2
value.

### Why the AlphaFold model is not optional

The obvious shortcut is to un-mature the existing CR2 in place: every immature
heavy atom has a counterpart in the mature chromophore, because the
imidazolinone ring *is* the triad backbone.

    Gly  N  CA  C        <- N1  CA1  C1
    Cys  N  CA  CB  C O  <- N2  CA2  CB2 C2 O2
    Gly  N  CA  C  O     <- N3  CA3  C3  O3

Drop the phenol, put SG on the tyrosine's CB->CG vector to keep chi1, rebuild
the first glycine's carbonyl O -- the one lost as water when the ring closed --
and you are done. It does not work, and the builder refuses to write the
result:

```
CLASH 0.21 A: GLY336 O ... GLY338 N
CLASH 1.35 A: GLY336 C ... GLY338 N
refusing to write: the new triad makes 6 contact(s) below 2.2 A ...
```

In the mature chromophore C1 already has three substituents -- CA1, N2 and the
ring N3 -- so the sp2 lobe the carbonyl oxygen has to occupy is the one N3 is
sitting in. The rebuilt O lands 0.21 A from N3. Opening the ring is a
conformational change of the triad inside a fixed barrel, not an atom
substitution; there is no placement of the oxygen that fixes it, and a 0.21 A
overlap is not something to hand a minimiser and hope. `--force-strained`
exists for inspecting the structure, not for simulating it.

AlphaFold folds the unreacted triad natively -- open ring, sp3 CA, real phi/psi
-- which is exactly and only what is taken from it.

### Running the AlphaFold Server job

The server needs a Google sign-in and cannot be driven from this container, so
this part is yours to click:

1. open https://alphafoldserver.com
2. upload `alphafold_tandem/af3_tdx_job_request_dark2.json` (its "upload JSON"
   button takes this file as-is -- `dialect: alphafoldserver`, `version: 3`)
3. the job reuses the TD run's model seed `425000176` and
   `useStructureTemplate: true`, so the TDX and TD predictions differ by the
   construct and nothing else
4. unzip the results into `alphafold_tandem/af3_tdx_server_run/`
5. `python3 alphafold_tandem/analyse_af.py` for the register diagnostics, and
   `build_tdx_structure.py --from-af3 <model_0.cif>` for the triad graft

Expect the same verdict as TD: excellent barrels, meaningless register, low
pLDDT across the 240-272 linker. That is the expected result, not a failed run,
and it is the reason the register comes from the crystal.

## 3. Running the simulation

`run_nvt.py` required **exactly two** CR2 residues in two places, which a
single-chromophore construct can never satisfy. It now requires at least one.
The only remaining two-chromophore assumption is `--restrain-interface`, which
tethers the two CR2 centroids and has nothing to tether in TDX; it now says so.

```bash
bash tandem_tdx_build/run_md_tdx.sh
```

| | |
|---|---|
| timestep | 1.0 fs |
| steps between frames | 2000 (2 ps) |
| frames | 500 |
| NVT steps | 1,000,000 = **1.000 ns** |
| NPT equilibration | 200,000 steps = 200 ps, matching v4's 200 ps |
| seed, force field | 20260618, amber14-all / tip3pfb -- v4's, deliberately |

Those four numbers only close on a 1 fs timestep: 500 x 2000 x 1 fs = 1 ns. The
TD trajectories ran at 2 fs, where the same frame stride and frame count give
2 ns instead; `DT=2.0 STEPS_PER_FRAME=1000 bash tandem_tdx_build/run_md_tdx.sh`
gives the same 500 frames over the same 1 ns at the TD timestep. Either way the
sampling matches v4's nanosecond thinned by two.

## 4. What is still TD-only

`terachem_site_energy_cd/scripts/prepare_production_frame.py` and
`prepare_candidate2_frame.py` require two CR2 residues for the partner-CR2
embedding, so QM/MM site energies on a structural TDX frame need the same
generalisation the MD driver just got. The existing TDX site-energy control
(`prepare_tdx_partner_variant.py`) sidesteps this by neutralising the partner
CR2's charges in a completed TD frame, which isolates the electrostatic part of
the shift; a structural TDX trajectory is the other half of that comparison.

## 5. Handing the browser steps to a local session

`TDX_LOCAL_AGENT_BRIEF.md` is the standalone task brief for a Claude Code
session on a machine with a browser: verify the sequence against Addgene
#250939, submit the AlphaFold job, graft, build and run. It repeats the two
constraints that matter -- the register is never imported, and the in-place
un-maturation is refused rather than forced.

## 6. Tests

`tests/test_tdx_build.py` guards both halves: that the sequence differs from
the deposit by exactly the one chromophore tyrosine, that the AlphaFold job
request is submittable and carries the TD seed, that the in-place build is
refused, and that the graft reproduces a known triad exactly (it is run against
a synthetic model built from the v4 structure itself, so the right answer is
the identity). It also pins the 1 fs x 2000 x 500 = 1 ns arithmetic.

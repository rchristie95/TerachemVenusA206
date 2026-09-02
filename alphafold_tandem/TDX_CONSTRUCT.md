# dVenus-TDX: sequence, structure, and what it takes to simulate it

2026-09-01. Companion to `ALPHAFOLD_REGISTER.md`, which is about the TD
register and applies here unchanged.

> **2026-09-02.** The Addgene check below was done and it changed the
> construct: #250939 **deletes** the barrel-2 chromophore tyrosine (510-aa
> insert), it does not mutate it to cysteine. Everything downstream was
> regenerated for the deposited sequence -- FASTA, AlphaFold job, the AF3 run in
> `af3_tdx_server_run/`, the builder and the tests. `python3
> alphafold_tandem/verify_tdx_against_addgene.py` rechecks the sequence offline.

## 1. The sequence

**One residue.** TDX is TD with the chromophore tyrosine of one barrel
**deleted**. The tyrosine is the phenol that becomes the chromophore's anionic
arm, and it is also the middle of the Gly-Tyr-Gly triad that cyclises to form
the imidazolinone; without it there is no triad to cyclise and nothing to
form. The barrel folds, docks, and holds no chromophore and no -1 e.

| | |
|---|---|
| parent | Addgene #250938, pdVenus-TD, 511 aa (`pdVenus-TD_Addgene250938_511aa.fasta`) |
| deposit | Addgene #250939, pdVenus-TDX, **510 aa** (`pdVenus-TDX_Addgene250939_deposit_510aa.fasta`) |
| chromophore triads in the parent | GYG at 66-68 and 338-340; tyrosines at **67** and **339** |
| lesion, dark barrel 2 (the deposit) | **Y339del** (barrel-local Y67del), GYG -> GG |
| lesion, dark barrel 1 | **Y67del** |
| length | **510 aa**, one shorter than TD |

```bash
python3 alphafold_tandem/make_tdx_sequence.py                  # dark barrel 2
python3 alphafold_tandem/make_tdx_sequence.py --dead-barrel 1  # dark barrel 1
```

writes `pdVenus-TDX_derived_510aa_dark{1,2}.fasta` and the matching
`af3_tdx_job_request_dark{1,2}.json`. Nothing is transcribed by hand: the parent
FASTA is read, the GYG triads are located, one character is removed -- and for
the deposited barrel the generator refuses to write unless the result is
byte-identical to the deposit FASTA.

`--mode substitute` still produces the Y339C chain (511 aa) this note used to
describe. That is Amber, the Vogel-lab non-fluorescent Venus (Koushik, Chen,
Thaler, Puhl & Vogel, *Biophys. J.* **91**, L99, 2006; Addgene #27798) -- a real
protein, and not the one #250939 contains. It is kept reproducible, named
`..._amberY67C`, and is not the default.

Barrel 2 is the default dark barrel because it is the one the deposit darkens
("The second dVenus lacks a chromophore"), and because the surviving chromophore
then keeps the residue numbering it has in every TD artefact (CR2 at PDB residue
66, sequence Tyr67).

### What is derived and what is verified -- checked 2026-09-02

`pdVenus-TD_Addgene250938_511aa.fasta` is the deposited TD sequence and is
verified. The TDX sequence here is derived from it, and on 2026-09-02 that
derivation was checked against **https://www.addgene.org/250939/**, whose only
sequence is an Addgene-verified NGS result for the whole 4450 bp plasmid. The
check overturned what this note previously assumed -- the deposit deletes the
chromophore tyrosine, it does not mutate it to cysteine -- so the evidence is
recorded here in full rather than summarised:

* The plasmid's single 559-aa ORF is the pRSETb 6His leader
  (`MRGSHHHHHH...MELLPVAT`, 49 aa, in
  `pdVenus-TDX_Addgene250939_prsetb_leader_49aa.fasta`) followed by a **510-aa**
  insert. Aligned against the 511-aa #250938 insert, `difflib` reports exactly
  one opcode: `delete ref[338:339]='Y'` -- residue **Y339**, the barrel-2
  chromophore tyrosine, and nothing else anywhere in the 510 residues.
* At the codon level barrel 1 reads `CTC GGC TAC GGC CTG CAG` = L-G-**Y**-G-L-Q,
  an intact triad; barrel 2 reads `CTG GGC GGG CTG CAG` = L-G-G-L-Q. There is no
  tyrosine codon and no cysteine codon in its place -- the codon is gone, in
  frame.
* Addgene's metadata agrees: insert 1549 bp and total vector 4450 bp for #250939
  against 1552 bp and 4453 bp for #250938. Exactly 3 bp shorter, which a point
  substitution cannot produce.
* The deposit's purpose line, "The second dVenus lacks a chromophore", confirms
  the C-terminal barrel is the dark one.

`python3 alphafold_tandem/verify_tdx_against_addgene.py` reruns the comparison
offline against the stored deposit translation and exits non-zero if the derived
dark-2 sequence ever stops matching it. The generator makes the same check
before it writes.

What this cost: the earlier Y339C reading (Amber) was a reasonable inference
from the literature -- Amber is the Vogel lab's non-fluorescent Venus and the
obvious way to darken a barrel -- but it is not this plasmid, and the difference
propagates. A deletion shortens the chain, shifts every downstream residue
number by one, and turns a three-residue triad into a two-residue one, which is
a different local backbone inside the dark barrel rather than a side-chain swap.
Sections 2 and 3 below are written for the deposited construct.

## 2. The structure

There is no crystal structure of the tandem construct, TD or TDX, so the fold
has to be predicted -- but *predicted* is not the same as *taken wholesale*.
`ALPHAFOLD_REGISTER.md` is the standing result, and section 1b there now
includes the TDX run: AF3 folds each barrel to 0.69 A and then reports 21-23 A
PAE **between** them, and the five models of a single seed put the second barrel
11.3 to 28.2 A from where the crystal puts it. AF settles the barrels and says
nothing about how they sit relative to each other, and the coupling depends
entirely on how they sit relative to each other.

So TDX is built the way the v4 TD structure was built: crystal register from
`tandem_dimer_v4_af3tails.pdb`, local conformation only from AlphaFold.

```bash
python3 tandem_tdx_build/build_tdx_structure.py \
    --from-af3 alphafold_tandem/af3_tdx_server_run/fold_dvenus_tdx_dark2_model_0.cif
```

writes `tandem_tdx_build/tandem_dimer_tdx_dark2.pdb`: **508 residues**, because
the one CR2 becomes the *two* standard residues the deleted-tyrosine triad
leaves. The pair is superposed on the 15 flanking barrel residues either side,
backbone atoms only, and the graft aborts above 1.5 A -- the same procedure
`tandem_v4_correct_seq/graft_af3_tails.py` used for the C-terminal tails. What
the run reports:

```
AF3 graft: anchored on 120 backbone atoms of 30 residues, superposition RMSD 0.21 A
surviving CR2: 1 residue, max coordinate change 0.000 A
all other residues: max coordinate change 0.000 A
chromophore separation CR2 CA2 <-> dead Gly(1) CA: 25.952 A (TD CR2-CR2 CA2-CA2 reference 26.334 A)
Gly(1) C - Gly(2) N: 1.321 A (peptide bond)
N-terminal junction C-N: 1.478 A (stretched, minimisation will take it out)
C-terminal junction C-N: 1.575 A (stretched, minimisation will take it out)
```

**The register is checked over every atom, not sampled.** "All other residues:
max coordinate change 0.000 A" is the statement that matters -- nothing outside
the two grafted residues moves, so the inter-barrel geometry is the crystal's by
construction rather than by inspection of one distance.

**The junctions are stretched and that is inherent.** The crystal leaves a slot
that a three-residue triad occupied, LEU335 C to LEU337 N = 5.145 A, and the AF3
Gly-Gly spans 5.546 A between the same flanking atoms. The 0.4 A excess has
nowhere to go but the two junction bonds, which come out at 1.48 and 1.58 A
against an ideal 1.33 A. That is the price of keeping the crystal register while
deleting a residue, it does not depend on the anchor window (1.54-1.59 A for
anchors from 4 to 25 residues), and it is small enough for the minimisation
before MD to take out locally without touching a register 26 A away. The builder
prints both junctions every run and refuses above 1.8 A, which would be a break
rather than a stretch.

### The in-place route: refused for Y->C, sound for the deletion

The shortcut is to collapse the existing CR2 in place rather than take geometry
from AlphaFold. For the *substitution* reading of the construct this was
impossible, and the reason is worth keeping because it is what sent us to
AlphaFold in the first place: the mature imidazolinone ring is
C1-N2-CA2-C2-N3, so C1 already carries three substituents -- CA1, N2 and the
ring N3 -- and the sp2 lobe the restored carbonyl oxygen needs is the one N3
occupies. The rebuilt O landed 0.21 A from N3 and the builder refused to write
it.

Deleting the tyrosine removes N2, and with it the problem. C1 keeps CA1 and N3,
the free lobe takes the carbonyl oxygen, and the 1.355 A C1-N3 ring bond becomes
the Gly-Gly peptide bond directly. So

```bash
python3 tandem_tdx_build/build_tdx_structure.py       # no --from-af3
```

now succeeds: no clashes, junctions untouched at the crystal's own 1.325 and
1.332 A, and a 1.355 A peptide bond between the glycines.

It is still not the production structure. Collapsing in place leaves the two
glycines exactly where the imidazolinone ring held them, with ring-sized bond
angles, in a barrel strand that is now one residue shorter -- a real
conformational change that nothing local to the crystal knows about. AlphaFold
folds the shortened strand natively, which is exactly and only what is taken
from it. Use the in-place build as a cross-check on the graft, not as the
starting structure.

### Running the AlphaFold Server job

Done for the deposited construct on 2026-09-02; `af3_tdx_server_run/` holds the
result. To repeat it, or to run the dark-barrel-1 variant:

1. open https://alphafoldserver.com (needs a Google sign-in)
2. upload `alphafold_tandem/af3_tdx_job_request_dark2.json` -- its "upload JSON"
   button takes the file as-is (`dialect: alphafoldserver`, `version: 3`), which
   is the point of generating it; do not retype the sequence into the web form
3. the job reuses the TD run's model seed `425000176` and
   `useStructureTemplate: true`, so the TDX and TD predictions differ by the
   construct and nothing else
4. unzip the results into `alphafold_tandem/af3_tdx_server_run/`
5. `python3 alphafold_tandem/analyse_af.py` for the register diagnostics, and
   `build_tdx_structure.py --from-af3 <model_0.cif>` for the graft

The verdict was the same as TD's: excellent barrels, meaningless register, a
linker below 30 pLDDT. That is the expected result, not a failed run, and it is
the reason the register comes from the crystal. Do not re-run seeds hoping it
improves -- section 3 of `ALPHAFOLD_REGISTER.md` shows the seed cannot probe it,
and the TDX run shows five models of *one* seed spanning 17 A.

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

## 5. The browser steps, and what they changed

`TDX_LOCAL_AGENT_BRIEF.md` was the task brief for a session with a browser:
verify the sequence against Addgene #250939, submit the AlphaFold job, graft,
build and run. Steps 1 and 2 were done on 2026-09-02 and step 1 changed the
construct, which is the case that brief was written to catch -- it said to
report a disagreement rather than work around it, and the disagreement was
real. The two standing constraints are unchanged: the inter-barrel register is
never imported from AlphaFold, and a refused build is not forced through with
`--force-strained`.

## 6. Tests

`tests/test_tdx_build.py`, 15 tests, guards both halves:

* the sequence differs from the #250938 parent by exactly one deleted
  chromophore tyrosine, for either dark barrel;
* the dark-barrel-2 sequence is byte-identical to the Addgene #250939 deposit;
* `--mode substitute` still reproduces the Amber Y339C chain and it is *not* the
  deposit;
* the AlphaFold job request is submittable and carries the TD seed;
* the graft reproduces a known pair exactly, against a synthetic model built
  from the v4 structure itself, so the right answer is the identity;
* the graft refuses a model with any residue between the two glycines -- a TD
  model, or the superseded Y->C reading -- and refuses a chain of the wrong
  length;
* the surviving chromophore and every atom outside the grafted pair are
  unmoved, and the output is 508 residues;
* the in-place collapse now produces a real peptide bond and no clashes;
* the production build from `af3_tdx_server_run/` passes all of the above;
* and the 1 fs x 2000 x 500 = 1 ns arithmetic still closes.

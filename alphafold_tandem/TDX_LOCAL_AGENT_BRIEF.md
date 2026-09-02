# Task brief: finish dVenus-TDX from a machine with a browser

For a Claude Code session running **locally**, where Chrome can reach
addgene.org and is signed in to the AlphaFold Server. The remote session that
prepared this branch could do neither: its egress policy allows GitHub and the
package registries only, and the AF Server needs a Google sign-in it cannot
hold. Everything that does not need a browser is already done and pushed.

    branch   claude/venus-tandem-tdx-structure-qk4pck
    context  alphafold_tandem/TDX_CONSTRUCT.md  <- read this first

There are two browser steps and one build. Do them in order; step 3 depends on
step 2, and step 1 can change what step 2 submits.

---

## Step 1 -- verify the sequence against Addgene

Open **https://www.addgene.org/250939/** (the TDX deposit; #250938 is the TD
parent already verified in `pdVenus-TD_Addgene250938_511aa.fasta`).

The TDX sequence on this branch is *derived*, not read off the deposit: it is
the TD sequence with the chromophore tyrosine of one barrel mutated to
cysteine, which is Amber, the Vogel lab's non-fluorescent Venus (barrel-local
Y67C). Check three things against the full protein/insert sequence on the
deposit page:

1. **Length is 511 aa** and the sequence is identical to #250938 apart from the
   chromophore.
2. **The mutated residue is the chromophore tyrosine and it becomes C.** In the
   511-aa numbering the two chromophore triads are `GYG` at 66-68 and 338-340,
   so the mutation is Y67C or Y339C. If the deposit instead *deletes* the
   tyrosine, or substitutes something other than cysteine, the derivation is
   wrong -- say so rather than working around it.
3. **Which barrel is dark**, N-terminal or C-terminal. Both variants are
   already generated; this decides which one is the real construct.

Do the comparison mechanically, not by eye:

```bash
python3 - <<'PY'
deposit = """<paste the deposit's protein sequence here>"""
deposit = "".join(deposit.split())
ours = "".join(l.strip() for l in
       open("alphafold_tandem/pdVenus-TDX_derived_511aa_dark2.fasta") if not l.startswith(">"))
print(len(deposit), len(ours))
print([(i+1, a, b) for i, (a, b) in enumerate(zip(deposit, ours)) if a != b])
PY
```

Record the outcome in `alphafold_tandem/TDX_CONSTRUCT.md` under "What is
derived and what is verified", replacing that paragraph with what the deposit
actually says and the date you checked it. If the mutation or the dark barrel
differs from the default, regenerate before step 2:

```bash
python3 alphafold_tandem/make_tdx_sequence.py --dead-barrel 1   # or 2
```

If Addgene shows only the plasmid map and not the insert protein sequence, the
"Sequences" tab (`https://www.addgene.org/250939/sequences/`) has the full
depositor sequence; translate the insert ORF rather than guessing.

## Step 2 -- run the AlphaFold Server job

`alphafold_tandem/AF_SERVER_TERMS_OF_USE.md` records the terms this project
accepted; nothing here changes them.

1. Open https://alphafoldserver.com
2. Upload **`alphafold_tandem/af3_tdx_job_request_dark2.json`** as-is (use the
   `dark1` file instead if step 1 says barrel 1 is the dark one). It is already
   in the server's own format -- `dialect: alphafoldserver`, `version: 3` -- and
   carries the TD run's model seed `425000176` with `useStructureTemplate:
   true`, so the TDX and TD predictions differ by the construct and nothing
   else. Do not retype the sequence into the web form.
3. Download the result and unzip it into
   **`alphafold_tandem/af3_tdx_server_run/`**, matching the layout of the
   existing `af3_server_run/` (the five `*_model_{0..4}.cif`, the
   `*_summary_confidences_*.json`, the job request, `templates/`).

**Expect the barrels to be excellent and the register to be meaningless** --
per-barrel pLDDT around 90, inter-barrel PAE near 20 A, pLDDT below 70 across
the linker. That is the established result for this construct
(`ALPHAFOLD_REGISTER.md`), not a failed run, and it is the reason step 3 takes
only the local triad from the model. Do not re-run seeds hoping the register
improves; section 3 of that note shows the seed cannot probe it (F = 0.30,
p = 0.93 over 7 seeds).

## Step 3 -- build the structure and check it

```bash
python3 tandem_tdx_build/build_tdx_structure.py \
    --from-af3 alphafold_tandem/af3_tdx_server_run/<...>_model_0.cif
```

Confirm in its output:

* `superposition RMSD` below ~0.6 A (it aborts above 1.5 A);
* `surviving CR2: 1 residue, max coordinate change 0.000 A`;
* `chromophore separation ... 26.334 A`, the TD CR2-CR2 value -- if this moved,
  the crystal register has been disturbed and something is wrong;
* `former ring bond Gly(1) C ... Gly(3) N` now **open** (> 2.5 A), which is the
  whole reason for using the AF triad;
* no `CLASH` lines. If any appear the script refuses to write, which is
  correct -- do not pass `--force-strained` to get past it.

Then:

```bash
python3 -m pytest tests/test_tdx_build.py -q     # 10 tests
python3 alphafold_tandem/analyse_af.py           # register diagnostics; edit ROOT first
```

`analyse_af.py` still points at a scratch directory from an earlier run; repoint
its `ROOT` at `af3_tdx_server_run/` before trusting its numbers, and add the TDX
row to the tables in `ALPHAFOLD_REGISTER.md` so the TD and TDX registers sit
side by side.

## Step 4 -- the simulation

```bash
bash tandem_tdx_build/run_md_tdx.sh
```

1 fs x 2000 steps per frame x 500 frames = 1.000 ns, after 200 ps NPT, at v4's
seed and force field. It needs the `run_nvt.py` on this branch: the old one
required exactly two CR2 residues and TDX has one. `--restrain-interface` is
not available for TDX and the script does not use it.

Note the timestep. Those three numbers (2000 steps, 500 frames, 1 ns) only
close at 1 fs; the TD trajectories ran at 2 fs, where the same stride and frame
count give 2 ns. `DT=2.0 STEPS_PER_FRAME=1000 bash tandem_tdx_build/run_md_tdx.sh`
gives the same 500 frames over the same 1 ns at the TD timestep, if matching TD
matters more than matching the stride.

---

## Do not

* **Do not import AlphaFold's inter-barrel register.** Only the three residues
  of the dead triad come from the model, superposed on the flanking barrel
  residues. AF3's register sits 11.7 A from crystal and the coupling depends
  entirely on the register.
* **Do not un-mature the chromophore in place** to skip the AF job. It is
  implemented, it is refused, and `TDX_CONSTRUCT.md` section 2 says why: the
  restored carbonyl oxygen has to occupy the lobe the ring nitrogen is in and
  lands 0.21 A from it.
* **Do not edit the derived FASTA by hand.** `make_tdx_sequence.py` regenerates
  it from the deposit parent; a hand edit loses the provenance that makes a
  TD-vs-TDX difference interpretable.

## Done when

Addgene is checked and the note says what it said; `af3_tdx_server_run/` is
populated; `tandem_dimer_tdx_dark{1,2}.pdb` builds from the model with no
clashes and the register intact; the tests pass; and the MD is running or its
first frames are written. Commit to the same branch.

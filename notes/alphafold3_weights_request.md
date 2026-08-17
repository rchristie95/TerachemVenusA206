# AlphaFold 3 model-parameters request — prepared draft

> **Status 2026-08-17: NOT PURSUED.** The weights were downloaded and then
> deleted unused. This machine cannot run the standard pipeline — the genetic
> databases need ~650 GB against 310 GB free, and neither the AlphaFold 3 source
> nor a JAX/CUDA stack is installed. The workaround (`--norun_data_pipeline`
> with precomputed MSAs) requires an external MSA source. Kept as a prepared
> draft in case the register question is revisited; see the last section for why
> it would not have settled δ.


**Form:** <https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform>

**Turnaround:** Google aim to respond in 2–3 business days. Two emails: an
immediate acknowledgement, then a download link on approval.

You submit this yourself — it accepts a licence in your name and asserts your
institutional eligibility.

---

## What you are agreeing to

From [WEIGHTS_TERMS_OF_USE.md](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md):

- **Eligible:** "non-commercial organizations (i.e., universities, non-profit
  organizations and research institutes, educational, journalism and government
  bodies)". A researcher at such an organisation qualifies provided they are not
  acting on behalf of a commercial organisation. University of Portsmouth /
  Surrey collaboration qualifies.
- **Must not:** carry out commercial activities *including research on behalf of
  commercial organisations*; share the weights outside your organisation; use
  outputs to train competing biomolecular structure-prediction models.
- **May:** publish, share and adapt AlphaFold 3 *outputs* (so predicted
  structures and anything derived from them can go in the paper), disclosing any
  modifications.
- Theoretical modelling only; not validated for clinical use.

Note the sharing clause: weights stay inside your organisation. If the Surrey
side needs them, they request separately.

---

## Draft text for the research-description field

> I am an academic researcher studying excitonic coupling between chromophores in
> tandem fluorescent-protein dimers, combining QM/MM transition-density
> calculations (DLPNO-STEOM-CCSD, TDDFT) with all-atom molecular dynamics.
>
> I intend to use AlphaFold 3 to predict the relative arrangement of the two
> β-barrel subunits in a tandem dimer of the yellow fluorescent protein Venus,
> in which two ~238-residue barrels are joined by a flexible 33-residue linker.
> The inter-barrel orientation is poorly constrained by the available
> crystallographic data — the crystal structure (PDB 1MYW) is of a van der Waals
> dimer rather than the linkered tandem construct — and it is the dominant
> uncertainty in our computed excitonic coupling. Generating an ensemble of
> predicted arrangements would let me test whether the register sampled by our
> molecular dynamics is representative, and compare against a recently published
> spectroscopic determination.
>
> This is non-commercial academic research. Outputs would be used for structural
> comparison and published as part of a peer-reviewed methods study; the model
> parameters would not be shared outside my institution.

Adjust the affiliation/PI details to match what the form asks for.

---

## Practical reality on this machine

**Your RTX 4080 (16 GB) is below the tested spec, but the target may just fit.**

- Officially supported: 1× A100 80 GB or 1× H100 80 GB.
- Documented lower-memory tiers: A100 40 GB up to 4,352 tokens; **V100 16 GB up
  to 1,280 tokens** with unified memory; P100 up to 1,024 tokens.
- Our target: 2 × 238 aa + 33-residue linker ≈ 509 residues, plus the chromophore
  as a modified residue (~1 token per atom) ≈ **570 tokens**. That sits inside the
  1,280-token envelope demonstrated on a 16 GB V100.

Unified-memory flags needed (spills to host RAM, significantly slower):

```bash
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
```

Two known frictions on this box: AF3 needs a working JAX/CUDA stack, and this
machine has a history of driver breakage after kernel bumps (the 4080 needs the
580.126.09 open module rebuilt and MOK-signed from `/usr/src`). Budget time for
the environment, not the inference.

Also note the CPU side: AF3's data pipeline needs genetic databases (~650 GB
download, needs fast local SSD). Currently 560 GB free — check before starting.

---

## What AF3 will and will not settle

**Will not:** the δ definition question. Their SI Note S4 says AlphaFold 3 "does
not consider any post-translational modifications", so it never builds the mature
chromophore — which is exactly why Cusick had to estimate δ from the OH→CB axis
of the tyrosine *precursor* after overlaying on the crystal. Running AF3
ourselves inherits that same proxy. See
[`cusick_comparison/results/delta_definition_offset.md`](../cusick_comparison/results/delta_definition_offset.md).

**Will:** the register question, which is the more valuable one. AF3 gives an
ensemble of candidate barrel–barrel arrangements — Cusick's Table S2 lists five,
spanning couplings of 29–43 cm⁻¹ — and our 1 ns (and 50 ns) MD cannot re-dock two
barrels, so it has only ever sampled the register it was built in. Comparing our
MD register against an independent AF3 ensemble is a real test of the leading
open uncertainty, and it is worth the request on those grounds even though it
will not resolve δ.

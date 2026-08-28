# The QM/MM embedding was the pathology (2026-08-28)

## Diagnosis

The residue-complete 12 A per-site cutoff is **charge-arbitrary and site-asymmetric**.
Between two structures of the same protein, site A's sphere lost Glu141/Glu403 and
gained Lys106/Arg429/Glu89 -- a **+3 e swing**. A unit charge at 8-12 A acting on the
chromophore's difference dipole is worth 800-3600 cm^-1, which covers the observed
+2550 cm^-1 shift in site A and explains why "which site is clean" flipped between
campaigns.

## Fix: full-system embedding

Every non-QM atom becomes a point charge; only the three severed-bond partners are
excluded (topological, not geometric -- a geometric cut deletes H-bond partners whose
membership flips with thermal bond length). Both flags already existed in
`prepare_candidate2_frame.py` and had never been back-ported here.

Frame 0, corrected 507-residue construct:

| scheme | charges | min/frame | site A err | site B err | GAP | gap err |
|---|---:|---:|---:|---:|---:|---:|
| full system | 74813 | 14.2 | 0 | 0 | **11** | 0 |
| union 15 A (symmetric) | 7056 | 5.6 | +142 | -194 | 347 | +336 |
| per-site 12 A (broken) | 2485 | 2.8 | +3139 | +459 | 2692 | +2681 |

(cm^-1). Full-system gives A/B MM net charge asymmetry of **exactly 0.0000 e**;
the per-site scheme gave 1.826 e.

## Why there is no cheap substitute -- three dead ends, all tested

1. **Neutralising the selection does nothing.** Forcing the truncated net charge to
   match full-system changed the field error by <1% at every radius (2381 -> 2367
   cm^-1 at 12 A). The monopole is NOT the dominant error.
2. **Enlarging the radius converges far too slowly.** Whole-residue selection leaves
   ~1134 cm^-1 of field error even at 30 A with 12967 charges; the error roughly
   halves per 18 A. A naive ATOM cutoff is worse -- it slices waters and charged
   residues, giving net charges like +11.3 e at 20 A and no improvement from 12 to 40 A.
3. **Symmetric union selection does not make the errors cancel.** The union is
   C2-symmetric but the structure is not (the two barrels differ by 1.2 A CA RMSD),
   so the site errors come out OPPOSITE in sign (+142, -194) and ADD in the gap.

**The decisive point:** with proper embedding the site-energy gap is ~10 cm^-1, which
is SMALLER than the truncation error of any affordable scheme. Truncation was
tolerable when the detuning was believed to be 576 cm^-1; it is not tolerable now.
That the artefact and the old "detuning" are the same order is itself the argument
that the old number was the artefact.

## Consequence to check

`<|Delta|>` = 576 cm^-1 underpins the paper's "detuning-dominated, localised"
conclusion. If the full-system ensemble confirms frame 0, that conclusion inverts.
ONE FRAME IS NOT A RESULT -- the 25-frame batch decides it.

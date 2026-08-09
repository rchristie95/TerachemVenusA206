# TeraChem site-energy/CD status

Starting commit: `e1c42f8cae66d1763b4fbe91e24fbe9288e99745`

Hardware: NVIDIA GeForce RTX 4080; driver 580.173.02

TeraChem: 1.97B-beta-251105; Python 3.12.12; OpenMM 8.4

## What completed

- Audited the 1,000-frame full-solvent rerun, topology, AMBER CR2 template,
  TeraChem binary, GPU, Python environment, and numerical settings into
  `results/run_manifest.json`.
- Prepared both central-frame sites as equivalent 44-atom CR2+Tyr203 models at
  charge `-1`, multiplicity `1`, with instantaneous link caps.
- Held the frame-0 MM atom inventory fixed, conserved charge at partially
  selected residues, and retained the complete partner CR2 with its audited
  RESP charge of `-0.999999999 e`.
- Ran seven-root wB97X-D3/6-31G* TDA energy calculations with no PCM. Both jobs
  converged: site A in 41.9 s and site B in 43.9 s. A second launcher pass
  reused both completed jobs without rerunning them.
- Parsed every root, oscillator strength, CI character, and electric transition
  dipole. The provisional bright `93 -> 95` candidates are A/root 1 at
  3.53379925 eV and B/root 2 at 3.29624849 eV.
- Added and passed five focused tests for the degenerate limit, zero-coupling
  localization/zero interaction CD, site-swap invariance, handedness sign
  reversal, and large-detuning mixing.

## Validation gate that did not pass

The available DCD hash is
`d854e18b90db301793f1550fcf28b8c1ccb4450732729ffe24deb1880c7fbab0`;
the coupling archive records
`5da1f8b2ce814dd04467b4d121c3ba70ccc8619852d045120d2e73086df53e5e`.
Frame identity with the retained `J_i` values is therefore unverified and the
workflow forbids a production join.

The provisional state candidates are also marked ambiguous because no
transition-density/NTO-overlap calculation was run for the corrected-charge
inputs. Their apparent detuning (1915.98 cm-1) is a smoke diagnostic only. The
archived frame-499 coupling was used solely to exercise the generalized
Hamiltonian/CD invariants; it is not paired scientific data.

The installed TeraChem manual documents electric transition dipoles but not
magnetic transition dipoles or rotational strengths. Accordingly the code can
reconstruct only relative interaction-induced exciton-chirality CD, not
absolute molar CD.

## Conclusion and remaining work

No conclusion about site-energy detuning, apparent CD peak separation, CD
strength, sign, or absorption is justified from this smoke test. The
deterministic 20-frame pilot, larger pilots, and all 1,000 frames remain
unlaunched. Restore the exact archived trajectory/topology, pass the hash gate,
generate overlap-based state diagnostics for both corrected-charge site jobs,
and validate a cold-start/restart pair before any expansion.

Focused tests: 5 passed. The broader selected suite gave 28 passed and 2
unrelated local-workspace failures: an untracked legacy coupling archive was
included by a repository glob, and the historical baseline lacks the NPZ path
expected by one migration test. Neither artifact was modified.

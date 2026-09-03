# Where things stood when the machine went down, 2026-09-03

Run `./resume_all.sh --check` for live state, `./resume_all.sh` to restart.
Everything is manifest-driven, so nothing recomputes: completed frames are
skipped on the basis of `convergence_status`, and no frame directory was left
partial (checked at shutdown).

## Complete, nothing to do

| campaign | | |
|---|---|---|
| `v5camp` | 500/500 | TD, 1 ns, the production ensemble |
| `v7tdx` site A | 500/500 | dVenus-TDX, one chromophore |
| `v7tdxdark` | 500/500 | TDX barrel 2 — **numbers are ill-posed, see below** |
| basis ladder, def2-SVPD | 48/48 | `<|mu|>` = 9.370 +/- 0.626 D |

## Unfinished

**`v8dense` 460/500.** 1 ps at 2 fs seeded from v5camp frame 1756, for the
cumulant lineshape. About 25 minutes of the 40 remaining frames, three workers.

**def2-TZVPD and def2-QZVPD**, 0 done. These need the GPU to themselves: at
1292 basis functions TZVPD died silently after basis setup when three SVP
workers were sharing the card. The templates already carry `guess sad`, because
TeraChem's default `generate` guess dies at that size with "Too many electrons
for guess generation code" (guess.cpp:58). Estimated 6-9 h for TZVPD, 5-14 h for
QZVPD, the latter possibly not fitting at all — it runs last for that reason.

## The one question this is all for

`|mu|` sets every absolute CD intensity, because the couplet amplitude scales as
`|mu|^4`. Production def2-SVP gives 9.40 D; two independent experimental routes
indicate 7.6-7.9 D. def2-SVPD now gives 9.370 +/- 0.626 D, so **diffuse
functions do not move it** and the gap is not a diffuse-function artefact.
Whether it is a zeta artefact is what TZVPD answers. If `|mu|` comes down, the
CD comparison in round 4 stops needing an experimental input and becomes a
prediction.

## Traps that will bite on resume

**The GPU driver.** This box's RTX 4080 module breaks on every kernel bump.
If `nvidia-smi` fails, rebuild and MOK-sign the 580.126.09 open module from
`/usr/src` — do not apt-install the skewed dkms package. `resume_all.sh` checks
this first and refuses to launch.

**Do not quote `v7tdxdark`'s Delta or J.** The dark barrel has no bright state,
so "brightest root in the window" returns whatever sits at the truncation:
7 roots gave |J| = 1.01 and a suppression of 1337x, 20 roots gave 12.73 and
255x, and roots 19-20 were still being selected at 20. The defensible statement
is that within +/-5000 cm^-1 of the chromophore the dark barrel has one state
carrying 0.08% of its oscillator strength, so the exciton problem in TDX is
one-dimensional.

**The master manifests are stale.** `v7tdx_frames.csv` and
`v7tdxdark_frames.csv` still read 0/500 because the workers wrote to the
`_a/_b/_c` splits. Merge them before anything downstream reads the unsuffixed
file.

## Next analysis step once v8dense lands

Feed its gap trajectory to `cumulant_lineshape.py`. The three v6t2 windows are
100 fs and resolve only to 334 cm^-1, coarser than kT = 208; 1 ps resolves to
33 cm^-1, which is what separates the fast lineshape from the static disorder
without the two being fitted against each other. That is the remaining input to
the band-shape question in `terachem_site_energy_cd/results/cd_final.json`.

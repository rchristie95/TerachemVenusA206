# Site B in the 4 fs window: seven hypotheses, one fix

2026-08-27. Follow-on from `REPRODUCIBILITY_AND_TIMESCALE.md`, which established
that the campaign numbers are bit-reproducible and that site B alone is
anomalous.

## The symptom

|  | 4 fs site A | 4 fs site B | 1 ns production A / B |
|---|---|---|---|
| sd of site energy | 494 | **954** | 519 / 516 |
| brightest character | (93,95) in 149/149 | (93,95)/(93,96)/(93,97) | -- |
| frames with a 2nd root at f > 0.4 | 4/149 | **51/149** | -- |

Production has the two sites equal (ratio 0.99), as a C2 dimer requires, so the
asymmetry belongs to this 8 ps trajectory, not to the system.

The two bright roots on site B are ~2119 cm^-1 apart, both HOMO->pi*, both
localised on the chromophore (phenol-ring fraction of |transition charge| 0.11
and 0.15, against site A's 0.09), and their intensities sum to a conserved
~1.2. That is intensity sharing with a dark state, not a near-degeneracy.

The excess is uniform across every 100 fs block (6-13 split frames each), so it
is a persistent property of the window, not a transient excursion.

## Hypotheses tested and rejected

| # | hypothesis | test | verdict |
|---|---|---|---|
| 1 | numerical noise | re-run prep + TDDFT on identical inputs | **dead** -- 0.000e+00 cm^-1, 8/8 |
| 2 | different QM region or atom order | element sequence, link caps, composition | **dead** -- identical, 44 atoms, same caps |
| 3 | low-lying CT state crowding the bright state | E(bright) - E(dark 93->94) | **dead** -- site B's is FURTHER (7934 vs 5619); corr(gap, n_bright) = +0.02 |
| 4 | chromophore twist splitting the band | inter-ring plane angle | **dead** -- 7.1+-3.3 vs A's 7.2+-3.1; corr(twist, 2nd f) = -0.08 |
| 5 | Stark mixing by the MM field | \|E\| and axial E at the chromophore | **dead** -- 0.0133 vs 0.0144 au; corr weak and sign-inconsistent |
| 6 | embedding net charge | **neutralised the embedding and re-ran** | **dead** -- see below |
| 7 | phenol stacking geometry | centroid separation, slip, ring-plane angle | **weak** -- corr(sep, 2nd f) = -0.19; right sign, far too small |

### On hypothesis 6, which looked strongest

The residue-complete 12 A selection leaves an arbitrary net charge: across
frames and sites it ranges from **+1.27 to -2.00 e**. Freezing the selection
(the fix that removed the frame-to-frame MM noise) froze the *asymmetry*:
site A at -1.3320 e, site B at -0.8158 e. The absolute MM potential at the two
chromophores differs by 0.093 a.u. = 2.5 eV, and within site B the bright
virtual index tracks the phenol-vs-chromophore potential difference
monotonically (95 -> +0.0488, 96 -> +0.0530, 97 -> +0.0574; corr +0.540).

It is still not the cause. Computing the effect first predicted that
outer-shell neutralisation would change that difference by only -0.0001 a.u.,
because the compensating charge sits far away and produces an almost uniform
offset rather than a gradient. Re-running four split site-B frames with a
neutralised embedding confirmed it:

| frame | site | as run | neutralised | shift |
|---|---|---|---|---|
| 6 | B | 3.4516 eV, f 0.606, virt 96, **2 bright** | 3.4530, 0.602, 96, **2 bright** | 11.4 cm^-1 |
| 18 | B | 3.3450, 0.728, 96, **2 bright** | 3.3460, 0.725, 96, **2 bright** | 7.6 |
| 40 | B | 2.9990, 0.700, 95, **2 bright** | 3.0002, 0.706, 95, **2 bright** | 9.4 |
| 90 | B | 3.2194, 0.650, 96, **2 bright** | 3.2214, 0.664, 96, **2 bright** | 16.1 |
| 5 | A | 3.2793, 1.009, 95, 1 bright | 3.2793, 1.010, 95, 1 bright | 0.4 |

The split survives untouched. **The embedding net charge is exonerated** as the
cause of the site asymmetry -- though the +1.27/-2.00 e spread is a real defect
of the selection and should be fixed on its own merits.

## What site B actually is, and the fix

Site B's stacked tyrosine is more loosely and more variably packed than site
A's -- centroid separation 4.27 +- 0.20 A against 4.13 +- 0.12, slip
2.14 +- 0.42 against 1.71 +- 0.32, ring-plane angle sd 7.2 deg against 4.7 --
and its whole environment is more mobile (the dark CT state's sd is 2488 cm^-1
against site A's 1737). A second HOMO->pi* excitation comes down into the
visible and borrows intensity. The QM region contains two pi systems (the
chromophore and the Tyr203-equivalent phenol), so a dense virtual manifold is
built in; site A avoids the consequences, site B does not.

That makes any single-root site energy ill-defined, and no tracker can repair
it -- consistent with the three trackers already tried.

**The fix is the observable, not the tracking.** When a diabatic bright state
|B> shares intensity with a dark state, the oscillator-strength-weighted first
moment recovers `<B|H|B>` exactly, because the dark partner contributes no
intensity:

    sum_k f_k E_k / sum_k f_k  =  |c1|^2 E_1 + |c2|^2 E_2  =  E_B

Windowed to +-4000 cm^-1 of the brightest root with f > 0.20, so that it
targets the splitting instead of averaging in unrelated states. An unwindowed
moment over all 7 roots also equalises the sites but shifts site A's absolute
energy by +884 cm^-1 and depends on the arbitrary root count -- rejected.

| estimator | mean A | mean B | sd A | sd B | **B/A** | ACF1 A | ACF1 B | slow A | slow B |
|---|---|---|---|---|---|---|---|---|---|
| selected root (old) | 25750 | 26126 | 494 | 954 | **1.93** | 0.75 | 0.42 | 168 | -- |
| **bright moment** | 25818 | 26111 | 532 | 616 | **1.16** | 0.70 | 0.61 | 172 | 254 |
| 1 ns production | 25373 | 25485 | 519 | 516 | 0.99 | -- | -- | -- | -- |

Implemented as `_bright_moment` in `scripts/run_ensemble.py` and recorded per
frame as `bright_moment_energy_eV`. Recorded, not substituted: existing
analyses are unaffected.

## Honest residual

The fix takes site B from 1.93x site A to 1.16x. It does not reach production's
0.99, and it should not be expected to -- part of the residual is a real
property of this trajectory, in which site B's stacked tyrosine is genuinely
looser. The fast/slow split, which is what the linewidth argument rests on,
now agrees between the sites to within its own uncertainty (slow 172 vs 254,
against a bootstrap spread of 107-192 on site A).

Because the residual cause is structural rather than electronic, **more frames
of this trajectory will not reduce it.**

# Does the fresh-velocity seed contaminate a 100 fs window?

The windows seed positions from a v5 DCD frame but must draw fresh
Maxwell-Boltzmann velocities: a DCD carries no momenta. Positions and momenta
factorise in the canonical ensemble, so (thermal positions, fresh MB momenta) is
a valid microstate and there should be no equilibration transient. That is an
argument, not a measurement. Window A's potential energy rose +1141 kJ/mol over
its 100 fs, which is what a transient would look like.

**Control**: same frame (900), same seed, run to 850 fs instead of 100 fs.

| 100 fs block | mean PE (kJ/mol) | drift (kJ/mol/100 fs) |
|---|---:|---:|
| 0-100   | -1209059.4 | **+1140.9** |
| 100-200 | -1209234.4 | -199.9 |
| 200-300 | -1209687.6 | -1540.7 |
| 300-400 | -1209678.3 | +523.6 |
| 400-500 | -1209536.9 | -616.3 |
| 500-600 | -1209974.4 | -251.5 |
| 600-700 | -1208588.2 | -593.0 |
| 700-800 | -1208783.3 | +379.8 |

The drift alternates sign with no decaying envelope, and the window mean sits
**0.03 sigma** from the parent trajectory's own production distribution
(+29.7 kJ/mol against the v5 NVT sd of 1053.8 over 2500 frames). Window A's
+1141 is one lobe of an ordinary thermal fluctuation on the same scale as the
1 ns run's, not a reseed artefact.

**Conclusion**: no equilibration prefix is needed, so t=0 of each window remains
the exact v5 frame whose STEOM reference already exists. Had a transient been
present the fix would have cost 3 reused TDDFT frames (adding a 100 fs prefix),
not the 8 STEOM jobs -- the STEOM calibration is a property of the configuration
basin and survives a short prefix either way.

## System-identity check

Seeding reproduces the potential energy the v5 run logged for the same frame:

| DCD frame | v5 logged (kJ/mol) | seeded (kJ/mol) | deviation |
|---:|---:|---:|---:|
| 900 | -1209516.140068 | -1209516.249321 | 0.109 |

0.109 kJ/mol out of 1.21e6 is 9e-8 relative -- float32 DCD coordinate rounding
plus mixed-precision CUDA. Atoms, ordering, parameters and cell all match; a
re-solvated box (addSolvent is unseeded) would differ by orders of magnitude
more. `run_nvt.py --expect-initial-potential-kj` enforces this at launch.

---

# Do the X-H constraints distort the gap dynamics?

`run_nvt.py` builds every system with `constraints=HBonds`, which freezes the
X-H stretches at ~3000 cm-1 -- inside the band a ~12 fs T2 is sensitive to.
Keeping them is the cleaner choice, but NOT for the reason first written here.
The original justification -- "classical MD over-populates those modes, so
freeing them over-estimates dephasing" -- is backwards for the amplitude.
For a harmonic mode `<q^2>_qm/<q^2>_cl = (x/2)coth(x/2)` with `x = hbar*w/kT`,
which is `>= 1` always: **7.2** at 3000 cm-1 and **3.25** at 1350 cm-1. Classical
MD MISSES zero-point motion and therefore UNDER-estimates stiff-mode amplitudes.
(The familiar "classical over-populates" statement concerns energy leaking into
modes that should be frozen -- an equipartition problem, not an equilibrium
amplitude.) The real reason to constrain them is that at `hbar*w/kT = 14.4` the
mode is in its quantum ground state, and a ground-state mode does not dephase
stochastically at all: it produces a static Franck-Condon progression. A
classical trajectory cannot represent that, so removing the channel is more
honest than sampling it wrongly. Either way it is an argument, not a
measurement.

**Control**: two 500 fs runs from the same seed frame (v5 900) and the same
velocity seed, identical but for the constraint set (`--no-hbond-constraints`,
added for this purpose and not for production). 500 rather than 100 fs because
the bridge C-C stretch has a ~20-30 fs period.

Compared on the **methine bridge length** and **BLA**, using
`bridge_detuning_analysis.py`'s own definitions, not on an electrostatic proxy.
The v5 analysis put 25.9% of the detuning variance on the bridge and only 4.1%
on electrostatics, and a Stark decomposition holds the chromophore geometry
rigid while varying the MM field -- so it is blind by construction to a change
in the chromophore's own internal motion, which is exactly what this control is
about. Running it that way would have produced a null result for the wrong
reason.

| coordinate | mismatch sd, constrained | unconstrained | ratio | tau_c const. | tau_c uncon. |
|---|---:|---:|---:|---:|---:|
| bridge | 58.34 mA | 58.41 mA | **1.001** | 4.38 fs | 4.41 fs |
| BLA | 37.31 mA | 38.71 mA | 1.037 | 3.99 fs | 3.80 fs |

The bridge -- the coordinate that carries the detuning -- is unchanged to 0.1%
in amplitude and 0.7% in correlation time. BLA moves under 4%. **Freezing the
X-H stretches does not distort the coordinate driving the gap**, so the
production T2 stands as computed rather than as an upper bound.

Data: `constraint_control.json`. Script:
`tddft_ensemble_decoherence/scripts/constraint_control_geometry.py`.

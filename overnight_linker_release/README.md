# Overnight linker-release experiment (2026-08-11 -> 08-12)

Tests the leading hypothesis of `notes/dimer_geometry_audit.md`: the tandem
was built in the 1MYW crystal register, which pins the 33-residue linker
(res 230-262) at 52.3 +/- 0.7 A end-to-end against ~38 A natural for a
disordered 33-mer, and 1 ns of MD cannot re-dock two beta-barrels — so no
alternative register has ever been sampled. Meanwhile both force fields move
the inter-chromophore axis angle *away* from the 131.3 deg that Nguyen's
limiting anisotropy requires (|cos alpha| = 0.660).

## Arms (both ~50 ns overnight, v3 ff19SB/OPC physics, 300 K NVT, 88.198 A box)

- `control/` — unbiased continuation from the v3 production end state.
  Null arm: does alpha stay pinned near ~95-105 deg with 50x more sampling?
- `release/` — steered harmonic on the linker end-to-end distance
  (CA229-CA263): ramp 52.7 -> 40 A over 2 ns, hold 2 ns, then k = 0 and
  unbiased for the rest of the night. If linker tension is what holds the
  crystal register, granting slack should let the interface move; the
  direction it moves (toward or away from 131 deg / |cos|=0.66, and what the
  chirality triple product does) is the result.

## Details that matter

- System rebuilt exactly as `run_nvt.py` v3: `amber19-all.xml` +
  `amber14/opc.xml` + retained generic CR2 template, then full CR2
  bonded/charge transplant from `anionic_build/monomer_solv.prmtop`
  (net charge verified 0).
- **Start-state trap**: `classical_relaxed.pdb` carries the *pre-NPT* CRYST1
  box (90.373 A) with unwrapped coordinates; re-imaging with that box creates
  ~20k fake overlaps and instant NaN. The true production box is 88.198 A;
  start positions/cell are the last production DCD frame
  (`v3_last_frame_xyz_ang.npy`, `v3_last_frame_box.json`).
- Metrics every 10 ps into `<arm>/<arm>_metrics.csv`: inter-CR2-long-axis
  angle alpha (OH minus imidazolinone-ring centroid — audit route 3, no
  dipoles needed), cos alpha, linker end-to-end, CR2 centroid separation,
  chirality triple product R_AB.(axA x axB), PE, T.
- Protein-only DCD every 25 ps, checkpoint every 100 ps.
- Throughput: ~280 ns/day solo, two concurrent arms share the RTX 4080.

## Morning

```bash
/home/robson/anaconda3/envs/TeraChem/bin/python /home/robson/PetaChem/overnight_linker_release/analyze_morning.py
```

Interpretation guide:
- release arm relaxes back to ~52 A and alpha returns to ~100 deg -> the
  register is held by the interface itself, not linker tension; hypothesis
  weakened.
- linker stays slack and alpha moves toward 131 deg (|cos| up toward 0.66)
  -> linker-tension hypothesis confirmed; the crystal register is an
  artifact of how the tandem was built, and the round-2 geometry needs the
  relaxed register.
- alpha moves but in the wrong direction / triple product flips -> new
  information either way; compare against the control arm before believing
  anything (the control shows what unbiased drift does on the same clock).

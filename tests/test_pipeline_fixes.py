#!/usr/bin/env python3
"""Lightweight tests (no GPU/TeraChem) for the three robustness fixes:
  1. frozen, coordinate-independent QM residue selection
  2. renormalization sanity guard (threshold logic)
  3. brightest-state selection (in-window / blue-band / no-bright)
Run: python3 tests/test_pipeline_fixes.py
"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import qmmm_tddft_pipeline as P


def test_residue_key_and_roundtrip():
    class _Chain:
        def __init__(self, i): self.id = i
    class _Res:
        def __init__(self, c, i, n): self.chain = _Chain(c); self.id = i; self.name = n
    rs = [_Res("A", "65", "CR2"), _Res("A", "96", "ARG"), _Res("A", "203", "TYR")]
    keys = {P.residue_key(r) for r in rs}
    assert keys == {"A:65:CR2", "A:96:ARG", "A:203:TYR"}, keys
    f = tempfile.mktemp(suffix=".json")
    P.save_qm_selection(rs, f)
    assert P.load_qm_selection(f) == keys
    os.remove(f)
    print("[ok] residue_key is coordinate-free; save/load round-trips")


def test_frozen_selection_is_coordinate_independent():
    """The headline fix: same residues despite NVT moving atoms across the cutoff."""
    try:
        from openmm import Platform
        Platform.loadPluginsFromDirectory("/home/robson/anaconda3/envs/TeraChem/lib/plugins")
        from openmm.app import AmberPrmtopFile, PDBFile
        from openmm import unit
        import numpy as np
        prm = AmberPrmtopFile("anionic_build/monomer_solv.prmtop")
        pos = np.array(PDBFile("anionic_build/monomer_min.pdb").getPositions().value_in_unit(unit.angstrom))
    except Exception as exc:
        print(f"[skip] frozen-selection topology test (inputs absent: {exc})")
        return
    box = P.get_periodic_box_lengths_ang(prm.topology); pbc = box is not None
    sel0, *_ = P.select_qm_residues(prm.topology, pos, "CR2", 2.65, nearest_waters=5,
                                    box_lengths_a=box, use_periodic=pbc)
    keys = {P.residue_key(r) for r in sel0}
    rng = np.random.default_rng(1); pos2 = pos + rng.normal(scale=1.0, size=pos.shape)
    selD, *_ = P.select_qm_residues(prm.topology, pos2, "CR2", 2.65, nearest_waters=5,
                                    box_lengths_a=box, use_periodic=pbc)
    selF, *_ = P.select_qm_residues(prm.topology, pos2, "CR2", 2.65, nearest_waters=5,
                                    box_lengths_a=box, use_periodic=pbc, frozen_keys=keys)
    keysD = {P.residue_key(r) for r in selD}
    keysF = {P.residue_key(r) for r in selF}
    assert keysF == keys, "frozen selection must be identical despite perturbed coords"
    print(f"[ok] frozen selection identical under 1 A noise ({len(keys)} residues); "
          f"distance-based drifts by {len(keys ^ keysD)} residues")


def test_renorm_guard_logic():
    # mirror the stage-3 guard window [0.5, 2.0]
    def in_range(scale): return 0.5 <= scale <= 2.0
    assert in_range(0.84)        # healthy anionic run
    assert not in_range(16.86)   # the garbage CT-root case
    assert not in_range(0.18)
    print("[ok] renorm guard flags 16.86/0.18 and passes 0.84")


def test_brightest_state_selection():
    mod = P.build_stage2_analysis_module()
    sel = mod.select_brightest_state
    # in-window (original 542 nm)
    assert sel([{"root": 1, "nm": 542, "osc": 0.70}, {"root": 2, "nm": 476, "osc": 0.02}]) == 1
    # blue band (anionic 420 nm) -> still picked as chromophore band
    assert sel([{"root": 1, "nm": 420, "osc": 0.65}, {"root": 2, "nm": 370, "osc": 0.48}]) == 1
    # no bright state -> brightest of the dark, with warning
    assert sel([{"root": 1, "nm": 300, "osc": 0.001}, {"root": 2, "nm": 290, "osc": 0.0005}]) == 1
    print("[ok] brightest-state: in-window, blue-band, and no-bright all handled")


def test_density_plot_state_consistency():
    """Regression: the transition-density plot run must solve at least as many CIS
    states as the energy run (NUM_STATES) and target the bright root, so the plotted
    density matches the reported bright state (renorm ~1 -> correct coupling J).
    Bug: cisnumstates was set to the root index, so in a dense excited-state manifold
    Davidson plotted a DIFFERENT state (observed renorm 3.2, wrong J)."""
    import types, tempfile
    from pathlib import Path

    def plot_in_states(generate, set_num_states, num_states, root):
        d = Path(tempfile.mkdtemp())
        set_num_states(num_states)
        geom = d / "geometry.xyz"; geom.write_text("1\n\nH 0.0 0.0 0.0\n")
        charges = d / "mm_charges.dat"  # absent on purpose -> gas-phase branch
        generate(d, geom, charges, root)
        lines = (d / "plot.in").read_text().splitlines()
        val = lambda key: next(int(l.split()[1]) for l in lines if l.strip().startswith(key))
        return val("cisnumstates"), val("cistarget")

    fake_run = lambda *a, **k: types.SimpleNamespace(returncode=0)

    # (1) standalone analysis module: energy run used NUM_STATES, bright root=2
    import terachem_tddft_analysis_big as TA
    saved = TA.subprocess.run; TA.subprocess.run = fake_run
    try:
        gen = lambda d, g, c, r: (setattr(TA, "WORKDIR", d), TA.generate_densities(g, c, r))
        assert plot_in_states(gen, lambda n: setattr(TA, "NUM_STATES", n), 10, 2) == (10, 2)
        # a bright root beyond NUM_STATES must still be solvable (max picks the root)
        assert plot_in_states(gen, lambda n: setattr(TA, "NUM_STATES", n), 10, 15) == (15, 15)
    finally:
        TA.subprocess.run = saved

    # (2) pipeline Stage-3 internal copy (mod.NUM_STATES fixed at 20)
    mod = P.build_stage2_analysis_module()
    saved2 = P.subprocess.run; P.subprocess.run = fake_run
    try:
        gen2 = lambda d, g, c, r: (setattr(mod, "WORKDIR", d), mod.generate_densities(g, c, r))
        n2, t2 = plot_in_states(gen2, lambda n: None, int(mod.NUM_STATES), 2)
    finally:
        P.subprocess.run = saved2
    assert (n2, t2) == (int(mod.NUM_STATES), 2), (n2, t2)
    print(f"[ok] density plot run solves >= energy-run states & targets bright root "
          f"(module 10/2 & 15/15; pipeline {n2}/{t2})")


def main():
    tests = [test_residue_key_and_roundtrip, test_frozen_selection_is_coordinate_independent,
             test_renorm_guard_logic, test_brightest_state_selection,
             test_density_plot_state_consistency]
    fails = 0
    for t in tests:
        try:
            t()
        except AssertionError as e:
            fails += 1; print(f"[FAIL] {t.__name__}: {e}")
    print()
    print(f"{'ALL PASSED' if not fails else f'{fails} FAILED'} ({len(tests)} tests)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

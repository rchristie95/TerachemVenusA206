#!/usr/bin/env python3
"""Apply the QM/MMpol polarization correction across a whole frame ensemble.

CPU-only post-processing of existing QM output, so it competes with nothing on
the GPU. The correction is PAIRED (same frame, with and without polarization),
which is far more powerful than comparing two trajectories: a 10-frame pilot put
the effect at -8.2% but only 0.54 SE from zero, and resolving it at 80% power
needs ~265 frames. This runs over every completed frame instead.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "terachem_site_energy_cd"))
sys.path.insert(0, str(ROOT / "terachem_site_energy_cd/scripts"))
from polarizable_embedding import polarization_shift_cm  # noqa: E402
from openmm.app import PDBFile  # noqa: E402
from mapping import reference_mapping, residue_by_key  # noqa: E402

CM_EV = 8065.54
PARSE = ROOT / "terachem_site_energy_cd/scripts/parse_terachem_tddft.py"


def bright_energy(qm_dir: Path, site: str, python: str, min_margin: float = 3.0):
    out = qm_dir / f"site_{site}/tddft.out"
    if not out.is_file() or "Job finished" not in out.read_text(errors="ignore"):
        return None
    js = qm_dir / f"site_{site}/energy_summary.json"
    if not js.is_file():
        subprocess.run([python, str(PARSE), str(out), "--json", str(js)],
                       capture_output=True)
    roots = json.loads(js.read_text())["roots"]
    ranked = sorted(roots, key=lambda r: -r["oscillator_strength"])
    if ranked[0]["oscillator_strength"] / max(ranked[1]["oscillator_strength"], 1e-12) < min_margin:
        return None
    return ranked[0]["energy_eV"] * CM_EV


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path,
                    default=ROOT / "terachem_site_energy_cd/results")
    ap.add_argument("--qm-prefix", default="v2_camb3lyp")
    ap.add_argument("--prep-prefix", default="v2_linkonly")
    ap.add_argument("--topology", type=Path,
                    default=ROOT / "tc_tandem_nvt_v2/solvated_protonated.pdb")
    ap.add_argument("--embedding-cache", type=Path,
                    default=ROOT / "terachem_site_energy_cd/results/v2_embedding_charges.npz")
    ap.add_argument("--probe", type=Path,
                    default=ROOT / "solvation_decoherence_test/steom_difference_probe_dipole_matched.npz")
    ap.add_argument("--polarizability-scale", type=float, default=1.0)
    ap.add_argument("--python", default="/home/robson/anaconda3/envs/TeraChem/bin/python")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    dq = np.load(args.probe)["atom_delta_q_e"]
    charges_all = np.load(args.embedding_cache)["charges_e"]
    pdb = PDBFile(str(args.topology))
    tyr_names = set(reference_mapping()[3])
    residues = list(pdb.topology.residues())

    cache: dict[tuple[str, str], tuple[np.ndarray, list]] = {}

    def context(prep_dir: Path, site: str):
        key = (prep_dir.name, site)
        if key in cache:
            return cache[key]
        meta = json.loads((prep_dir / "preparation.json").read_text())["sites"][site]
        cr2 = residue_by_key(pdb.topology, tuple(meta["cr2_key"]), "CR2")
        tyr = residue_by_key(pdb.topology, tuple(meta["stacked_tyr_key"]), "TYR")
        cr2_idx = [a.index for a in cr2.atoms()]
        tyr_idx = [a.index for a in tyr.atoms() if a.name in tyr_names]
        qm_idx = set(cr2_idx) | set(tyr_idx)
        link = {int(c["mm_partner_index_zero_based"]) for c in meta["link_caps"]}
        q_ground = charges_all[np.array(cr2_idx + tyr_idx, dtype=int)]
        syms = [a.element.symbol if a.element is not None else "C"
                for r in residues for a in r.atoms()
                if a.index not in qm_idx and a.index not in link]
        cache[key] = (q_ground, syms)
        return cache[key]

    rows = []
    for qm_dir in sorted(args.results_dir.glob(f"{args.qm_prefix}_frame_*")):
        frame = int("".join(c for c in qm_dir.name if c.isdigit())[-4:])
        prep = args.results_dir / f"{args.prep_prefix}_frame_{frame:04d}"
        if not (prep / "preparation.json").is_file():
            continue
        ea = bright_energy(qm_dir, "A", args.python)
        eb = bright_energy(qm_dir, "B", args.python)
        if ea is None or eb is None:
            continue
        shift = {}
        for site in ("A", "B"):
            q_ground, syms = context(prep, site)
            g = (qm_dir / f"site_{site}/geometry.xyz").read_text().splitlines()
            xyz = np.array([[float(v) for v in ln.split()[1:4]]
                            for ln in g[2:2 + 44]])[:41]
            raw = np.loadtxt(qm_dir / f"site_{site}/mm_charges.dat", skiprows=2)
            shift[site] = polarization_shift_cm(
                xyz, q_ground, dq, raw[:, 1:4], syms[:len(raw)],
                polarizability_scale=args.polarizability_scale)[0]
        rows.append((frame, ea - eb, shift["A"], shift["B"],
                     (ea - eb) + shift["A"] - shift["B"]))
        print(f"  frame {frame:>4}: static {rows[-1][1]:+9.1f}  "
              f"pol A {shift['A']:+9.1f}  B {shift['B']:+9.1f}  "
              f"-> {rows[-1][4]:+9.1f} cm^-1", flush=True)

    if not rows:
        raise SystemExit("no completed frames found")

    arr = np.array(rows, float)
    static, polarised = arr[:, 1], arr[:, 4]
    change = np.abs(polarised) - np.abs(static)
    n = len(arr)
    se = change.std(ddof=1) / np.sqrt(n)
    print()
    print(f"n = {n} frames")
    print(f"  mean|D| static     = {np.abs(static).mean():8.1f} cm^-1")
    print(f"  mean|D| polarised  = {np.abs(polarised).mean():8.1f} cm^-1  "
          f"({100 * (np.abs(polarised).mean() / np.abs(static).mean() - 1):+.1f}%)")
    print(f"  paired change      = {change.mean():+8.1f} +/- {se:.1f} cm^-1 "
          f"({abs(change.mean() / se):.2f} SE, "
          f"{'significant' if abs(change.mean() / se) > 2 else 'NOT significant'})")
    print(f"  sd static / polarised = {static.std(ddof=1):.1f} / {polarised.std(ddof=1):.1f}")
    print(f"  |D| > 2|J| (65.63): static {int((np.abs(static) > 65.63).sum())}/{n}, "
          f"polarised {int((np.abs(polarised) > 65.63).sum())}/{n}")
    if abs(change.mean()) > 0:
        need = int(np.ceil((2.8 * change.std(ddof=1) / abs(change.mean())) ** 2))
        print(f"  frames needed to resolve this effect at 80% power: {need}")

    if args.output:
        np.savez(args.output, frame=arr[:, 0], detuning_static_cm=static,
                 pol_shift_a_cm=arr[:, 2], pol_shift_b_cm=arr[:, 3],
                 detuning_polarised_cm=polarised)
        print(f"  wrote {args.output}")


if __name__ == "__main__":
    main()

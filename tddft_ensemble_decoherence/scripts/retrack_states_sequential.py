#!/usr/bin/env python3
"""Re-assign TDDFT roots with a temporally-continuous tracker.

THE DEFECT. select_state_by_overlap.py classifies each frame INDEPENDENTLY
against a single STATIC reference density. As the geometry evolves, every
candidate's overlap with that fixed reference degrades, and whenever two roots
come close the pick flips. On the 8 ps/4 fs campaign that gives root switches in
38% (site A) and 46% (site B) of 4 fs steps, with |dE| averaging 1047 cm^-1
across a switch against 414 cm^-1 without one. The resulting "site energy
fluctuation" is tracking noise, not dynamics -- its gap sd (1256 cm^-1 over
456 fs) exceeds the 1 ns production ensemble's (690 cm^-1), which is impossible
for a shorter window.

THE FIX. Track identity as a PATH, not a per-frame pick:
  * reference is the PREVIOUS frame's selected state, not a fixed one, so the
    reference co-evolves with the geometry (the same idea as ORCA's root homing);
  * hysteresis -- keep the incumbent root unless a challenger beats it by
    `--hysteresis`, which stops flapping between near-degenerate states;
  * an anchor check against frame 0 flags genuine drift, so a co-evolving
    reference cannot silently walk onto a different state.

Comparison uses the retained 41-atom `atom_transition_q` fingerprint rather than
the 241k-point grid: atom ordering is fixed across frames, so a cosine on the
charge vector needs no alignment and is ~4 orders cheaper. Sign is arbitrary, so
|cosine| is used, as in the original.

Recomputes NOTHING -- all seven root densities per frame are already on disk.

Writes retracked_states.json next to the manifest.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def load_roots(site_dir: Path) -> list[dict]:
    out = []
    for p in sorted(site_dir.glob("root*_density.npz")):
        d = np.load(p)
        out.append({"root": int(d["root"]), "energy_eV": float(d["energy_eV"]),
                    "q": np.asarray(d["atom_transition_q"], float)})
    return out


def cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return 0.0 if na == 0 or nb == 0 else abs(float(a @ b) / (na * nb))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", default="short")
    ap.add_argument("--hysteresis", type=float, default=0.05)
    ap.add_argument("--anchor-warn", type=float, default=0.5)
    ap.add_argument("--out", type=Path, default=ROOT / "results_retracked_states.json")
    a = ap.parse_args()

    fdir = ROOT / "frames" / a.window
    frames = sorted(fdir.glob("frame_*"), key=lambda p: int(p.name.split("_")[1]))
    EV = 8065.543937
    rec = {"site_A": [], "site_B": []}
    for site in ("site_A", "site_B"):
        prev_q = anchor_q = None
        incumbent = None
        for fr in frames:
            sd = fr / site
            if not sd.is_dir() or not any(sd.glob("root*_density.npz")):
                continue
            cands = load_roots(sd)
            if not cands:
                continue
            if prev_q is None:                       # frame 0: keep the original pick
                sel = json.loads((sd / "state_assignment.json").read_text())
                pick = next(c for c in cands if c["root"] == sel["selected_root"])
                anchor_q = pick["q"].copy()
            else:
                scored = sorted(((cos(c["q"], prev_q), c) for c in cands),
                                key=lambda t: t[0], reverse=True)
                best_s, best = scored[0]
                inc = next(((s, c) for s, c in scored if c["root"] == incumbent), None)
                # hysteresis: keep the incumbent unless clearly beaten
                pick = inc[1] if (inc and best_s - inc[0] < a.hysteresis) else best
            anch = cos(pick["q"], anchor_q)
            rec[site].append({"frame": int(fr.name.split("_")[1]),
                              "root": pick["root"],
                              "energy_cm": pick["energy_eV"] * EV,
                              "anchor_overlap": anch})
            prev_q = pick["q"].copy(); incumbent = pick["root"]

    # report
    for site in ("site_A", "site_B"):
        r = rec[site]
        roots = np.array([x["root"] for x in r]); e = np.array([x["energy_cm"] for x in r])
        anch = np.array([x["anchor_overlap"] for x in r])
        sw = int((np.diff(roots) != 0).sum())
        print(f"{site}: {len(r)} frames, {sw} root switches ({100*sw/max(1,len(r)-1):.0f}% of steps)")
        print(f"   sd(E) {e.std(ddof=1):7.1f} cm^-1   mean |dE| per 4 fs step {np.abs(np.diff(e)).mean():7.1f}")
        print(f"   anchor overlap min {anch.min():.3f} (warn below {a.anchor_warn})")
        if anch.min() < a.anchor_warn:
            print("   WARNING: reference has drifted from frame 0 -- inspect before use")
    a.out.write_text(json.dumps(rec, indent=1))
    print(f"\nwrote {a.out.name}")


if __name__ == "__main__":
    main()

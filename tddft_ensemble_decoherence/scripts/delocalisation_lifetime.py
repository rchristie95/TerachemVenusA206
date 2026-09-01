#!/usr/bin/env python3
"""How long does an FP-dimer exciton stay delocalised at 300 K?

Answered by propagating the measured two-site Hamiltonian, with no bath model,
no T2 and no Lindblad operator:

    H(t) = [[Delta(t)/2, J(t)], [J(t), -Delta(t)/2]]

Delta(t) and J(t) are the QM/MM values at 1 fs resolution -- not a
parameterisation of them. This matters because the classical-bath machinery the
manuscript currently uses has no valid regime here: modes with hbar*w < kT have
period > 160 fs and are quasi-static over any plausible coherence, while modes
fast enough to dephase lie above 834 cm^-1 where they are quantum-frozen and
produce vibronic structure. Direct propagation sidesteps all of that.

Window v6t2b is the natural experiment. Its first frame sits at |Delta| = 4
cm^-1 against 2|J| = 54 -- essentially perfect degeneracy, a genuine 50/50
exciton -- and it is the only configuration in the campaign where delocalisation
is possible at all.

Two curves are reported and they answer different questions:

  ADIABATIC  minor-site weight of the instantaneous eigenstate. What the exciton
             would look like if it followed H(t).
  PROPAGATED minor-site weight of the actual wavefunction. What it does.

hbar/2|J| = 61 fs is the timescale on which the exciton can respond; Delta sweeps
away in ~2 fs. That is the sudden limit, in which the wavefunction should be
FROZEN in its initial superposition while the Hamiltonian races off -- so the
two curves are expected to diverge, and the difference between them is the
physics.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EV = 8065.544005
HBAR_EV_FS = 0.6582119569
WINDOWS = {"v6t2a": 900, "v6t2b": 1400, "v6t2c": 1900}
SUBSTEPS = 200          # per 1 fs frame interval, with H linearly interpolated


def load(window: str):
    with (ROOT / "manifests" / f"{window}_frames.csv").open(newline="") as h:
        rows = [r for r in csv.DictReader(h) if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))
    d = np.array([float(r["site_A_energy_eV"]) - float(r["site_B_energy_eV"])
                  for r in rows]) * EV
    j = np.array([float(r["J_screened_cm-1"]) for r in rows])
    return d, j


def minor_weight(c: np.ndarray) -> float:
    p = np.abs(c) ** 2
    return float(min(p[0], p[1]) / p.sum())


def adiabatic_minor(delta: float, j: float) -> float:
    """Minor-site weight of an eigenstate of [[d/2, J], [J, -d/2]]."""
    return 0.5 * (1.0 - abs(delta) / np.sqrt(delta**2 + 4.0 * j**2))


def propagate(d_cm: np.ndarray, j_cm: np.ndarray):
    """Exact-in-each-substep propagation with H linearly interpolated."""
    n = len(d_cm)
    # start in the instantaneous eigenstate that is most delocalised at t0
    d0, j0 = d_cm[0], j_cm[0]
    h0 = np.array([[d0 / 2.0, j0], [j0, -d0 / 2.0]]) / EV
    _, vecs = np.linalg.eigh(h0)
    psi = vecs[:, 0].astype(np.complex128)
    prop_minor = [minor_weight(psi)]
    coherence = [abs(psi[0] * np.conj(psi[1]))]
    dt = 1.0 / SUBSTEPS
    for k in range(n - 1):
        for s in range(SUBSTEPS):
            f = (s + 0.5) / SUBSTEPS
            d = (1 - f) * d_cm[k] + f * d_cm[k + 1]
            j = (1 - f) * j_cm[k] + f * j_cm[k + 1]
            h = np.array([[d / 2.0, j], [j, -d / 2.0]]) / EV
            w, v = np.linalg.eigh(h)
            psi = v @ (np.exp(-1j * w * dt / HBAR_EV_FS) * (v.conj().T @ psi))
        prop_minor.append(minor_weight(psi))
        coherence.append(abs(psi[0] * np.conj(psi[1])))
    return np.array(prop_minor), np.array(coherence)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path,
                    default=ROOT / "results/delocalisation_lifetime.json")
    args = ap.parse_args()

    out = {}
    print("=== delocalisation lifetime from the measured H(t), no bath model ===\n")
    for w, seed in WINDOWS.items():
        d, j = load(w)
        adia = np.array([adiabatic_minor(x, y) for x, y in zip(d, j)])
        prop, coh = propagate(d, j)
        t = np.arange(1, len(d) + 1, dtype=float)   # frame 0 is t = 1 fs

        deloc = np.abs(d) < 2.0 * j
        first_exit = int(t[np.argmax(~deloc)]) if (~deloc).any() and not deloc.all() else None
        # time for the PROPAGATED weight to fall below half its initial value
        half = prop[0] / 2.0
        below = np.flatnonzero(prop <= half)
        t_half = float(t[below[0]]) if len(below) else None
        hbar_over_2J = HBAR_EV_FS / (2.0 * j.mean() / EV)

        out[w] = {
            "seed_frame": seed,
            "delta_t0_cm": float(abs(d[0])),
            "J_mean_cm": float(j.mean()),
            "hbar_over_2J_fs": float(hbar_over_2J),
            "adiabatic_minor_t0": float(adia[0]),
            "propagated_minor_t0": float(prop[0]),
            "frames_delocalised": int(deloc.sum()),
            "first_exit_fs": first_exit,
            "propagated_half_life_fs": t_half,
            "adiabatic_minor_mean": float(adia.mean()),
            "propagated_minor_mean": float(prop.mean()),
            "propagated_minor": prop.tolist(),
            "adiabatic_minor": adia.tolist(),
            "coherence": coh.tolist(),
        }

        print(f"  {w} (seed {seed}): |Delta|(t0) = {abs(d[0]):6.1f} cm-1, "
              f"<J> = {j.mean():.2f}, hbar/2|J| = {hbar_over_2J:.0f} fs")
        print(f"    minor-site weight at t0: adiabatic {adia[0]:.3f}, propagated {prop[0]:.3f}")
        print(f"    |Delta| < 2|J| in {int(deloc.sum())}/{len(d)} frames"
              + (f", first exit at t = {first_exit} fs" if first_exit else ""))
        print(f"    {'t (fs)':>8s} " + " ".join(f"{int(x):>6d}" for x in (1,2,3,5,10,20,50,100) if x <= t[-1]))
        idx = [int(x) - 1 for x in (1,2,3,5,10,20,50,100) if x <= t[-1]]
        print(f"    {'adiabatic':>8s} " + " ".join(f"{adia[i]:6.3f}" for i in idx))
        print(f"    {'PROPAGATED':>8s} " + " ".join(f"{prop[i]:6.3f}" for i in idx))
        print(f"    mean over the window: adiabatic {adia.mean():.3f}, "
              f"propagated {prop.mean():.3f}  (ratio {prop.mean()/max(adia.mean(),1e-9):.1f}x)")
        print()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

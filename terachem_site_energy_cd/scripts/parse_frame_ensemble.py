#!/usr/bin/env python3
"""Extract a per-frame exciton ensemble from completed TeraChem site directories.

For every frame directory holding converged site A and site B jobs, pull the
bright-state energy, the transition-dipole VECTOR, and the chromophore centroid,
so the pair can be fed to a vibronic/CD model without any surrogate sampling.

Bright-state selection is by maximum oscillator strength, with the runner-up
margin reported; frames whose margin falls below --min-margin are dropped rather
than averaged, because a fragmented bright state makes the site energy
ill-defined by up to ~1000 cm^-1.

Dipoles and centroids come from each site's own geometry.xyz, which the frame
preparer writes in imaged lab coordinates, so both sites already share a frame
and the triple product R_AB . (mu_A x mu_B) is meaningful. The centroid
separation is checked against a plausible range to catch periodic-image errors.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

CM_PER_EV = 8065.54
REPO = Path(__file__).resolve().parents[2]
PARSER = Path(__file__).resolve().parent / "parse_terachem_tddft.py"
CR2_ATOM_COUNT = 29          # geometry.xyz layout: CR2(29) + Tyr phenol(12) + caps(3)

# Excited-state phase is arbitrary per job, so a raw transition dipole may come
# back as +mu or -mu at random. Left unfixed the chirality triple product
# R_AB . (mu_A x mu_B) flips sign frame to frame and the ensemble CD couplet
# averages to zero. Convention: orient every dipole along the chromophore long
# axis, from the imidazolinone ring centroid toward the phenolate oxygen.
# Indices into the 29-atom CR2 block (see cr2_atom_names in the STEOM probe):
IMIDAZOLINONE_RING = (3, 4, 5, 7, 8)     # N2, CA2, C2, N3, CA3
PHENOLATE_OXYGEN = 22                    # OH


def _bright(site_dir: Path, python: str) -> dict | None:
    out = site_dir / "tddft.out"
    if not out.is_file() or "Job finished" not in out.read_text(errors="ignore"):
        return None
    summary = site_dir / "energy_summary.json"
    if not summary.is_file():
        subprocess.run([python, str(PARSER), str(out), "--json", str(summary)],
                       capture_output=True)
    if not summary.is_file():
        return None
    roots = json.loads(summary.read_text())["roots"]
    ranked = sorted(roots, key=lambda r: -r["oscillator_strength"])
    best, runner = ranked[0], ranked[1]
    if best["transition_dipole_au"] is None:
        return None
    return {
        "energy_cm": best["energy_eV"] * CM_PER_EV,
        "mu_au": np.asarray(best["transition_dipole_au"], float),
        "oscillator_strength": best["oscillator_strength"],
        "margin": best["oscillator_strength"] / max(runner["oscillator_strength"], 1e-12),
        "root": best["root"],
    }


def _cr2_coords(site_dir: Path) -> np.ndarray:
    lines = (site_dir / "geometry.xyz").read_text().splitlines()
    return np.array([[float(v) for v in ln.split()[1:4]]
                     for ln in lines[2:2 + CR2_ATOM_COUNT]])


def _long_axis(cr2: np.ndarray) -> np.ndarray:
    """Unit vector from the imidazolinone ring centroid to the phenolate oxygen."""
    axis = cr2[PHENOLATE_OXYGEN] - cr2[list(IMIDAZOLINONE_RING)].mean(axis=0)
    return axis / np.linalg.norm(axis)


def _phase_align(mu: np.ndarray, cr2: np.ndarray) -> tuple[np.ndarray, float]:
    """Fix the arbitrary excited-state phase; return (oriented mu, |cos angle|)."""
    axis = _long_axis(cr2)
    projection = float(np.dot(mu, axis))
    cosine = abs(projection) / np.linalg.norm(mu)
    return (mu if projection >= 0.0 else -mu), cosine


def collect(results_dir: Path, prefix: str, python: str, min_margin: float = 3.0,
            separation_range=(15.0, 40.0)) -> dict:
    frames, dropped, suspect = [], [], []
    for frame_dir in sorted(results_dir.glob(f"{prefix}_frame_*")):
        if any(tag in frame_dir.name for tag in ("diffuse", "gas", "selection")):
            continue
        a = _bright(frame_dir / "site_A", python)
        b = _bright(frame_dir / "site_B", python)
        if a is None or b is None:
            continue
        margin = min(a["margin"], b["margin"])
        if margin < min_margin:
            dropped.append((frame_dir.name, round(margin, 2)))
            continue
        cr2_a = _cr2_coords(frame_dir / "site_A")
        cr2_b = _cr2_coords(frame_dir / "site_B")
        a["mu_au"], cos_a = _phase_align(a["mu_au"], cr2_a)
        b["mu_au"], cos_b = _phase_align(b["mu_au"], cr2_b)
        r_a, r_b = cr2_a.mean(axis=0), cr2_b.mean(axis=0)
        separation = float(np.linalg.norm(r_b - r_a))
        if not separation_range[0] <= separation <= separation_range[1]:
            suspect.append((frame_dir.name, round(separation, 2)))
            continue
        frames.append({
            "frame": int("".join(ch for ch in frame_dir.name if ch.isdigit())[-4:]),
            "e_a_cm": a["energy_cm"], "e_b_cm": b["energy_cm"],
            "mu_a_au": a["mu_au"], "mu_b_au": b["mu_au"],
            "r_a_ang": r_a, "r_b_ang": r_b,
            "separation_ang": separation,
            "triple_product": float(np.dot(r_b - r_a, np.cross(a["mu_au"], b["mu_au"]))),
            "margin": margin,
            "axis_cosine": min(cos_a, cos_b),
        })
    return {"frames": frames, "dropped_low_margin": dropped,
            "dropped_bad_separation": suspect}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path,
                    default=REPO / "terachem_site_energy_cd/results")
    ap.add_argument("--prefix", default="hi_camb3lyp")
    ap.add_argument("--min-margin", type=float, default=3.0)
    ap.add_argument("--python", default="/home/robson/anaconda3/envs/TeraChem/bin/python")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    data = collect(args.results_dir, args.prefix, args.python, args.min_margin)
    frames = data["frames"]
    if not frames:
        raise SystemExit(f"no usable frames for prefix {args.prefix!r}")

    detuning = np.array([f["e_a_cm"] - f["e_b_cm"] for f in frames])
    separation = np.array([f["separation_ang"] for f in frames])
    triple = np.array([f["triple_product"] for f in frames])
    print(f"prefix {args.prefix}: {len(frames)} usable frames "
          f"({len(data['dropped_low_margin'])} dropped low-margin, "
          f"{len(data['dropped_bad_separation'])} dropped bad separation)")
    print(f"  detuning   mean {detuning.mean():+8.1f}  sd {detuning.std(ddof=1):7.1f}  "
          f"mean|D| {np.abs(detuning).mean():7.1f} cm^-1")
    print(f"  separation mean {separation.mean():7.2f} +/- {separation.std(ddof=1):.2f} A")
    cosines = np.array([f["axis_cosine"] for f in frames])
    print(f"  triple product sign: {int((triple > 0).sum())} positive / "
          f"{int((triple < 0).sum())} negative  (a consistent sign is required "
          f"for a well-defined couplet handedness)")
    print(f"  phase-alignment |cos| to chromophore long axis: "
          f"min {cosines.min():.3f}  mean {cosines.mean():.3f}  "
          f"(near 0 would make the sign convention ambiguous)")
    if data["dropped_low_margin"]:
        print(f"  low-margin frames: {data['dropped_low_margin']}")
    if data["dropped_bad_separation"]:
        print(f"  suspect separations: {data['dropped_bad_separation']}")

    if args.output:
        np.savez(args.output,
                 frame=np.array([f["frame"] for f in frames]),
                 e_a_cm=np.array([f["e_a_cm"] for f in frames]),
                 e_b_cm=np.array([f["e_b_cm"] for f in frames]),
                 mu_a_au=np.array([f["mu_a_au"] for f in frames]),
                 mu_b_au=np.array([f["mu_b_au"] for f in frames]),
                 r_a_ang=np.array([f["r_a_ang"] for f in frames]),
                 r_b_ang=np.array([f["r_b_ang"] for f in frames]),
                 triple_product=triple, separation_ang=separation,
                 margin=np.array([f["margin"] for f in frames]))
        print(f"  wrote {args.output}")


if __name__ == "__main__":
    main()

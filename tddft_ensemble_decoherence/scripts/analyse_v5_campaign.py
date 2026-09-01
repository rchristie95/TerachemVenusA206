#!/usr/bin/env python3
"""Merge the three v5 worker manifests and report the campaign result.

The three workers hold disjoint interleaved slices of the same 500-row manifest
(save_manifest rewrites the whole csv, so they cannot share one). This merges
them back into manifests/v5camp_frames.csv in frame order, then reports:

  * the coupling and detuning ensembles, against v4full and the rigid-fit value
    the campaign exists to replace;
  * the STEOM -> TDDFT site-energy calibration from the four matched frames,
    which is the point of running STEOM at all: it tests whether the single
    constant shift (DELTA_STEOM_EV) the pipeline applies is defensible, or
    whether the shift tracks the methine-bridge coordinate that carries most of
    the site-energy variance.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
MANIFESTS = ROOT / "manifests"
STEOM = ROOT / "steom_v5"
DELTA_STEOM_EV = -1.06778287       # the single constant the pipeline applies
EV_CM = 8065.543937
RIGID_FIT_J = 32.82                # static-density production value being replaced


def merge() -> tuple[list[str], list[dict]]:
    rows: list[dict] = []
    cols: list[str] | None = None
    for tag in "abc":
        p = MANIFESTS / f"v5camp_{tag}_frames.csv"
        with p.open() as fh:
            rd = csv.DictReader(fh)
            cols = cols or list(rd.fieldnames or [])
            rows.extend(rd)
    rows.sort(key=lambda r: int(r["frame_index"]))
    out = MANIFESTS / "v5camp_frames.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return cols, rows


def stats(name: str, x: np.ndarray) -> str:
    if x.size == 0:
        return f"  {name:<28s} (no data)"
    sem = x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else float("nan")
    return (f"  {name:<28s} {x.mean():9.3f} +/- {sem:6.3f}   "
            f"sd {x.std(ddof=1) if x.size > 1 else float('nan'):8.3f}   n={x.size}")


def ensemble(rows: list[dict]) -> dict:
    done = [r for r in rows if r["convergence_status"] == "complete"]
    if not done:
        print("no completed frames yet")
        return {}
    J = np.array([float(r["J_screened_cm-1"]) for r in done])
    Jv = np.array([float(r["J_vacuum_cm-1"]) for r in done])
    Jp = np.array([float(r["J_pda_screened_cm-1"]) for r in done])
    EA = np.array([float(r["site_A_energy_cm-1"]) for r in done])
    EB = np.array([float(r["site_B_energy_cm-1"]) for r in done])
    D = EA - EB
    print(f"\n=== v5 campaign ensemble, {len(done)}/{len(rows)} frames ===")
    print(stats("J screened (cm^-1)", J))
    print(stats("J vacuum   (cm^-1)", Jv))
    print(stats("J PDA scr. (cm^-1)", Jp))
    print(stats("|Delta|    (cm^-1)", np.abs(D)))
    print(stats("Delta      (cm^-1)", D))
    print(f"  {'TDC/PDA ratio':<28s} {np.mean(J / Jp):9.3f}")
    print(f"\n  rigid-fit static density   J = {RIGID_FIT_J:.2f} cm^-1 "
          f"-> per-frame recalculation moves it by "
          f"{J.mean() - RIGID_FIT_J:+.2f} ({100 * (J.mean() / RIGID_FIT_J - 1):+.1f} %)")

    # Two-state delocalisation on the INSTANTANEOUS gap, reported per frame
    # rather than from mean J against mean |D|. The measure is the weight the
    # eigenvector puts on its minor site,
    #     w = (1/2) (1 - |D| / sqrt(D^2 + 4 J^2)),
    # which runs 0 (fully localised) to 1/2 (fully delocalised). Do NOT use
    # sin^2(theta) with theta = arctan2(2J, D): arctan2 returns ~pi whenever D is
    # negative, so that expression reports ~1 for a strongly localised frame
    # purely because site B happened to be the higher of the two.
    mix = 0.5 * (1.0 - np.abs(D) / np.sqrt(D ** 2 + 4 * J ** 2))
    print(f"  {'minor-site weight':<28s} {mix.mean():9.5f}  "
          f"median {np.median(mix):.5f}  max {mix.max():.5f}   (1/2 = full)")
    print(f"  {'frames with |Delta| < 2|J|':<28s} {int(np.sum(np.abs(D) < 2 * J)):9d}"
          f"  ({100 * np.mean(np.abs(D) < 2 * J):.1f} %)")
    return {"n": len(done), "J_mean": float(J.mean()),
            "J_sem": float(J.std(ddof=1) / np.sqrt(J.size)),
            "J_sd": float(J.std(ddof=1)),
            "abs_delta_mean": float(np.abs(D).mean()),
            "abs_delta_sd": float(np.abs(D).std(ddof=1)),
            "mixing_mean": float(mix.mean()), "mixing_max": float(mix.max())}


def steom_calibration(rows: list[dict]) -> dict:
    """TDDFT -> STEOM site-energy shift at the four matched frames."""
    by_frame = {int(r["frame_index"]): r for r in rows
                if r["convergence_status"] == "complete"}
    pts = []
    for d in sorted(STEOM.glob("frame_*")):
        frame = int(d.name.split("_")[1])
        for site in "AB":
            res = d / f"site_{site}" / "steom_result.json"
            if not res.exists():
                continue
            j = json.loads(res.read_text())
            if j.get("status") != "complete" or "energy_eV" not in j:
                continue
            row = by_frame.get(frame)
            if row is None:
                continue
            td = float(row[f"site_{site}_energy_eV"])
            pts.append({"frame": frame, "site": site, "tddft_eV": td,
                        "steom_eV": j["energy_eV"], "shift_eV": j["energy_eV"] - td,
                        "steom_mu_D": j.get("mu_D"), "steom_f": j.get("f")})
    if not pts:
        print("\n=== STEOM calibration: no matched points yet ===")
        return {}
    s = np.array([p["shift_eV"] for p in pts])
    print(f"\n=== STEOM -> TDDFT calibration, {len(pts)} matched site-frames ===")
    print(f"  {'frame':>6s} {'site':>5s} {'TDDFT eV':>10s} {'STEOM eV':>10s} "
          f"{'shift eV':>10s} {'|mu| D':>8s}")
    for p in sorted(pts, key=lambda q: (q["frame"], q["site"])):
        mu = p["steom_mu_D"]
        print(f"  {p['frame']:6d} {p['site']:>5s} {p['tddft_eV']:10.4f} "
              f"{p['steom_eV']:10.4f} {p['shift_eV']:+10.4f} "
              f"{(f'{mu:8.3f}' if mu else '       -')}")
    print(stats("shift (eV)", s))
    print(f"  spread of the shift        {s.std(ddof=1) * EV_CM:9.1f} cm^-1"
          if s.size > 1 else "")
    print(f"  pipeline constant          {DELTA_STEOM_EV:+.5f} eV; "
          f"measured mean {s.mean():+.5f} eV "
          f"({(s.mean() - DELTA_STEOM_EV) * EV_CM:+.1f} cm^-1 from it)")
    if s.size > 2:
        # The shift is only usable as a constant if it does not track the site
        # energy itself; a slope means the STEOM correction is state-dependent.
        td = np.array([p["tddft_eV"] for p in pts])
        A = np.vstack([td - td.mean(), np.ones_like(td)]).T
        slope, icpt = np.linalg.lstsq(A, s, rcond=None)[0]
        resid = s - (A @ [slope, icpt])
        print(f"  shift vs TDDFT energy      slope {slope:+.4f} eV/eV, "
              f"residual sd {resid.std(ddof=1) * EV_CM:.1f} cm^-1")
    return {"n": len(pts), "shift_mean_eV": float(s.mean()),
            "shift_sd_eV": float(s.std(ddof=1)) if s.size > 1 else None,
            "pipeline_constant_eV": DELTA_STEOM_EV, "points": pts}


def emit_coupling_samples(rows: list[dict], out_dir: Path) -> None:
    """Write the campaign in the shape make_tandem_panels.py already reads.

    That generator (arXiv-2605.00027v1/make_tandem_panels.py) takes a
    coupling_samples.csv with `frame`, `J_cm`, `J_pda_cm`, `angle_deg`,
    `separation_A` plus a coupling_distribution.json carrying `mean`/`std`, and
    is the source of the Figure 5 panels that FIGURE_DATA_AUDIT.md flags as
    built by the superseded rigid-fit route. Emitting this shape lets the figure
    be regenerated from per-frame densities without touching the generator.

    Two deliberate differences from the rigid-fit file:

      * `frame` is the 0-based production index, not the DCD frame number. The
        generator computes its time axis as `frame / frame.max()`, so feeding it
        DCD numbers 900-2896 would compress 1 ns of trajectory into the last
        69% of the axis.
      * `aln_A_rms` / `aln_B_rms` are absent. They recorded how well one static
        density was fitted onto each frame, and nothing is fitted here -- the
        density is recomputed. The generator reads only the columns it uses, so
        their absence is not an error, and inventing a value would be.
    """
    done = sorted((r for r in rows if r["convergence_status"] == "complete"),
                  key=lambda r: int(r["frame_index"]))
    if not done:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "coupling_samples.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "J_cm", "J_pda_cm", "angle_deg", "separation_A"])
        for i, r in enumerate(done):
            w.writerow([i, r["J_screened_cm-1"], r["J_pda_screened_cm-1"],
                        r["dipole_angle_deg"], r["separation_A"]])
    J = np.array([float(r["J_screened_cm-1"]) for r in done])
    (out_dir / "coupling_distribution.json").write_text(json.dumps({
        "mean": float(J.mean()),
        "std": float(J.std(ddof=1)) if J.size > 1 else 0.0,
        "n": int(J.size),
        "source": "v5 campaign, per-frame transition densities, full-system embedding",
    }, indent=2))
    print(f"  wrote {csv_path} and coupling_distribution.json ({len(done)} frames)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", type=Path, default=ROOT / "results/v5camp_summary.json")
    ap.add_argument(
        "--emit-coupling-samples",
        type=Path,
        nargs="?",
        const=ROOT / "results/coupling_tandem_v5",
        help="Also write coupling_samples.csv / coupling_distribution.json for "
             "make_tandem_panels.py (Figure 5).",
    )
    args = ap.parse_args()
    _, rows = merge()
    out = {"ensemble": ensemble(rows), "steom": steom_calibration(rows)}
    if args.emit_coupling_samples:
        print()
        emit_coupling_samples(rows, args.emit_coupling_samples)
    args.json.parent.mkdir(exist_ok=True)
    args.json.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()

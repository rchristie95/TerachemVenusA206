#!/usr/bin/env python3
"""Build the TDX control for an existing TD site-energy frame.

Nguyen's TDX construct deletes the tyrosine of one chromophore-forming triad,
so the partner beta-barrel is still present and docked but its chromophore is
absent -- and, decisively, the partner is no longer an anionic residue carrying
-1 e at ~25 A. The measured TD-minus-TDX absorption red shift (35.3 cm^-1 for
dVenus) is therefore the sum of that electrostatic change and any excitonic
shift, and this script isolates the electrostatic part.

Everything except the partner CR2's 29 MM point charges is carried over
byte-for-byte from the completed TD calculation: the same QM geometry, the same
TeraChem input, the same 69,374-charge embedding field at the same coordinates.
Only the charge column of the partner-CR2 block changes, so the difference in
excitation energy is attributable to that block alone.

Two models are provided:

  neutral  add +1/29 e to each partner-CR2 atom, so the residue is neutral but
           keeps its charge distribution shape.  This is the faithful TDX
           model: the residue is still there, it is simply not an anion.
  zero     set all 29 charges to zero, removing the monopole and every higher
           multipole.  An upper bound, and the difference between the two
           models separates the monopole term from the rest.

The partner CR2 is located geometrically -- by matching the opposite site's QM
chromophore coordinates into this site's MM field -- rather than by trusting an
index, and the match is checked to be exactly 29 contiguous rows summing to
-1 e before anything is written.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

OPPOSITE = {"A": "B", "B": "A"}
CR2_ATOMS = 29


def read_geometry(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    return np.array(
        [[float(x) for x in line.split()[1:4]] for line in lines[2 : 2 + count]]
    )


def read_field(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    header = lines[1]
    rows = [line.split() for line in lines[2 : 2 + count]]
    charges = np.array([float(r[0]) for r in rows])
    coords = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
    return charges, coords, header


def write_field(path: Path, charges, coords, header: str) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(charges)}\n")
        handle.write(f"{header}\n")
        for q, c in zip(charges, coords):
            handle.write(f"{q: .8f} {c[0]: .8f} {c[1]: .8f} {c[2]: .8f}\n")


def locate_partner_cr2(field_coords: np.ndarray, partner_cr2: np.ndarray) -> list[int]:
    """Rows of this site's MM field that are the partner chromophore's atoms."""
    tree = cKDTree(field_coords)
    found = set()
    for point in partner_cr2:
        found.update(tree.query_ball_point(point, 0.05))
    index = sorted(found)
    if len(index) != CR2_ATOMS:
        raise RuntimeError(
            f"matched {len(index)} MM rows to the partner chromophore, expected {CR2_ATOMS}"
        )
    if index != list(range(index[0], index[-1] + 1)):
        raise RuntimeError("partner-CR2 MM rows are not contiguous")
    return index


def build_site(td_site: Path, tdx_site: Path, partner_cr2: np.ndarray, model: str) -> dict:
    tdx_site.mkdir(parents=True, exist_ok=True)
    for name in ("geometry.xyz", "tddft.in"):
        shutil.copy2(td_site / name, tdx_site / name)

    charges, coords, header = read_field(td_site / "mm_charges.dat")
    index = locate_partner_cr2(coords, partner_cr2)
    before = float(charges[index].sum())
    if abs(before + 1.0) > 1e-6:
        raise RuntimeError(f"partner CR2 net charge is {before:.6f} e, expected -1")

    modified = charges.copy()
    if model == "neutral":
        modified[index] += 1.0 / CR2_ATOMS
    elif model == "zero":
        modified[index] = 0.0
    else:
        raise ValueError(model)
    after = float(modified[index].sum())
    if abs(after) > 1e-9:
        raise RuntimeError(f"partner CR2 net charge after modification is {after:.3e} e")

    untouched = np.ones(len(charges), bool)
    untouched[index] = False
    if not np.array_equal(charges[untouched], modified[untouched]):
        raise RuntimeError("charges outside the partner chromophore changed")

    write_field(
        tdx_site / "mm_charges.dat",
        modified,
        coords,
        f"{header}; TDX control: partner CR2 charges {model}"
        f" (was {before:+.6f} e, now {after:+.6f} e)",
    )
    return {
        "partner_cr2_rows": [index[0], index[-1]],
        "net_charge_before_e": before,
        "net_charge_after_e": after,
        "model": model,
        "separation_A": float(
            np.linalg.norm(coords[index].mean(0) - read_geometry(td_site / "geometry.xyz")[:CR2_ATOMS].mean(0))
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("td_frame_dir", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--model", choices=("neutral", "zero"), default="neutral")
    args = ap.parse_args()

    record = {"td_frame_dir": str(args.td_frame_dir), "model": args.model, "sites": {}}
    for site in ("A", "B"):
        td_site = args.td_frame_dir / f"site_{site}"
        if not (td_site / "mm_charges.dat").is_file():
            raise SystemExit(f"missing TD baseline at {td_site}")
        partner = read_geometry(args.td_frame_dir / f"site_{OPPOSITE[site]}" / "geometry.xyz")[:CR2_ATOMS]
        record["sites"][site] = build_site(
            td_site, args.output_dir / f"site_{site}", partner, args.model
        )
    (args.output_dir / "tdx_preparation.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the 88-atom two-chromophore QM region for the supermolecular J test.

Every coupling in this project is a Coulomb integral between two SEPARATELY
computed monomer transition densities. This emits the inputs for the one test
that can fail independently of every unit convention: put both chromophores in a
single QM region and read the excitonic splitting out of the eigenvalues.

Four job directories per frame, all sharing one geometry construction so nothing
differs between them except what is meant to:

    super_vac_true   88 atoms, charge -2, NO point charges, true geometry
    super_vac_far    the identical 88 atoms with chromophore B translated 200 A
    super_emb        88 atoms, charge -2, MM field excluding BOTH regions
    monoA_emb        44 atoms, charge -1, that same field, partner ABSENT
    monoB_emb        44 atoms, charge -1, that same field, partner ABSENT

`super_vac_far` is the point of the design. It is the non-interacting reference
computed through the identical code path and basis, so its two bright roots ARE
the isolated monomer energies: Delta_vac comes free, basis-set superposition
cancels by construction, and whatever splitting survives at 200 A beyond
Delta_vac IS the numerical floor. At frame 1832 the predicted signal is
0.006087 eV, so the floor has to be measured before any result can be claimed.

The monomer references deliberately do NOT follow the production convention of
retaining the partner CR2 at -1 e in the field. In the supermolecule the partner
is QM, so the monomer reference must have it absent entirely or the identity
Omega = sqrt(Delta^2 + 4J^2) cannot hold. Delta from these will therefore not
equal the production 0.17 cm^-1, and that is correct.

Atom ordering is 82 physical then 6 link for the supermolecule (41_A, 41_B, then
3_A, 3_B), keeping link atoms last as compact_tddft_density.py requires.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from openmm.app import PDBFile
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_production_frame import (                             # noqa: E402
    SITE_KEYS,
    cap_along_severed_bond,
    image_near,
    load_or_build_charges,
    read_dcd_frame,
    reference_mapping,
    residue_by_key,
    severed_bond,
    unwrap,
    write_field,
)

FAR_TRANSLATION_A = 200.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--topology", type=Path, required=True)
    p.add_argument("--trajectory", type=Path, required=True)
    p.add_argument("--frame-index", type=int, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--embedding-cache", type=Path, required=True)
    p.add_argument("--site-b-offset", type=int, default=0)
    return p.parse_args()


def build_site(pdb, positions, box, keys, cr2_names, cr2_symbols,
               tyr_names, tyr_symbols):
    """One 44-atom region: 41 physical coords/symbols, 3 caps, and its indices."""
    cr2_residue = residue_by_key(pdb.topology, keys["cr2"], "CR2")
    tyr_residue = residue_by_key(pdb.topology, keys["tyr"], "TYR")
    cr2_by_name = {a.name: a for a in cr2_residue.atoms()}
    if set(cr2_by_name) != set(cr2_names):
        raise RuntimeError("CR2 atom inventory mismatch")
    cr2_coords = unwrap(
        np.asarray([positions[cr2_by_name[n].index] for n in cr2_names]), box)

    tyr_by_name = {a.name: a for a in tyr_residue.atoms()}
    tyr_atoms = list(tyr_residue.atoms())
    tyr_all = image_near(
        np.asarray([positions[a.index] for a in tyr_atoms]),
        cr2_coords.mean(axis=0), box)
    tyr_coord_by_name = dict(zip((a.name for a in tyr_atoms), tyr_all))
    tyr_coords = np.asarray([tyr_coord_by_name[n] for n in tyr_names])

    qm_indices = {a.index for a in cr2_residue.atoms()}
    qm_indices |= {tyr_by_name[n].index for n in tyr_names}
    cuts = (
        (severed_bond(pdb.topology, qm_indices, "N1", cr2_residue),
         1.01, cr2_coords[cr2_names.index("N1")]),
        (severed_bond(pdb.topology, qm_indices, "C3", cr2_residue),
         1.09, cr2_coords[cr2_names.index("C3")]),
        (severed_bond(pdb.topology, qm_indices, "CG", tyr_residue),
         1.09, tyr_coord_by_name["CG"]),
    )
    caps, diagnostics, partners = [], [], set()
    for cut, length, parent in cuts:
        cap, diag = cap_along_severed_bond(*cut, positions, box, length, parent)
        caps.append(cap)
        diagnostics.append(diag)
        partners.add(cut[1].index)
    return {
        "physical_coords": np.vstack([cr2_coords, tyr_coords]),
        "physical_symbols": list(cr2_symbols) + list(tyr_symbols),
        "cap_coords": np.asarray(caps),
        "qm_indices": qm_indices,
        "link_partner_indices": partners,
        "link_caps": diagnostics,
    }


def write_xyz(path: Path, symbols, coords, comment: str) -> None:
    with path.open("w") as handle:
        handle.write(f"{len(symbols)}\n{comment}\n")
        for s, c in zip(symbols, coords):
            handle.write(f"{s:<2s} {c[0]: .10f} {c[1]: .10f} {c[2]: .10f}\n")


def build_field(pdb, positions, box, charges, excluded, anchor):
    """MM field excluding `excluded` atom indices, imaged near `anchor`."""
    keep_idx, keep_xyz = [], []
    counts = {"protein": 0, "water": 0, "ions": 0}
    for residue in pdb.topology.residues():
        atoms = [a for a in residue.atoms() if a.index not in excluded]
        if not atoms:
            continue
        all_atoms = list(residue.atoms())
        imaged = image_near(
            np.asarray([positions[a.index] for a in all_atoms]), anchor, box)
        by_index = dict(zip((a.index for a in all_atoms), imaged))
        for atom in atoms:
            keep_idx.append(atom.index)
            keep_xyz.append(by_index[atom.index])
        name = residue.name
        key = ("water" if name in ("HOH", "WAT")
               else "ions" if len(all_atoms) == 1 else "protein")
        counts[key] += len(atoms)
    return np.asarray(keep_idx, dtype=int), np.asarray(keep_xyz), counts


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdb = PDBFile(str(args.topology))
    positions, box, dcd_meta = read_dcd_frame(args.trajectory, args.frame_index)
    charges = load_or_build_charges(pdb, args.embedding_cache)
    cr2_names, cr2_symbols, _, tyr_names, tyr_symbols = reference_mapping()

    def offset_keys(site):
        out = {}
        for role, value in SITE_KEYS[site].items():
            if role == "kind":          # not a residue key; filter before unpack
                continue
            chain, num = value
            out[role] = (chain,
                         str(int(num) + (args.site_b_offset if site == "B" else 0)))
        return out

    keys = {s: offset_keys(s) for s in ("A", "B")}
    sites = {s: build_site(pdb, positions, box, keys[s], cr2_names,
                           cr2_symbols, tyr_names, tyr_symbols)
             for s in ("A", "B")}

    # ---- supermolecule geometry: 82 physical then 6 link -------------------
    a, b = sites["A"], sites["B"]
    # image B near A so the pair is contiguous rather than split by the box
    shift = np.round((a["physical_coords"].mean(0) - b["physical_coords"].mean(0))
                     / box) * box
    for key in ("physical_coords", "cap_coords"):
        b[key] = b[key] + shift
    super_symbols = (a["physical_symbols"] + b["physical_symbols"]
                     + ["H"] * 3 + ["H"] * 3)
    super_coords = np.vstack([a["physical_coords"], b["physical_coords"],
                              a["cap_coords"], b["cap_coords"]])
    separation = float(np.linalg.norm(a["physical_coords"].mean(0)
                                      - b["physical_coords"].mean(0)))

    far_coords = super_coords.copy()
    n_a, n_b = len(a["physical_symbols"]), len(b["physical_symbols"])
    far_coords[n_a:n_a + n_b] += FAR_TRANSLATION_A          # B physical
    far_coords[-3:] += FAR_TRANSLATION_A                    # B caps

    for tag, coords in (("super_vac_true", super_coords),
                        ("super_vac_far", far_coords)):
        d = args.output_dir / tag
        d.mkdir(exist_ok=True)
        write_xyz(d / "geometry.xyz", super_symbols, coords,
                  f"frame {args.frame_index} {tag}; 82 physical + 6 link; "
                  f"charge -2; vacuum")

    # ---- embedded jobs -----------------------------------------------------
    both_qm = a["qm_indices"] | b["qm_indices"]
    both_partners = a["link_partner_indices"] | b["link_partner_indices"]
    anchor = super_coords[:n_a + n_b].mean(axis=0)

    jobs = {
        "super_emb": (super_symbols, super_coords, both_qm | both_partners, -2),
        "monoA_emb": (a["physical_symbols"] + ["H"] * 3,
                      np.vstack([a["physical_coords"], a["cap_coords"]]),
                      both_qm | both_partners, -1),
        "monoB_emb": (b["physical_symbols"] + ["H"] * 3,
                      np.vstack([b["physical_coords"], b["cap_coords"]]),
                      both_qm | both_partners, -1),
    }
    summary = {
        "frame_index": args.frame_index,
        "topology": str(args.topology.resolve()),
        "trajectory": str(args.trajectory.resolve()),
        "box_A": [float(v) for v in box],
        "dcd": dcd_meta,
        "centroid_separation_A": separation,
        "far_translation_A": FAR_TRANSLATION_A,
        "site_keys": {s: {r: list(v) for r, v in keys[s].items()} for s in keys},
        "link_caps": {s: sites[s]["link_caps"] for s in sites},
        "jobs": {},
    }
    for tag, (symbols, coords, excluded, charge) in jobs.items():
        d = args.output_dir / tag
        d.mkdir(exist_ok=True)
        # Both monomer jobs exclude BOTH QM regions from the field: the partner
        # is absent entirely, matching the supermolecule where it is QM. This is
        # NOT the production convention of retaining it at -1 e.
        idx, xyz, counts = build_field(pdb, positions, box, charges, excluded,
                                       anchor)
        q = charges[idx]
        write_field(d / "mm_charges.dat", q, xyz)
        write_xyz(d / "geometry.xyz", symbols, coords,
                  f"frame {args.frame_index} {tag}; charge {charge}")
        summary["jobs"][tag] = {
            "qm_atoms": len(symbols),
            "qm_net_charge_e": charge,
            "mm_point_charges": int(len(q)),
            "mm_net_charge_e": float(q.sum()),
            "mm_component_atom_counts": counts,
            "partner_present_in_field": False,
        }
    for tag in ("super_vac_true", "super_vac_far"):
        summary["jobs"][tag] = {"qm_atoms": len(super_symbols),
                                "qm_net_charge_e": -2,
                                "mm_point_charges": 0,
                                "vacuum": True}
    (args.output_dir / "preparation.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(f"[super] frame {args.frame_index}: {len(super_symbols)} QM atoms, "
          f"centroid separation {separation:.2f} A")
    for tag, info in summary["jobs"].items():
        extra = ("vacuum" if info.get("vacuum")
                 else f"{info['mm_point_charges']} charges, "
                      f"net {info['mm_net_charge_e']:+.6f} e")
        print(f"  {tag:<14} {info['qm_atoms']:3d} atoms, charge "
              f"{info['qm_net_charge_e']:+d}, {extra}")


if __name__ == "__main__":
    main()

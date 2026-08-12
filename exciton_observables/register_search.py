#!/usr/bin/env python3
"""Does ANY rigid dimer arrangement satisfy all the experimental constraints?

The joint analysis leaves one structural question: two independent observables
(the limiting anisotropy, |cos alpha| = 0.660 +/- 0.061, and the excitonic part
of the absorption red shift, 0.600 +/- 0.125) agree on alpha ~ 130 deg, while
every geometry we have -- crystal 110.4 deg, MD 102 deg -- sits ~20 deg short.
notes/lattice_dimer_scan.py only ever tested the six crystallographic operators,
a vanishing subset of the arrangements a tethered tandem can adopt.

This searches the full rigid-body space. Each trial is built so that the two
chromophore separation and the inter-dipole angle are satisfied by construction,
then filtered on sterics and scored on the remaining observables:

  alpha           125-136 deg     (anisotropy + red shift)
  separation      23-29 A         (sets J, and the crystal sits at 25.4)
  J               20-50 cm^-1     (transition-density ensemble is 32.8 +/- 1.6)
  triple product  negative        (CD couplet handedness, every frame so far)
  sterics         no heavy atom below 2.4 A (the biological dimer has none)
  interface       >= 40 atoms in the [2.4, 4.5) A shell. CALIBRATED against the
                  biological dimer, which scores 72; a lattice brush scores 22.
  linker span     <= 60 A         Cterm(A res 230) -> Nterm(B res 0), against
                  54.2 A in the crystal register, which the filter must admit

A hit would be a new structural model for the tandem that reconciles the
geometry with experiment. No hit, over a dense search, is itself a result: it
would mean no rigid two-barrel arrangement can do it, and the discrepancy has
to live in the photophysical model rather than the structure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_CM = 219474.63
EPSILON = 1.77
TDC_OVER_PDA = 1.1872

PDB = Path("/home/robson/PetaChem/1MYW.pdb")
DENSITY = Path("/home/robson/PetaChem/neo_model/orca_steom/steom_transdens_specnorm_oldframe.npz")
MONOMER = Path("/home/robson/PetaChem/tc_simple_old/classical_relaxed.pdb")

CLASH_A = 2.4   # calibrated: the 1MYW biological dimer has 0 heavy atoms below 2.4 A
CONTACT_A = 4.5  # shell [CLASH_A, CONTACT_A): 72 for the biological dimer, 22 for a lattice brush
VOXEL_A = 0.6


def kabsch(source, target):
    sc, tc = source - source.mean(0), target - target.mean(0)
    u, _, vt = np.linalg.svd(sc.T @ tc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, target.mean(0) - rot @ source.mean(0)


def read_chain_a():
    heavy, cr2, ca = [], {}, {}
    for line in open(PDB):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if line[17:20].strip() == "HOH" or line[76:78].strip().upper() == "H":
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        heavy.append(xyz)
        name, resi = line[12:16].strip(), int(line[22:26])
        if line[17:20].strip() == "CR2":
            cr2[name] = np.array(xyz)
        if name == "CA":
            ca[resi] = np.array(xyz)
    return np.array(heavy), cr2, ca


def cr2_reference():
    atoms = {}
    for line in open(MONOMER):
        if line[:6].strip() in ("ATOM", "HETATM") and line[17:20].strip() == "CR2":
            if line[76:78].strip().upper() == "H":
                continue
            atoms[line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    d = np.load(DENSITY)
    pts, q = d["pts_ang"], d["q"]
    mu = ((pts - pts.mean(0)) * q[:, None]).sum(0) * ANGSTROM_TO_BOHR
    return atoms, mu, pts.mean(0)


def random_rotations(n, rng):
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1),
    ], 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=4_000_000)
    ap.add_argument("--batch", type=int, default=20000)
    ap.add_argument("--alpha-lo", type=float, default=125.0)
    ap.add_argument("--alpha-hi", type=float, default=136.0)
    ap.add_argument("--sep-lo", type=float, default=23.0)
    ap.add_argument("--sep-hi", type=float, default=29.0)
    ap.add_argument("--min-contacts", type=int, default=40)
    ap.add_argument("--max-linker", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=20260812)
    ap.add_argument("--out", type=Path, default=Path("exciton_observables/register_search.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    heavy, cr2_a, ca = read_chain_a()
    ref_atoms, ref_mu, ref_origin = cr2_reference()
    shared = sorted(set(ref_atoms) & set(cr2_a))
    rot0, tr0 = kabsch(np.array([ref_atoms[n] for n in shared]),
                       np.array([cr2_a[n] for n in shared]))
    mu_a = rot0 @ ref_mu
    origin_a = rot0 @ ref_origin + tr0
    nterm, cterm = min(ca), max(ca)
    print(f"chain A: {len(heavy)} heavy atoms, residues {nterm}-{cterm}, "
          f"|mu| = {np.linalg.norm(mu_a):.3f} a.u.")

    # Occupancy grids for steric screening.
    lo = heavy.min(0) - (CONTACT_A + args.sep_hi + 60.0)
    hi = heavy.max(0) + (CONTACT_A + args.sep_hi + 60.0)
    shape = np.ceil((hi - lo) / VOXEL_A).astype(int) + 1
    print(f"occupancy grid {tuple(shape)} at {VOXEL_A} A "
          f"({np.prod(shape)/1e6:.0f} M voxels)")
    from scipy import ndimage
    occ = np.zeros(shape, dtype=bool)
    idx = np.floor((heavy - lo) / VOXEL_A).astype(int)
    occ[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    dist = ndimage.distance_transform_edt(~occ, sampling=VOXEL_A)
    clash_grid = dist < CLASH_A
    contact_grid = (dist >= CLASH_A) & (dist < CONTACT_A)
    del dist, occ

    mu_hat_a = mu_a / np.linalg.norm(mu_a)
    cos_lo = np.cos(np.radians(args.alpha_hi))
    cos_hi = np.cos(np.radians(args.alpha_lo))
    ca_cterm, ca_nterm = ca[cterm], ca[nterm]

    hits = []
    kept_geom = 0
    done = 0
    while done < args.trials:
        n = min(args.batch, args.trials - done)
        done += n
        rot = random_rotations(n, rng)
        mu_b = np.einsum("nij,j->ni", rot, mu_a)
        cosang = mu_b @ mu_hat_a / np.linalg.norm(mu_a)
        keep = (cosang >= cos_lo) & (cosang <= cos_hi)
        if not keep.any():
            continue
        rot, mu_b, cosang = rot[keep], mu_b[keep], cosang[keep]
        m = len(rot)

        direction = rng.normal(size=(m, 3))
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        radius = rng.uniform(args.sep_lo, args.sep_hi, m)
        r_vec = direction * radius[:, None]

        # Triple product (CD handedness) and coupling: both from this geometry.
        mu_hat_b = mu_b / np.linalg.norm(mu_b, axis=1)[:, None]
        triple = np.sum(r_vec * np.cross(np.broadcast_to(mu_hat_a, (m, 3)), mu_hat_b), axis=1)
        rb = r_vec * ANGSTROM_TO_BOHR
        rr = np.linalg.norm(rb, axis=1)
        rh = rb / rr[:, None]
        jdd = (mu_b @ mu_a) - 3.0 * (rh @ mu_a) * np.sum(mu_b * rh, axis=1)
        j_cm = jdd / (rr ** 3 * EPSILON) * HARTREE_TO_CM * TDC_OVER_PDA

        good = (triple < 0) & (j_cm > 20.0) & (j_cm < 50.0)
        if not good.any():
            continue
        rot, r_vec, cosang, j_cm, triple = (rot[good], r_vec[good], cosang[good],
                                            j_cm[good], triple[good])
        m = len(rot)
        kept_geom += m

        # Translation implied by putting chromophore B at origin_a + r_vec.
        trans = origin_a + r_vec - np.einsum("nij,j->ni", rot, origin_a)

        # Linker span first: cheap and rejects most.
        nterm_b = np.einsum("nij,j->ni", rot, ca_nterm) + trans
        span = np.linalg.norm(nterm_b - ca_cterm, axis=1)
        ok = span <= args.max_linker
        if not ok.any():
            continue
        rot, trans, cosang, j_cm, triple, r_vec, span = (
            rot[ok], trans[ok], cosang[ok], j_cm[ok], triple[ok], r_vec[ok], span[ok])

        for i in range(len(rot)):
            b = heavy @ rot[i].T + trans[i]
            gi = np.floor((b - lo) / VOXEL_A).astype(int)
            if (gi < 0).any() or (gi >= shape).any():
                continue
            cl = clash_grid[gi[:, 0], gi[:, 1], gi[:, 2]].sum()
            if cl:
                continue
            ct = int(contact_grid[gi[:, 0], gi[:, 1], gi[:, 2]].sum())
            if ct < args.min_contacts:
                continue
            hits.append({
                "alpha_deg": float(np.degrees(np.arccos(np.clip(cosang[i], -1, 1)))),
                "separation_A": float(np.linalg.norm(r_vec[i])),
                "J_cm": float(j_cm[i]),
                "triple_product": float(triple[i]),
                "contacts": ct,
                "linker_span_A": float(span[i]),
                "rotation": rot[i].tolist(),
                "translation": trans[i].tolist(),
            })
        if done % 400000 == 0:
            print(f"  {done:>9,} trials | {kept_geom:>7,} passed geometry | "
                  f"{len(hits):>5} full hits", flush=True)

    print(f"\n{done:,} trials, {kept_geom:,} passed alpha/sep/J/handedness, "
          f"{len(hits)} also passed sterics + linker")
    if hits:
        hits.sort(key=lambda h: -h["contacts"])
        print(f"\n{'alpha':>8}{'sep':>7}{'J':>8}{'contacts':>10}{'linker':>8}{'triple':>9}")
        for h in hits[:15]:
            print(f"{h['alpha_deg']:8.1f}{h['separation_A']:7.2f}{h['J_cm']:8.2f}"
                  f"{h['contacts']:10d}{h['linker_span_A']:8.1f}{h['triple_product']:9.2f}")
    args.out.write_text(json.dumps(
        {"settings": vars(args) | {"out": str(args.out)},
         "trials": done, "passed_geometry": kept_geom,
         "n_hits": len(hits), "hits": hits[:200]}, indent=2, default=str) + "\n")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

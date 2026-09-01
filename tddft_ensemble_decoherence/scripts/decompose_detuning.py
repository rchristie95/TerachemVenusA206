#!/usr/bin/env python3
"""Per-residue decomposition of the site-energy detuning Delta = E_A - E_B.

First-order electrochromic (Stark) shift. The excitation energy of a chromophore
in the field of an external point charge q_k at R_k shifts by

    dE = q_k * sum_a  dq_a / |r_a - R_k|

where dq_a are ATOM-CENTRED DIFFERENCE CHARGES (excited minus ground state),
obtained by Voronoi-partitioning the STEOM difference density onto the 41
physical QM atoms and rescaled so the atomic set reproduces the grid difference
dipole exactly. This is a monopole decomposition, not a point-dipole one, so it
stays valid for charges close to the chromophore.

Same machinery as the validated TD-TDX electrostatic analysis, which gave
15.6 +/- 4.0 cm^-1 of the 35.3 cm^-1 measured red shift.

Each MM atom contributes to BOTH sites, so its contribution to the gap is the
difference of its two shifts. Barrel-1 residue i pairs with barrel-2 residue
i+270; that pairing is what makes the result actionable, because a large paired
contribution means the two structurally equivalent positions are NOT
electrostatically equivalent, which is exactly what a mutation could fix.
"""
import json, sys
from pathlib import Path
import os

import numpy as np
from openmm.app import PDBFile

ROOT = Path(__file__).resolve().parent.parent
HARTREE_CM = 219474.6313632
# Which ensemble to decompose. The frame ids in the output are that
# window's own DCD frame numbers; do NOT correlate them against another
# window's frames by integer index -- different trajectories reuse the same
# integers for unrelated configurations.
WINDOW = os.environ.get('DETUNING_WINDOW', 'v5camp')
CACHE_TAG = 'v5' if WINDOW == 'v5camp' else 'v4'
TOPOLOGY = (Path('/home/robson/PetaChem/tc_tandem_v5_campaign/solvated_protonated.pdb')
            if WINDOW == 'v5camp' else
            Path('/home/robson/PetaChem/tc_tandem_v4_highcadence/solvated_protonated.pdb'))
BOHR_PER_ANG = 1.8897261246257702
BARREL_OFFSET = 270          # site A CR2 66 / Tyr 202 -> site B CR2 336 / Tyr 472
PROBE = Path('/home/robson/PetaChem/solvation_decoherence_test/'
             'steom_difference_probe_dipole_matched.npz')


def load_static():
    # Topology must match WINDOW: v4 and v5 are different builds with different
    # atom counts (74,857 vs 75,079), so a mismatched PDB makes the MM index
    # map silently wrong before the charge-count guard catches it.
    pdb = PDBFile(str(TOPOLOGY))
    atoms = list(pdb.topology.atoms())
    charges = np.load(ROOT / f'cache/embedding_charges_{CACHE_TAG}.npz')['charges_e']
    probe = np.load(PROBE, allow_pickle=True)
    return pdb, atoms, charges, probe['atom_delta_q_e'].astype(float)


def qm_indices(pdb, cr2_key, tyr_key, tyr_names):
    cr2 = tyr = None
    for r in pdb.topology.residues():
        k = (r.chain.id, r.id)
        if k == tuple(cr2_key) and r.name == 'CR2': cr2 = r
        if k == tuple(tyr_key) and r.name == 'TYR': tyr = r
    if cr2 is None or tyr is None:
        raise RuntimeError(f'residue lookup failed for {cr2_key} {tyr_key}')
    idx = [a.index for a in cr2.atoms()]
    by_name = {a.name: a for a in tyr.atoms()}
    idx += [by_name[n].index for n in tyr_names]
    return cr2, tyr, idx


def frame_site(frame_dir, site):
    """QM coords (41,3) in Angstrom and the MM charge/position table as written."""
    g = (frame_dir / f'site_{site}' / 'geometry.xyz').read_text().splitlines()
    xyz = np.array([[float(x) for x in l.split()[1:4]] for l in g[2:2 + 44]])[:41]
    raw = (frame_dir / f'site_{site}' / 'mm_charges.dat').read_text().splitlines()
    tab = np.array([[float(x) for x in l.split()] for l in raw[2:] if l.strip()])
    return xyz, tab[:, 0], tab[:, 1:4]


def main(frames, out_json):
    pdb, atoms, charges, dq = load_static()
    tyr_names = ['CG','CD1','CD2','CE1','CE2','CZ','OH','HD1','HD2','HE1','HE2','HH']
    prep0 = json.load(open(sorted(
        (ROOT / f'frames/{WINDOW}').glob('frame_*/preparation.json'))[0]))

    site_info = {}
    for s in 'AB':
        cr2_key = prep0['sites'][s]['cr2_key']; tyr_key = prep0['sites'][s]['stacked_tyr_key']
        cr2, tyr, qidx = qm_indices(pdb, cr2_key, tyr_key, tyr_names)
        links = [c['mm_partner_index_zero_based'] for c in prep0['sites'][s]['link_caps']]
        excluded = set(qidx) | set(links)
        mm_idx = np.array([a.index for a in atoms if a.index not in excluded])
        site_info[s] = dict(qidx=qidx, mm_idx=mm_idx, cr2=cr2, tyr=tyr)
        print(f'  site {s}: {len(qidx)} QM atoms, {len(links)} link partners, {mm_idx.size} MM charges')

    # residue label per topology atom
    res_of = np.empty(len(atoms), dtype=object)
    for r in pdb.topology.residues():
        for a in r.atoms():
            res_of[a.index] = (r.name, int(r.id) if str(r.id).isdigit() else -1)

    per_res = {}          # paired-position key -> list of per-frame gap contributions
    totals = []
    for fi in frames:
        fd = ROOT / f'frames/{WINDOW}/frame_{fi:04d}'
        if not (fd / 'site_A' / 'mm_charges.dat').exists():
            continue
        shift = {}
        for s in 'AB':
            xyz, q_written, pos = frame_site(fd, s)
            mm_idx = site_info[s]['mm_idx']
            if q_written.size != mm_idx.size:
                raise RuntimeError(f'frame {fi} site {s}: {q_written.size} charges vs {mm_idx.size} expected')
            # VALIDATION: the written charges must equal the cached charges in
            # topology order, EXCEPT at the partner chromophore. The cache stores
            # CR2 as 0.0 because load_or_build_charges deletes it by name;
            # --retain-partner-cr2-charges then restores its AMBER RESP values in
            # the written field. Permit exactly that difference and check it is
            # the physical one, rather than waving the mismatch through.
            diff = ~np.isclose(q_written, charges[mm_idx], atol=1e-6)
            if diff.any():
                names = {atoms[i].residue.name for i in mm_idx[diff]}
                restored = q_written[diff].sum()
                if (names != {'CR2'} or not np.allclose(charges[mm_idx][diff], 0.0)
                        or abs(restored + 1.0) > 1e-3):
                    raise RuntimeError(
                        f'frame {fi} site {s}: MM charge mapping mismatch at '
                        f'{int(diff.sum())} atoms in {sorted(names)}, '
                        f'restored charge {restored:+.4f} e (expected -1 on CR2)')
            d = np.linalg.norm(pos[:, None, :] - xyz[None, :, :], axis=2) * BOHR_PER_ANG
            shift[s] = q_written * (dq[None, :] / d).sum(axis=1) * HARTREE_CM
        gap_total = shift['A'].sum() - shift['B'].sum()
        totals.append((fi, gap_total))
        # attribute: every MM atom of site A, matched to the same atom in site B's list
        posA = {int(i): j for j, i in enumerate(site_info['A']['mm_idx'])}
        posB = {int(i): j for j, i in enumerate(site_info['B']['mm_idx'])}
        for topo_i in set(posA) & set(posB):
            g = shift['A'][posA[topo_i]] - shift['B'][posB[topo_i]]
            name, rid = res_of[topo_i]
            if name in ('HOH','WAT','TIP3','TIP3P'): key = ('solvent', -1)
            elif name in ('NA','CL','Na+','Cl-'):    key = ('ion', -1)
            else:
                base = rid - BARREL_OFFSET if rid > BARREL_OFFSET else rid
                key = (name, base)
            per_res.setdefault(key, []).append(g)
        print(f'  frame {fi:5d}: gap from electrostatics = {gap_total:+9.1f} cm^-1', flush=True)

    summary = []
    nf = len(totals)
    for key, vals in per_res.items():
        v = np.array(vals)
        summary.append(dict(residue=key[0], position=key[1],
                            mean_cm=float(v.sum() / nf), rms_cm=float(np.sqrt((v ** 2).sum() / nf)),
                            n=int(v.size)))
    summary.sort(key=lambda d: -abs(d['rms_cm']))
    json.dump(dict(frames=[int(f) for f, _ in totals],
                   gap_totals={str(f): float(g) for f, g in totals},
                   contributions=summary), open(out_json, 'w'), indent=2)
    print(f'\nwrote {out_json}: {nf} frames, {len(summary)} residue groups')
    return summary, totals


if __name__ == '__main__':
    frames = [int(x) for x in sys.argv[1:-1]] if len(sys.argv) > 2 else [0, 26, 48]
    main(frames, sys.argv[-1])

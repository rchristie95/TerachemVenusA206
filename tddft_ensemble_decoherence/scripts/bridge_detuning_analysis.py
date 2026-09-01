#!/usr/bin/env python3
"""What sets the site-energy detuning in a fluorescent-protein dimer?

Tests two hypotheses against the full-system-embedded QM/MM frames of the v5
campaign (500 frames over 1 ns; --window v4full reproduces the earlier 92-frame
analysis):
  (a) the electrostatic environment differs between the two sites
      -> first-order electrochromic model, results_detuning_decomposition.json
  (b) the chromophores' own internal geometry differs
      -> bond/dihedral descriptors computed here

Writes results_bridge_detuning.json and Fig_Bridge_Detuning.pdf.
"""
import argparse
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
WINDOW = 'v5camp'
CA2, CB2, CG2, CD1, N2, CZ, OH = 4, 13, 14, 15, 3, 17, 22   # 0-based, geometry.xyz
TWO_J = 54.1     # 2|J| from the 500-frame per-frame-density ensemble


def geom(fd, s):
    L = (fd / f'site_{s}' / 'geometry.xyz').read_text().splitlines()
    return np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + 44]])


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def descriptors(X):
    d = lambda i, j: float(np.linalg.norm(X[i] - X[j]))
    return dict(bridge=d(CA2, CB2) + d(CB2, CG2),
                BLA=d(CA2, CB2) - d(CB2, CG2),
                CO=d(CZ, OH),
                tau=dihedral(X[N2], X[CA2], X[CB2], X[CG2]),
                phi=dihedral(X[CA2], X[CB2], X[CG2], X[CD1]))


def main():
    rows = []
    for p in sorted(glob.glob(str(ROOT / f'frames/{WINDOW}/frame_*/frame_summary.json'))):
        d = json.load(open(p))
        A, B = d.get('site_A') or {}, d.get('site_B') or {}
        if A.get('status') != 'selected' or B.get('status') != 'selected':
            continue
        fd = Path(p).parent
        dA, dB = descriptors(geom(fd, 'A')), descriptors(geom(fd, 'B'))
        r = dict(frame=int(d['frame_index']), gap=A['energy_cm-1'] - B['energy_cm-1'],
                 bA=dA['bridge'], bB=dB['bridge'])
        for k in dA:
            a, b = dA[k], dB[k]
            if k in ('tau', 'phi'):
                a, b = abs((a + 90) % 180 - 90), abs((b + 90) % 180 - 90)
            r['d_' + k] = a - b
        rows.append(r)

    g = np.array([r['gap'] for r in rows])
    db = np.array([r['d_bridge'] for r in rows])
    bA = np.array([r['bA'] for r in rows]); bB = np.array([r['bB'] for r in rows])
    keys = [k for k in rows[0] if k.startswith('d_')]
    corr = {k: float(np.corrcoef([r[k] for r in rows], g)[0, 1]) for k in keys}

    # The electrostatic decomposition must come from the SAME window. Frame ids
    # are per-trajectory DCD indices, so intersecting two windows by integer
    # index matches unrelated configurations: v4full's 92 frames (ids 0-1999)
    # and v5camp's 500 (ids 900-2896) share 8 integers that correspond to
    # nothing physical, and correlating on them returned r = -0.465 where the
    # matched analysis gives +0.10.
    dec_path = ROOT / (f'results_detuning_decomposition_v5.json' if WINDOW == 'v5camp'
                       else 'results_detuning_decomposition.json')
    dec = json.load(open(dec_path))
    qm = {r['frame']: r['gap'] for r in rows}
    fr = [int(k) for k in dec['gap_totals'] if int(k) in qm]
    if len(fr) < 0.5 * len(rows):
        raise RuntimeError(
            f'{dec_path.name} overlaps only {len(fr)} of {len(rows)} {WINDOW} frames; '
            f'it was computed for a different window and must not be correlated '
            f'against these gaps by frame number')
    es = np.array([dec['gap_totals'][str(i)] for i in fr])
    r_es = float(np.corrcoef(es, [qm[i] for i in fr])[0, 1])
    print(f'  electrostatic model: {len(fr)} frames matched from {dec_path.name}')

    X = np.column_stack([[r[k] for r in rows] for k in keys] + [np.ones(len(rows))])
    beta, *_ = np.linalg.lstsq(X, g, rcond=None)
    R2 = float(1 - ((g - X @ beta) ** 2).sum() / ((g - g.mean()) ** 2).sum())

    slope, icpt = np.polyfit(db, g, 1)
    need_mA = TWO_J / abs(slope) * 1000
    out = dict(
        n_frames=len(rows), gap_sd_cm=float(g.std(ddof=1)),
        bridge_A_mean_A=float(bA.mean()), bridge_A_sd_A=float(bA.std(ddof=1)),
        bridge_B_mean_A=float(bB.mean()), bridge_B_sd_A=float(bB.std(ddof=1)),
        bridge_diff_mean_A=float(db.mean()), bridge_diff_sd_A=float(db.std(ddof=1)),
        sensitivity_cm_per_A=float(slope),
        gap_sd_predicted_from_bridge_cm=float(abs(slope) * db.std(ddof=1)),
        correlations=corr, r_electrostatic=r_es, R2_all_geometry=R2,
        two_J_cm=TWO_J, required_bridge_match_mA=float(need_mA),
        reduction_factor=float(db.std(ddof=1) * 1000 / need_mA),
        single_site_thermal_sd_mA=float(bA.std(ddof=1) * 1000),
    )
    json.dump(out, open(ROOT / 'results_bridge_detuning.json', 'w'), indent=2)

    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
    ax[0].scatter(db * 1000, g, s=20, alpha=0.75, edgecolor='none', color='#1f77b4')
    xs = np.linspace(db.min(), db.max(), 50)
    ax[0].plot(xs * 1000, slope * xs + icpt, 'k-', lw=1.3)
    ax[0].axhspan(-TWO_J, TWO_J, color='0.80', zorder=0)
    ax[0].annotate(r'$|\Delta|<2|J|$: delocalisation window',
                   xy=(db.max() * 1000 * 0.55, 0), xytext=(db.max() * 1000 * 0.20, 1450),
                   fontsize=8, color='0.25',
                   arrowprops=dict(arrowstyle='->', color='0.35', lw=0.9))
    ax[0].set_xlabel(r'bridge-length mismatch $b_A-b_B$ (m$\AA$)')
    ax[0].set_ylabel(r'site-energy gap $E_A-E_B$ (cm$^{-1}$)')
    ax[0].set_title(f"$r={corr['d_bridge']:+.3f}$, {corr['d_bridge']**2*100:.0f}% of variance", fontsize=9)
    labels = ['bridge\nlength', 'all geometry\njointly', 'electrostatics\n(Stark)']
    vals = [corr['d_bridge'] ** 2 * 100, R2 * 100, r_es ** 2 * 100]
    ax[1].bar(labels, vals, color=['#1f77b4', '#4c9bd4', '#d62728'])
    for i, v in enumerate(vals):
        ax[1].text(i, v + 1.5, f'{v:.1f}%', ha='center', fontsize=9)
    ax[1].set_ylabel('variance of $E_A-E_B$ explained (%)')
    ax[1].set_ylim(0, 70)
    fig.tight_layout()
    fig.savefig(ROOT.parent / 'manuscript' / 'Fig_Bridge_Detuning.pdf')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()

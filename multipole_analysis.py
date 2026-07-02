#!/usr/bin/env python3
r"""
multipole_analysis.py  --  Multipole decomposition of the excitonic coupling.

Reviewer item 3 (R4 explicit, R1 conceptual): decompose the factor-of-~5.6
enhancement of the full Transition Density Coupling (J_TDC = 74.38 cm^-1) over
the point-dipole approximation (J_PDA = 13.31 cm^-1) into multipole orders --
dipole-dipole, dipole-quadrupole, quadrupole-quadrupole (and higher) -- to show
which terms drive the near-field enhancement at the 27.6 A centroid separation,
and where the PDA breaks down.

Method (primitive Cartesian multipole expansion of the Coulomb interaction):
  The exact coupling is  J = sum_ij q_i q_j / |R + a_i - b_j|  (atomic units),
  with a_i / b_j the transition-density grid points expressed about each site's
  centroid and R the inter-centroid vector. Taylor-expanding 1/|R+s| in
  s = a_i - b_j and grouping by total order n gives

      J = sum_n sum_{p+q=n} [(-1)^q / (p! q!)] * M^A_(p) (x) M^B_(q) : T^(n)(R)

  where  M_(p)_{a..} = sum_i q_i a_{i,a}...   are primitive moments
  (p=1 dipole, p=2 second moment / quadrupole, p=3 octupole) and
  T^(n) = grad^n (1/R) are the interaction tensors. The n=2,p=q=1 term is
  algebraically identical to the PDA (verified in the self-test), so the higher
  orders are exactly the corrections the PDA misses.

Reuses coupling_core.py (read_dx, get_super_matrices_with_pymol,
apply_pymol_matrix, calculate_coupling, ...) for the real-data path, identical
to qmmm_tddft_pipeline.py Stage 3.

Run the built-in self-test (no PyMOL / GPU needed):
    python multipole_analysis.py --self-test

Real data:
    python multipole_analysis.py --workdir tc_tddft_old_current \
        --monomer tc_simple_old/classical_relaxed.pdb --dimer venus_dimer.pdb

Outputs (in --out): multipole_analysis.csv, Fig_Multipole_Decomposition.pdf
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

from coupling_core import (
    HARTREE_TO_CM,
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    read_dx,
    calculate_coupling,
    get_super_matrices_with_pymol,
    apply_pymol_matrix,
    transition_dipole_au,
    oscillator_to_dipole_au,
    autodetect_workdir_and_candidates,
    select_target_state_and_density,
    print_excited_state_table,
)


# --------------------------------------------------------------------------- #
# Interaction tensors  T^(n) = grad^n (1/R)   (R in Bohr; result in 1/Bohr^(n+1))
# --------------------------------------------------------------------------- #
def interaction_tensors(R_vec):
    """Return T2 (3,3), T3 (3,3,3), T4 (3,3,3,3) for f = 1/R."""
    R = float(np.linalg.norm(R_vec))
    n = np.asarray(R_vec, float)
    d = np.eye(3)

    # T2_ab = (3 n_a n_b - R^2 d_ab) / R^5
    T2 = (3 * np.outer(n, n) - R**2 * d) / R**5

    # T3_abc = -(15 n_a n_b n_c - 3 R^2 (d_ab n_c + d_ac n_b + d_bc n_a)) / R^7
    nnn = np.einsum("a,b,c->abc", n, n, n)
    dn = (np.einsum("ab,c->abc", d, n)
          + np.einsum("ac,b->abc", d, n)
          + np.einsum("bc,a->abc", d, n))
    T3 = -(15 * nnn - 3 * R**2 * dn) / R**7

    # T4_abcd = [105 n_a n_b n_c n_d
    #            - 15 R^2 (d_ab n_c n_d + d_ac n_b n_d + d_ad n_b n_c
    #                      + d_bc n_a n_d + d_bd n_a n_c + d_cd n_a n_b)
    #            + 3 R^4 (d_ab d_cd + d_ac d_bd + d_ad d_bc)] / R^9
    nnnn = np.einsum("a,b,c,e->abce", n, n, n, n)
    dnn = (np.einsum("ab,c,e->abce", d, n, n)
           + np.einsum("ac,b,e->abce", d, n, n)
           + np.einsum("ae,b,c->abce", d, n, n)
           + np.einsum("bc,a,e->abce", d, n, n)
           + np.einsum("be,a,c->abce", d, n, n)
           + np.einsum("ce,a,b->abce", d, n, n))
    dd = (np.einsum("ab,ce->abce", d, d)
          + np.einsum("ac,be->abce", d, d)
          + np.einsum("ae,bc->abce", d, d))
    T4 = (105 * nnnn - 15 * R**2 * dnn + 3 * R**4 * dd) / R**9
    return T2, T3, T4


def primitive_moments(pts_bohr, q, center_bohr):
    """Primitive Cartesian transition moments about `center` (atomic units)."""
    a = pts_bohr - center_bohr
    M0 = float(np.sum(q))
    M1 = np.einsum("i,ia->a", q, a)
    M2 = np.einsum("i,ia,ib->ab", q, a, a)
    M3 = np.einsum("i,ia,ib,ic->abc", q, a, a, a)
    return M0, M1, M2, M3


def multipole_terms(ptsA_ang, qA, ptsB_ang, qB):
    """
    Decompose the coupling (Hartree, in vacuum) by multipole order.

    Returns dict with named contributions:
      dip-dip (n=2), dip-quad (n=3), quad-quad + dip-oct (n=4 split).
    """
    a = ptsA_ang * ANGSTROM_TO_BOHR
    b = ptsB_ang * ANGSTROM_TO_BOHR
    cA = a.mean(axis=0)
    cB = b.mean(axis=0)
    # The pairwise separation vector is (cA + a_i) - (cB + b_j) = D0 + (a_i - b_j)
    # with D0 = cA - cB, so the interaction tensors must be evaluated at D0.
    # (Even-order tensors are sign-insensitive, so dip-dip still equals the PDA,
    # but odd-order T3 flips sign -- getting D0 wrong breaks dip-quad convergence.)
    R_vec = cA - cB

    _, A1, A2, A3 = primitive_moments(a, qA, cA)
    _, B1, B2, B3 = primitive_moments(b, qB, cB)
    T2, T3, T4 = interaction_tensors(R_vec)

    # coefficient (-1)^q / (p! q!)
    def coef(p, q):
        return ((-1) ** q) / (math.factorial(p) * math.factorial(q))

    # n=2 : dipole-dipole (p=q=1)  -> equals the PDA
    dd = coef(1, 1) * np.einsum("a,b,ab->", A1, B1, T2)

    # n=3 : dipole-quadrupole  (p=1,q=2) + (p=2,q=1)
    dq = (coef(1, 2) * np.einsum("a,bc,abc->", A1, B2, T3)
          + coef(2, 1) * np.einsum("ab,c,abc->", A2, B1, T3))

    # n=4 : quadrupole-quadrupole (p=q=2)
    qq = coef(2, 2) * np.einsum("ab,cd,abcd->", A2, B2, T4)
    # n=4 : dipole-octupole (p=1,q=3)+(p=3,q=1)
    do = (coef(1, 3) * np.einsum("a,bcd,abcd->", A1, B3, T4)
          + coef(3, 1) * np.einsum("abc,d,abcd->", A3, B1, T4))

    return {
        "dip-dip (PDA)": float(dd),
        "dip-quad": float(dq),
        "quad-quad": float(qq),
        "dip-oct": float(do),
    }


def exact_coupling_bruteforce(ptsA_ang, qA, ptsB_ang, qB):
    """Exact O(NM) Coulomb double sum (Hartree, vacuum). For validation / small sets."""
    a = ptsA_ang * ANGSTROM_TO_BOHR
    b = ptsB_ang * ANGSTROM_TO_BOHR
    diff = a[:, None, :] - b[None, :, :]
    r = np.sqrt(np.sum(diff**2, axis=2))
    return float(np.sum((qA[:, None] * qB[None, :]) / r))


# --------------------------------------------------------------------------- #
# Self-test: reconstruct an exact small-cluster coupling from the multipoles.
# --------------------------------------------------------------------------- #
def self_test():
    rng = np.random.default_rng(7)
    # Two compact (~2 A) charge-neutral clusters separated by ~28 A along an
    # arbitrary axis -- mimics two chromophore transition densities.
    def neutral_cluster(center):
        pts = center + rng.normal(0, 1.0, size=(8, 3))
        q = rng.normal(0, 1.0, size=8)
        q -= q.mean()  # enforce charge neutrality (transition density)
        return pts, q

    axis = np.array([1.0, 0.6, 0.3]); axis /= np.linalg.norm(axis)
    ptsA, qA = neutral_cluster(np.zeros(3))
    ptsB, qB = neutral_cluster(28.0 * axis)

    exact = exact_coupling_bruteforce(ptsA, qA, ptsB, qB)
    terms = multipole_terms(ptsA, qA, ptsB, qB)
    cumulative = sum(terms.values())

    print("Self-test: multipole reconstruction of an exact two-cluster coupling")
    print(f"  exact (brute-force)        : {exact*HARTREE_TO_CM: .4f}  (scaled)")
    running = 0.0
    for name, val in terms.items():
        running += val
        print(f"  + {name:14s} = {val*HARTREE_TO_CM: .4f}   cumulative {running*HARTREE_TO_CM: .4f}")
    rel_err = abs(cumulative - exact) / (abs(exact) + 1e-30)
    print(f"  cumulative / exact         : {cumulative/exact:.4f}  (rel. error {rel_err:.3%})")
    # Dipole-dipole term must equal the PDA computed independently.
    pda = pda_reference(ptsA, qA, ptsB, qB)
    print(f"  dip-dip vs independent PDA : {terms['dip-dip (PDA)']*HARTREE_TO_CM:.4f} "
          f"vs {pda*HARTREE_TO_CM:.4f}")
    ok = rel_err < 0.05 and abs(terms["dip-dip (PDA)"] - pda) < 1e-6 * (abs(pda) + 1)
    print("  RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def pda_reference(ptsA_ang, qA, ptsB_ang, qB):
    """Independent point-dipole estimate (Hartree, vacuum), matching Stage 3."""
    a = ptsA_ang * ANGSTROM_TO_BOHR
    b = ptsB_ang * ANGSTROM_TO_BOHR
    cA, cB = a.mean(axis=0), b.mean(axis=0)
    muA = transition_dipole_au(ptsA_ang, qA, origin_angstrom=cA * BOHR_TO_ANGSTROM)
    muB = transition_dipole_au(ptsB_ang, qB, origin_angstrom=cB * BOHR_TO_ANGSTROM)
    R = cB - cA
    Rn = np.linalg.norm(R)
    Rhat = R / Rn
    return float((np.dot(muA, muB) - 3 * np.dot(muA, Rhat) * np.dot(muB, Rhat)) / Rn**3)


# --------------------------------------------------------------------------- #
# Plot / IO
# --------------------------------------------------------------------------- #
def write_csv(terms, cumulative, tdc, pda, enhancement, path):
    with open(path, "w") as f:
        f.write("term,contribution_cm,cumulative_cm\n")
        run = 0.0
        for name, val in terms.items():
            run += val
            f.write(f"{name},{val:.4f},{run:.4f}\n")
        f.write(f"multipole_sum,,{cumulative:.4f}\n")
        f.write(f"full_TDC,,{tdc:.4f}\n")
        f.write(f"PDA_reference,,{pda:.4f}\n")
        f.write(f"TDC_over_PDA,,{enhancement:.4f}\n")


def plot_decomposition(terms, cumulative, tdc, pda, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(terms.keys())
    vals = [terms[n] for n in names]
    run = np.cumsum(vals)

    fig, ax = plt.subplots(figsize=(7, 4.4))
    x = np.arange(len(names))
    ax.bar(x, vals, color="#4C72B0", alpha=0.85, label="order contribution")
    ax.plot(x, run, "o-", color="#C44E52", lw=2, label="cumulative multipole")
    ax.axhline(tdc, color="k", ls="--", lw=1.6, label=f"full TDC = {tdc:.1f}")
    ax.axhline(pda, color="0.5", ls=":", lw=1.6, label=f"PDA = {pda:.1f}")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel(r"coupling contribution (cm$^{-1}$)")
    enh = tdc / pda if pda else float("nan")
    ax.set_title(fr"Multipole decomposition (TDC/PDA $\approx$ {enh:.1f}$\times$)")
    ax.legend(fontsize=8, frameon=True)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Real-data driver (reuses the Stage 3 flow)
# --------------------------------------------------------------------------- #
def run_real(args):
    workdir, candidates, tried = autodetect_workdir_and_candidates(args.workdir)
    if not candidates:
        print(f"[!] No excited-state candidates. Tried: {tried}")
        sys.exit(1)
    print_excited_state_table(candidates)
    state, density_file, _, mode = select_target_state_and_density(
        candidates, workdir, args.density_mode, requested_root=args.root)
    mu_target = oscillator_to_dipole_au(state["ev"], state["osc"])
    print(f"    - Root {state['root']} ({mode}): f={state['osc']:.4f}, |mu|={mu_target:.4f} a.u.")

    pts_opt, q_opt = read_dx(density_file, threshold=args.thresh, stride=args.grid_stride)
    if pts_opt.size == 0:
        print("[!] Empty transition density.")
        sys.exit(1)

    # Renormalise to the oscillator-strength dipole (Stage 3 convention).
    if np.isfinite(mu_target):
        local_origin = np.mean(pts_opt, axis=0)
        dip_mag = np.linalg.norm(np.dot(q_opt, pts_opt - local_origin)) / BOHR_TO_ANGSTROM
        if dip_mag > 1e-6:
            q_opt = q_opt * (mu_target / dip_mag)

    matrix_A, matrix_B, aln_A, aln_B, err = get_super_matrices_with_pymol(args.monomer, args.dimer)
    if err:
        print(f"[!] PyMOL error: {err}")
        sys.exit(1)
    pts_A = apply_pymol_matrix(pts_opt, matrix_A)
    pts_B = apply_pymol_matrix(pts_opt, matrix_B)

    # Full TDC (vacuum) and PDA, then per-order decomposition.
    tdc_vac = calculate_coupling(
        pts_A, q_opt, pts_B, q_opt, backend=args.backend,
        gpu_chunk=args.gpu_chunk, opencl_chunk=args.opencl_chunk,
        opencl_platform=args.opencl_platform, opencl_device=args.opencl_device)
    tdc = tdc_vac * HARTREE_TO_CM / args.epsilon
    pda = pda_reference(pts_A, q_opt, pts_B, q_opt) * HARTREE_TO_CM / args.epsilon
    terms_ha = multipole_terms(pts_A, q_opt, pts_B, q_opt)
    terms = {k: v * HARTREE_TO_CM / args.epsilon for k, v in terms_ha.items()}
    return terms, tdc, pda


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--self-test", action="store_true", help="Run the numerical validation and exit.")
    p.add_argument("--workdir", type=Path, default=Path("tc_tddft_old_current"))
    p.add_argument("--monomer", type=Path, default=Path("tc_simple_old/classical_relaxed.pdb"))
    p.add_argument("--dimer", type=Path, default=Path("venus_dimer.pdb"))
    p.add_argument("--root", type=int, default=None)
    p.add_argument("--density-mode", choices=["signed", "auto", "abs"], default="auto")
    p.add_argument("--thresh", type=float, default=1e-7)
    p.add_argument("--grid-stride", type=int, default=1)
    p.add_argument("--epsilon", type=float, default=1.77)
    p.add_argument("--backend", default="auto", choices=["auto", "gpu", "opencl"])
    p.add_argument("--gpu-chunk", type=int, default=10000)
    p.add_argument("--opencl-chunk", type=int, default=20000)
    p.add_argument("--opencl-platform", type=int, default=None)
    p.add_argument("--opencl-device", type=int, default=None)
    p.add_argument("--out", type=Path, default=Path("multipole_out"))
    args = p.parse_args(argv)

    if args.self_test:
        sys.exit(self_test())

    args.out.mkdir(parents=True, exist_ok=True)
    terms, tdc, pda = run_real(args)
    cumulative = sum(terms.values())
    enhancement = tdc / pda if pda else float("nan")

    print("\n--- MULTIPOLE DECOMPOSITION (cm^-1) ---")
    run = 0.0
    for name, val in terms.items():
        run += val
        print(f"  {name:14s} = {val: .3f}   cumulative {run: .3f}")
    print(f"  multipole sum   = {cumulative: .3f}")
    print(f"  full TDC        = {tdc: .3f}")
    print(f"  PDA reference   = {pda: .3f}")
    print(f"  TDC / PDA       = {enhancement: .2f}x")

    csv_path = args.out / "multipole_analysis.csv"
    fig_path = args.out / "Fig_Multipole_Decomposition.pdf"
    write_csv(terms, cumulative, tdc, pda, enhancement, csv_path)
    try:
        plot_decomposition(terms, cumulative, tdc, pda, fig_path)
    except Exception as exc:
        print(f"[!] Plotting failed: {exc}")
        fig_path = None
    print(f"\n  csv    : {csv_path}")
    if fig_path:
        print(f"  figure : {fig_path}")


if __name__ == "__main__":
    main()

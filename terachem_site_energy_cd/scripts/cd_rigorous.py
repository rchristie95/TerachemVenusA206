#!/usr/bin/env python3
"""Absolute CD of the Venus tandem dimer with no fitted or borrowed lineshape.

Supersedes scripts/cd_from_simulation.py, which had four aggregate knobs doing
real work.  Every one is removed here:

  REMOVED  the multiplicative detuning rescale (k_delta = 0.295).  It multiplied
           each frame's instantaneous Delta by a constant so the ensemble spread
           matched an assumed slow component -- a construction, and the single
           largest lever in the figure.  Replaced by an exact deconvolution:
           the v5 gap, the v5 mean site energy and the within-window residuals
           are each Gaussian (Shapiro p = 0.64, 0.36, 0.11), so
           sigma_slow^2 = sigma_total^2 - sigma_fast^2 is exact, and the static
           disorder is integrated by Gauss-Hermite quadrature rather than by
           deforming individual frames.  Legitimate because the detuning is
           measured to be uncorrelated with every geometric amplitude in the
           model: |r| <= 0.12 against J, T, cos(alpha) and separation over 500
           frames, the largest being 2.6 sigma and worth 1.4% of variance.

  REMOVED  the fitted homogeneous HWHM of 340 cm^-1.  It came from a second-order
           cumulant on the classical gap ACF, which assigns to homogeneous
           dephasing the same nuclear motion that the vibronic progression below
           already carries -- a double count.  The manuscript's own position is
           that no classical mode qualifies as a bath: everything below kT has a
           period > 160 fs and is quasi-static, everything above 834 cm^-1 is
           quantum-frozen and makes Franck-Condon structure.  There is therefore
           NO homogeneous width in this model.  The radiative contribution
           (~3 ns) is 0.002 cm^-1 and is neglected.

  REMOVED  the single effective vibronic mode (S = 0.253, omega = 1369.7 cm^-1),
           imported from a separate electronic-structure calculation.  The mode
           set is now read off the measured per-site energy spectral density of
           the three 1 fs windows: the variance in each resolved band gives
           lambda_i = sigma_i^2 / 2kT and hence S_i = lambda_i / omega_i.  Bands
           at or below kT are not modes at all on this timescale and are moved
           into the static disorder, which is where the resolution limit puts
           them honestly.

  REMOVED  the equal-|mu| homodimer assumption.  The exciton transition dipole is
           now the vector sum mu_k = sum_n c^k_{A,n} F_n mu_A + (A -> B), so the
           interference term is exact and the 6% difference between |mu_A| and
           |mu_B| is carried rather than averaged away.

  REPLACED sigma_site_slow = 159 cm^-1, which differenced a 1 ns run against a
           592 fs run from a DIFFERENT trajectory.  Both halves of the split now
           share one lineage (v5camp 1 ns, v6t2 3 x 100 fs at 1 fs).

What remains taken from outside the simulation: the band origin, in the
`anchored` case only.  The `abinitio` case takes it from the eight matched
STEOM/TDDFT site pairs and is reported alongside.

Read-only apart from results/cd_rigorous.{json,npz}.
"""
from __future__ import annotations

import csv
import itertools
import math
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))
from absolute_cd import ELLIPTICITY_PER_DELTA_EPSILON, DEBYE_PER_AU  # noqa: E402

TDDFT = REPO / "tddft_ensemble_decoherence"
ENSEMBLE = ROOT / "results" / "ensembles" / "ens_v5camp.npz"
OUT = ROOT / "results"

KT_CM = 208.5                      # 300 K
C_CM_S = 29979245800.0
STEOM_SHIFT_EV = -0.9745050587499999
CM_PER_EV = 8065.543937

EXP_CENTER_CM = 19322.0
EXP_HALF_SPLIT = 130.79
DIFF_SIGN = -1.0
EXP_BANDS = tuple((DIFF_SIGN * a * 1e5, c, w) for a, c, w in (
    (+3.59, EXP_CENTER_CM - EXP_HALF_SPLIT, 293.51),
    (-2.27, EXP_CENTER_CM + EXP_HALF_SPLIT, 439.25)))
EPS_TD_EXP = 184400.0

R_TO_DELTA_EPS = 1.0 / 2.296e-39   # R  = 2.296e-39 int (Delta_eps / nubar) dnubar
D_TO_EPS = 1.0 / 9.184e-39         # D  = 9.184e-39 int (eps        / nubar) dnubar
DEBYE_ESU_CM = 1.0e-18

GRID = np.arange(15000.0, 25000.0, 2.0)
N_QUANTA = 3                       # 0,1,2 per mode; S_i <= 0.13 so P(n=2) < 1%
N_HERMITE = 24                     # Gauss-Hermite nodes for the static detuning
BAND_EDGES = (400.0, 900.0, 1500.0, 2100.0, 3000.0)
S_MIN = 5e-3                       # modes below this carry <0.3% of the variance


# --------------------------------------------------------------------------
# 1.  measured disorder decomposition, from the frame data of both campaigns
# --------------------------------------------------------------------------
def site_series():
    """Per-site site-energy trajectories from the three 1 fs windows."""
    out = []
    for w in "abc":
        rows = [r for r in csv.DictReader(
            (TDDFT / f"manifests/v6t2{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
        rows.sort(key=lambda r: int(r["frame_index"]))
        for site in ("A", "B"):
            out.append(np.array([float(r[f"site_{site}_energy_cm-1"])
                                 for r in rows]))
    return out


def disorder_and_modes(series, sigma_site_total):
    """Split the measured fluctuation into static disorder and FC modes.

    Returns (sigma_site_static, [(omega_i, S_i)], diagnostics).
    """
    n = series[0].size
    win = np.hanning(n)
    freqs = np.fft.rfftfreq(n, 1e-15) / C_CM_S
    power = np.zeros(freqs.size)
    for s in series:
        power += np.abs(np.fft.rfft((s - s.mean()) * win)) ** 2
    power /= power.sum()

    var_fast = float(np.mean([s.var(ddof=1) for s in series]))
    var_slow = float(sigma_site_total ** 2 - var_fast)

    modes, quasi_static_var, diag = [], 0.0, []
    edges = (0.0,) + BAND_EDGES + (float(freqs[-1]) + 1.0,)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (freqs >= lo) & (freqs < hi)
        if not m.any():
            continue
        frac = float(power[m].sum())
        if frac < 1e-3:
            continue
        omega = float((freqs[m] * power[m]).sum() / power[m].sum())
        var_i = frac * var_fast
        lam_i = var_i / (2.0 * KT_CM)
        if omega <= 2.0 * KT_CM:          # not a mode on the absorption timescale
            quasi_static_var += var_i
            diag.append(dict(band=[lo, hi], omega_cm=omega, var_fraction=frac,
                             treatment="quasi-static"))
        elif lam_i / omega >= S_MIN:
            modes.append((omega, lam_i / omega))
            diag.append(dict(band=[lo, hi], omega_cm=omega, var_fraction=frac,
                             lambda_cm=lam_i, huang_rhys=lam_i / omega,
                             treatment="Franck-Condon"))
        else:                              # negligible displacement; fold into static
            quasi_static_var += var_i
            diag.append(dict(band=[lo, hi], omega_cm=omega, var_fraction=frac,
                             huang_rhys=lam_i / omega, treatment="negligible"))
    sigma_static = float(np.sqrt(var_slow + quasi_static_var))
    return sigma_static, modes, dict(
        sigma_site_total_cm=float(sigma_site_total),
        sigma_site_fast_cm=float(np.sqrt(var_fast)),
        sigma_site_slow_cm=float(np.sqrt(max(var_slow, 0.0))),
        sigma_site_quasi_static_cm=float(np.sqrt(quasi_static_var)),
        sigma_site_static_cm=sigma_static,
        lambda_site_fast_cm=float(var_fast / (2.0 * KT_CM)),
        spectral_peak_cm=float(freqs[1:][power[1:].argmax()]),
        resolution_cm=float(freqs[1]),
        bands=diag)


# --------------------------------------------------------------------------
# 2.  multi-mode one-particle Holstein dimer
# --------------------------------------------------------------------------
def vibrational_basis(modes, n_quanta):
    """Occupation tuples, their energies, and multi-mode FC overlaps <n|0>."""
    occ = list(itertools.product(range(n_quanta), repeat=len(modes)))
    energy = np.array([sum(n * w for n, (w, _) in zip(o, modes)) for o in occ])
    fc = np.ones(len(occ))
    for k, (_, s) in enumerate(modes):
        amp = np.exp(-s / 2.0)
        for idx, o in enumerate(occ):
            n = o[k]
            fc[idx] *= amp * (-np.sqrt(s)) ** n / np.sqrt(float(math.factorial(n)))
    return energy, fc


def frame_sticks(e_a, e_b, j_cm, mu_a, mu_b, triple, vib_e, fc):
    """Exact exciton states of one configuration: (nu_k, R_k, D_k)."""
    nv = vib_e.size
    dim = 2 * nv
    h = np.zeros((dim, dim))
    h[:nv, :nv] = np.diag(e_a + vib_e)
    h[nv:, nv:] = np.diag(e_b + vib_e)
    block = j_cm * np.outer(fc, fc)
    h[:nv, nv:] = block
    h[nv:, :nv] = block.T
    nu, vec = np.linalg.eigh(h)

    a_k = fc @ vec[:nv, :]                       # site-A vibronic amplitude
    b_k = fc @ vec[nv:, :]
    mu_k = np.outer(a_k, mu_a) + np.outer(b_k, mu_b)      # exact vector sum
    d_k = np.einsum("ki,ki->k", mu_k, mu_k) * DEBYE_ESU_CM ** 2
    chirality = a_k * b_k          # == einsum("n,nk,m,mk->k", fc, ca, fc, cb)
    r_k = -np.pi * nu * chirality * triple
    return nu, r_k, d_k


def accumulate(nu, weight, grid, out):
    """Deposit sticks onto the grid by linear interpolation (area conserving)."""
    dx = grid[1] - grid[0]
    pos = (nu - grid[0]) / dx
    lo = np.floor(pos).astype(int)
    frac = pos - lo
    ok = (lo >= 0) & (lo < grid.size - 1)
    np.add.at(out, lo[ok], weight[ok] * (1.0 - frac[ok]) / dx)
    np.add.at(out, lo[ok] + 1, weight[ok] * frac[ok] / dx)


def gaussian_convolve(spectrum, grid, sigma):
    dx = grid[1] - grid[0]
    half = int(np.ceil(5.0 * sigma / dx))
    x = np.arange(-half, half + 1) * dx
    k = np.exp(-0.5 * (x / sigma) ** 2)
    k /= k.sum()
    return np.convolve(spectrum, k, mode="same")


# --------------------------------------------------------------------------
def main() -> None:
    d = np.load(ENSEMBLE)
    e_a, e_b, j = d["e_a"], d["e_b"], d["j"]
    triple = d["triple"]
    mu_a_vec, mu_b_vec = d["mu_a_vec"], d["mu_b_vec"]

    mean_raw = 0.5 * (e_a + e_b)
    sigma_site_total = float(np.sqrt(0.5 * (e_a - e_b).var(ddof=1)))

    series = site_series()
    sigma_static, modes, diag = disorder_and_modes(series, sigma_site_total)
    sigma_gap_static = np.sqrt(2.0) * sigma_static
    sigma_mean_static = sigma_static / np.sqrt(2.0)

    print("=== measured disorder decomposition (no fitting) ===")
    print(f"  per-site total (v5, 1 ns)        {diag['sigma_site_total_cm']:7.1f} cm^-1")
    print(f"  per-site fast (<100 fs, v6t2)    {diag['sigma_site_fast_cm']:7.1f}")
    print(f"  per-site slow (>100 fs)          {diag['sigma_site_slow_cm']:7.1f}")
    print(f"  + quasi-static (<= 2kT) part     {diag['sigma_site_quasi_static_cm']:7.1f}")
    print(f"  => static inhomogeneous          {diag['sigma_site_static_cm']:7.1f}")
    print(f"     applied as gap {np.sqrt(2)*sigma_static:.1f} (quadrature) and "
          f"mean {sigma_static/np.sqrt(2):.1f} (convolution), once each")
    print(f"  reorganisation energy, fast part {diag['lambda_site_fast_cm']:7.1f}")
    print(f"  spectral peak {diag['spectral_peak_cm']:.0f} cm^-1, "
          f"resolution {diag['resolution_cm']:.0f} cm^-1")
    print("  Franck-Condon modes read off the measured spectral density:")
    for w, s in modes:
        print(f"    omega = {w:7.1f} cm^-1   S = {s:.4f}")
    print(f"  homogeneous width: none (radiative 0.002 cm^-1, neglected)")

    vib_e, fc = vibrational_basis(modes, N_QUANTA)
    print(f"  vibrational basis {vib_e.size} states/site, "
          f"sum |<n|0>|^2 = {np.sum(fc**2):.5f}")

    nodes, weights = np.polynomial.hermite_e.hermegauss(N_HERMITE)
    weights = weights / weights.sum()

    shift_anchor = EXP_CENTER_CM - mean_raw.mean()
    centre_abinitio = mean_raw.mean() + STEOM_SHIFT_EV * CM_PER_EV

    # The per-frame mean site energy must NOT be carried through: the v5 frames
    # already sample the FULL mean-energy distribution (sd 387 cm^-1), fast part
    # included, so keeping it and then convolving the measured static width on
    # top counts the same disorder twice -- and over-broadens the band top, which
    # is exactly what sets a bisignate couplet's peak separation.  The mean
    # energy is therefore treated symmetrically with the gap: every frame sits at
    # the common centre, and the measured STATIC mean-energy spread is applied
    # once, by convolution.  Only the geometry-dependent quantities (J, T, the
    # two dipole vectors) stay per-frame.
    cases = {}
    for name, sig_gap in (("rigorous", sigma_gap_static),
                          ("degenerate_control", 0.0)):
        d_eps = np.zeros_like(GRID)
        eps = np.zeros_like(GRID)
        base = np.full_like(mean_raw, mean_raw.mean() + shift_anchor)
        for idx in range(j.size):
            for node, wt in zip(nodes, weights):
                delta = sig_gap * node
                nu, r_k, d_k = frame_sticks(
                    base[idx] + 0.5 * delta, base[idx] - 0.5 * delta, j[idx],
                    mu_a_vec[idx], mu_b_vec[idx], triple[idx], vib_e, fc)
                accumulate(nu, R_TO_DELTA_EPS * r_k * nu * wt, GRID, d_eps)
                accumulate(nu, D_TO_EPS * d_k * nu * wt, GRID, eps)
        d_eps /= j.size
        eps /= j.size
        d_eps = gaussian_convolve(d_eps, GRID, sigma_mean_static)
        eps = gaussian_convolve(eps, GRID, sigma_mean_static)
        cases[name] = (ELLIPTICITY_PER_DELTA_EPSILON * d_eps, eps)

    def metrics(theta, eps=None):
        hi, lo = int(np.argmax(theta)), int(np.argmin(theta))
        m = dict(peak_to_peak=float(theta.max() - theta.min()),
                 separation_cm=abs(float(GRID[hi] - GRID[lo])),
                 lobe_ratio=float(abs(theta.max() / theta.min())),
                 first_moment=float(np.trapezoid(
                     theta * (GRID - EXP_CENTER_CM) / GRID, GRID)),
                 nu_pos_cm=float(GRID[hi]), nu_neg_cm=float(GRID[lo]))
        if eps is not None:
            half = eps.max() / 2.0
            i = np.where(eps >= half)[0]
            m.update(eps_max=float(eps.max()),
                     eps_fwhm_cm=float(GRID[i[-1]] - GRID[i[0]]) if i.size else np.nan,
                     lambda_max_nm=float(1e7 / GRID[int(np.argmax(eps))]),
                     g_abs_peak=float(np.abs(theta).max()
                                      / ELLIPTICITY_PER_DELTA_EPSILON / eps.max()))
        return m

    exp_theta = sum(a / (1.0 + ((GRID - c) / w) ** 2) for a, c, w in EXP_BANDS)
    results = {k: metrics(*v) for k, v in cases.items()}
    results["experiment"] = metrics(exp_theta)
    results["experiment"]["eps_max"] = EPS_TD_EXP
    results["experiment"]["g_abs_peak"] = float(
        np.abs(exp_theta).max() / ELLIPTICITY_PER_DELTA_EPSILON / EPS_TD_EXP)

    print("\n  case                 pk-pk theta    sep   lobe   1st moment"
          "     eps_max   FWHM   g_abs")
    for k in ("rigorous", "degenerate_control", "experiment"):
        m = results[k]
        print(f"  {k:<20} {m['peak_to_peak']:10.3e} {m['separation_cm']:6.0f} "
              f"{m['lobe_ratio']:6.2f} {m['first_moment']:12.3e} "
              f"{m.get('eps_max', float('nan')):10.3e} "
              f"{m.get('eps_fwhm_cm', float('nan')):6.0f} {m['g_abs_peak']:8.2e}")

    r, e = results["rigorous"], results["experiment"]
    print(f"\n  predicted / observed   peak-to-peak {r['peak_to_peak']/e['peak_to_peak']:.2f}x"
          f",  separation {r['separation_cm']/e['separation_cm']:.2f}x"
          f",  g_abs {r['g_abs_peak']/e['g_abs_peak']:.2f}x")
    print(f"  band maximum {r['lambda_max_nm']:.1f} nm, absorption FWHM "
          f"{r['eps_fwhm_cm']:.0f} cm^-1")
    print(f"  ab initio origin would place it at {1e7/centre_abinitio:.1f} nm "
          f"({centre_abinitio - EXP_CENTER_CM:+.0f} cm^-1 from experiment)")

    payload = dict(disorder=diag,
                   modes=[dict(omega_cm=w, huang_rhys=s) for w, s in modes],
                   quadrature=dict(hermite_nodes=N_HERMITE, quanta=N_QUANTA,
                                   basis_states_per_site=int(vib_e.size),
                                   sigma_gap_static_cm=float(sigma_gap_static),
                                   sigma_mean_static_cm=float(sigma_mean_static)),
                   origin=dict(abinitio_centre_cm=float(centre_abinitio),
                               experimental_centre_cm=EXP_CENTER_CM),
                   cases=results)
    (OUT / "cd_rigorous.json").write_text(json.dumps(payload, indent=2))
    np.savez(OUT / "cd_rigorous.npz", grid_cm=GRID, experiment_theta=exp_theta,
             **{f"{k}_theta": v[0] for k, v in cases.items()},
             **{f"{k}_eps": v[1] for k, v in cases.items()})
    print(f"\n  wrote {OUT/'cd_rigorous.json'}")


if __name__ == "__main__":
    main()

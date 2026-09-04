#!/usr/bin/env python3
"""Quantum-corrected cumulant lineshape from the full v8dense C(t). No boundary.

Every previous treatment of the Venus band shape in this repo has cut the
measured fluctuation in two at a frequency chosen by hand:

  cd_rigorous.py      everything below 2kT is "not a mode on this timescale"
                      and becomes static disorder; everything above becomes a
                      Franck-Condon mode with S_i = sigma_i^2 / (2 kT omega_i).
  v8dense_cumulant.py the same cut at kT, reported as a variance fraction.

The cut is the problem. It is not a physical boundary -- a mode at 150 cm^-1 is
neither strictly static nor strictly quantum -- and the answer moves when the
cut moves, which is why the partition never settled. The second-order cumulant
removes the need for one. A mode contributes

    Re g(t) ~ (1/2) sigma_i^2 t^2      when omega_i t << 1   (static limit)
    Re g(t) ~ sigma_i^2 tau_i t        when omega_i t >> 1   (dephasing limit)

continuously and by itself, with no instruction from us. The whole spectral
density enters at once and each frequency is weighted by its own hbar*omega/kT.

Physics, with hbar = 1 and all frequencies in fs^-1:

    S(omega) = int dt e^{i omega t} C_cl(t)          classical, from the record
    J(omega) = (beta omega / 2) S(omega)             spectral density
    lambda   = (1/pi) int_0^inf domega J/omega = sigma^2 / (2 kT)   [checked]

    g(t) = (1/pi) int_0^inf domega (J/omega^2)
             [ coth(beta omega / 2) (1 - cos omega t) + i (sin omega t - omega t) ]

The classical cumulant is the coth -> 2/(beta omega) limit of the same integral
with the imaginary part dropped, so both are computed here from ONE spectral
density and differ only by the quantum statistics.

Im g is the part that matters and the part no Gaussian model has: it produces
the Stokes shift 2*lambda and the blue vibronic tail. The measured Venus band is
asymmetric (blue half-width / red = 1.78), and that asymmetry is the observable
this script can be tested against.

The one boundary that remains is the record length, and it is set by the
observation window rather than by us: 1 ps resolves to 33 cm^-1 and nothing
slower is present in C(t). That part is supplied as a static Gaussian from the
1 ns v5 ensemble, sigma_slow^2 = sigma_1ns^2 - sigma_1ps^2, and is measured to
be consistent with zero.

Read-only apart from results/v8dense_quantum_cumulant.{json,npz}.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent

C_CM_S = 2.99792458e10
TWO_PI_C = 2.0 * np.pi * 2.99792458e-5      # cm^-1 -> fs^-1
KT_CM = 208.5                                # 300 K
DT_FS = 2.0

# experiment, digitised from Nguyen Fig 2e (peak, FWHM, asymmetry) and the
# Venus emission maximum at 528 nm (Nagai 2002) for the Stokes shift
EXP_PEAK_CM = 19341.0
EXP_FWHM_CM = 972.0
EXP_ASYM = 1.78
EXP_EMISSION_NM = 528.0
# blue-wing intensities of the digitised band, relative to its maximum
EXP_WING = {"+800": 0.420, "+1200": 0.364, "+1600": 0.264, "+2000": 0.143}
EXP_STOKES_CM = EXP_PEAK_CM - 1.0e7 / EXP_EMISSION_NM


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def dense_series():
    """The 1 ps / 2 fs site energies, in cm^-1, ordered and gap-checked."""
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/v8dense_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    rows.sort(key=lambda r: int(r["frame_index"]))
    idx = np.array([int(r["frame_index"]) for r in rows])
    if not np.all(np.diff(idx) == 1):
        raise RuntimeError("v8dense record has gaps; C(t) would be wrong")
    a = np.array([float(r["site_A_energy_cm-1"]) for r in rows])
    b = np.array([float(r["site_B_energy_cm-1"]) for r in rows])
    return a, b


def production_sigma():
    """Per-site sigma over the 1 ns v5 ensemble -- the slow reference."""
    rows = [r for w in "abc"
            for r in csv.DictReader((ROOT / f"manifests/v5camp_{w}_frames.csv").open())
            if r["convergence_status"] == "complete"]
    a = np.array([float(r["site_A_energy_cm-1"]) for r in rows])
    b = np.array([float(r["site_B_energy_cm-1"]) for r in rows])
    return a, b


# --------------------------------------------------------------------------
# spectral density
# --------------------------------------------------------------------------
def power_spectrum(x_cm, dt_fs=DT_FS):
    """Two-sided PSD S(omega) in fs^-1 with (1/2pi) int S domega = sigma^2.

    Single Hann taper on the full record: 33 cm^-1 resolution, chi^2_2 noise per
    bin. The noise does not matter because every use of S below is an integral
    over omega, and the Parseval check confirms the normalisation exactly.
    """
    x = (np.asarray(x_cm, float) - np.mean(x_cm)) * TWO_PI_C     # fs^-1
    n = x.size
    w = np.hanning(n)
    xw = x * w
    f = np.fft.rfft(xw) * dt_fs                    # approximates int dt e^{iwt}
    s = np.abs(f) ** 2 / (dt_fs * np.sum(w ** 2))  # fs^-1
    omega = 2.0 * np.pi * np.fft.rfftfreq(n, dt_fs)  # fs^-1
    return omega, s, float(np.sum(xw ** 2) / np.sum(w ** 2))     # + tapered var


def cumulant_g(t_fs, omega, s, beta_fs, quantum=True):
    """g(t) by direct spectral quadrature. omega[0] = 0 handled analytically.

    Re integrand -> S(0) t^2 / beta and Im -> 0 as omega -> 0, both finite, so
    the DC bin is included at its limit rather than dropped.
    """
    t = np.asarray(t_fs, float)[:, None]
    w = np.asarray(omega, float)[None, :]
    sv = np.asarray(s, float)[None, :]
    nz = w > 0.0

    if quantum:
        # coth(beta w / 2) / w, finite-safe
        x = np.where(nz, beta_fs * w / 2.0, 1.0)
        kernel = np.where(nz, np.cosh(x) / np.sinh(x) / np.where(nz, w, 1.0), 0.0)
    else:
        kernel = np.where(nz, 2.0 / (beta_fs * np.where(nz, w, 1.0) ** 2), 0.0)

    re_int = sv * kernel * (1.0 - np.cos(w * t))
    re_int = np.where(nz, re_int, sv * t ** 2 / beta_fs)     # omega -> 0 limit

    if quantum:
        im_int = np.where(nz, sv * (np.sin(w * t) - w * t) / np.where(nz, w, 1.0),
                          0.0)
    else:
        im_int = np.zeros_like(re_int)

    pref = beta_fs / (2.0 * np.pi)
    g = pref * (np.trapezoid(re_int, omega, axis=1)
                + 1j * np.trapezoid(im_int, omega, axis=1))
    return g


def lineshape(t_fs, g, grid_cm, centre_cm, sigma_static_cm=0.0):
    """I(nu) = Re int_0^inf dt exp[i(nu - nu_vert + lambda) t - g(t)].

    centre_cm is the VERTICAL (mean) transition; Im g supplies the shift from
    it, so the caller passes nu_vert and lambda is already inside g.
    """
    t = np.asarray(t_fs, float)
    decay = np.exp(-g)
    if sigma_static_cm > 0.0:                       # static Gaussian, in time
        decay = decay * np.exp(-0.5 * (sigma_static_cm * TWO_PI_C * t) ** 2)
    dw = (np.asarray(grid_cm, float) - centre_cm) * TWO_PI_C
    phase = np.exp(1j * dw[:, None] * t[None, :])
    return np.real(np.trapezoid(phase * decay[None, :], t, axis=1)) / np.pi


def emission_lineshape(t_fs, g, lam_fs, grid_cm, centre_cm, sigma_static_cm=0.0):
    """Fluorescence from the SAME g(t):  exp[i(w - w_eg + 2 lambda) t - g*(t)].

    Needed because 2*lambda is NOT the peak-to-peak Stokes shift. A mode with
    hbar*omega >> kT and S < 1 has its 0-0 line as the peak of BOTH bands, so it
    contributes its full lambda to the first moment and nothing at all to the
    separation of the maxima. Comparing lambda directly against the measured
    515/528 nm separation therefore compares different quantities; this computes
    the separation the model actually predicts.
    """
    t = np.asarray(t_fs, float)
    decay = np.exp(-np.conj(g) + 2.0j * lam_fs * t)
    if sigma_static_cm > 0.0:
        decay = decay * np.exp(-0.5 * (sigma_static_cm * TWO_PI_C * t) ** 2)
    dw = (np.asarray(grid_cm, float) - centre_cm) * TWO_PI_C
    phase = np.exp(1j * dw[:, None] * t[None, :])
    return np.real(np.trapezoid(phase * decay[None, :], t, axis=1)) / np.pi


def first_moment(grid, y):
    """<nu> of a (possibly slightly negative-tailed) lineshape."""
    w = np.clip(y, 0.0, None)
    return float(np.trapezoid(w * grid, grid) / np.trapezoid(w, grid))


def band_metrics(grid, y):
    """Peak, FWHM and blue/red half-width ratio of a lineshape."""
    y = np.asarray(y, float)
    i = int(y.argmax())
    half = 0.5 * y[i]
    left = np.where(y[:i] <= half)[0]
    right = np.where(y[i:] <= half)[0]
    if left.size == 0 or right.size == 0:
        return dict(peak_cm=float(grid[i]), fwhm_cm=float("nan"),
                    asymmetry=float("nan"))
    lo = float(np.interp(half, [y[left[-1]], y[left[-1] + 1]],
                         [grid[left[-1]], grid[left[-1] + 1]]))
    j = i + right[0]
    hi = float(np.interp(half, [y[j], y[j - 1]], [grid[j], grid[j - 1]]))
    pk = float(grid[i])
    return dict(peak_cm=pk, fwhm_cm=hi - lo, red_hwhm_cm=pk - lo,
                blue_hwhm_cm=hi - pk,
                asymmetry=(hi - pk) / (pk - lo) if pk > lo else float("nan"))


def acf(x):
    """Unbiased autocorrelation. Reliable only at short lag -- see tau_c below."""
    x = np.asarray(x, float) - np.mean(x)
    n = x.size
    f = np.fft.rfft(x, 2 * n)
    c = np.fft.irfft(f * np.conj(f), 2 * n)[:n].real
    return c / np.arange(n, 0, -1)


def tau_correlation(c, dt_fs=DT_FS):
    """tau_c integrated to the FIRST ZERO CROSSING, not to the record end.

    The running integral of this C(t) does not converge: it reaches 6.1 fs at
    20 fs, 22.5 fs at 500 fs and -96.6 fs at 1 ps, because beyond the decay the
    unbiased estimator is a random walk whose step size grows as 1/sqrt(N-k).
    Integrating to the end of the record therefore reports wherever that walk
    happened to be. Truncating at the first zero crossing is the standard fix
    and is what the initial decay actually supports.
    """
    cn = c / c[0]
    zero = np.where(cn <= 0.0)[0]
    stop = int(zero[0]) if zero.size else cn.size // 2
    lag = np.arange(stop + 1) * dt_fs
    return float(np.trapezoid(cn[: stop + 1], lag)), stop * dt_fs


def implied_acf(s, n, dt_fs=DT_FS):
    """C(t) corresponding exactly to the estimated PSD, for the closure test."""
    full = np.concatenate([s, s[-2:0:-1]])
    return np.fft.irfft(full, n=full.size).real[:n] / dt_fs


def band_table(omega, s):
    """Fraction of the variance in each frequency band, one-sided."""
    nu = omega / TWO_PI_C
    w = s.copy()
    w[1:-1] *= 2.0
    tot = float(np.trapezoid(w, omega))
    rows = []
    for lo, hi in ((0.0, KT_CM), (KT_CM, 2 * KT_CM), (2 * KT_CM, 900.0),
                   (900.0, 1500.0), (1500.0, 2200.0), (2200.0, float(nu[-1]) + 1)):
        m = (nu >= lo) & (nu < hi)
        if not m.any():
            continue
        rows.append(dict(lo_cm=lo, hi_cm=min(hi, float(nu[-1])),
                         variance_fraction=float(np.trapezoid(w[m], omega[m]) / tot)))
    return rows, tot


def solve_static(t, g, grid, centre, target_fwhm, hi=900.0):
    """sigma_static that brings the cumulant band to target_fwhm, by bisection.

    The naive quadrature estimate sqrt(FWHM_exp^2 - FWHM_cum^2)/2.3548 is wrong
    here -- the cumulant band is not Gaussian, so widths do not add in
    quadrature -- and it overstates the requirement by ~20%.
    """
    lo = 0.0
    if band_metrics(grid, lineshape(t, g, grid, centre, hi))["fwhm_cm"] < target_fwhm:
        return float("nan")
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        w = band_metrics(grid, lineshape(t, g, grid, centre, mid))["fwhm_cm"]
        lo, hi = (mid, hi) if w < target_fwhm else (lo, mid)
    return 0.5 * (lo + hi)


def wing_profile(grid, y, offsets=(800.0, 1200.0, 1600.0, 2000.0)):
    """Blue-wing intensity relative to the maximum, peak-aligned.

    The half-max metrics (FWHM, blue/red) are blind to a stiff mode: with
    S ~ 0.19 at 1550 cm^-1 the sideband sits well below half height. The wing is
    where those modes actually appear, so it is the observable that tests them.
    """
    y = np.asarray(y, float)
    pk = grid[y.argmax()]
    ymax = y.max()
    return {f"+{int(o)}": float(np.interp(pk + o, grid, y) / ymax) for o in offsets}


def band_decomposition(t, omega, s, beta_fs, grid, centre):
    """Which part of J(omega) actually sets the visible band shape.

    lambda integrates the whole spectral density, but the FWHM and the blue/red
    asymmetry are half-max properties: a mode at 1550 cm^-1 with S = 0.19 puts
    its sideband far outside the half-max crossings and cannot move them. The
    two observables therefore probe different parts of the same C(t).
    """
    nu = omega / TWO_PI_C
    rows = []
    for name, keep in (("full", np.ones_like(nu, bool)),
                       ("omega < 417", nu < 2 * KT_CM),
                       ("417-900", (nu >= 2 * KT_CM) & (nu < 900)),
                       ("900-1500", (nu >= 900) & (nu < 1500)),
                       ("1500-2200", (nu >= 1500) & (nu < 2200)),
                       ("all except 1500-2200", ~((nu >= 1500) & (nu < 2200)))):
        ss = np.where(keep, s, 0.0)
        sig = float(np.sqrt(2 * np.trapezoid(ss, omega) / (2 * np.pi)) / TWO_PI_C)
        lam = float(np.trapezoid(0.5 * beta_fs * ss[omega > 0],
                                 omega[omega > 0]) / np.pi) / TWO_PI_C
        m = band_metrics(grid, lineshape(t, cumulant_g(t, omega, ss, beta_fs, True),
                                         grid, centre))
        rows.append(dict(part=name, sigma_cm=sig, lambda_cm=lam,
                         fwhm_cm=m["fwhm_cm"], asymmetry=m["asymmetry"]))
    return rows


def main() -> None:
    a_dense, b_dense = dense_series()
    a_prod, b_prod = production_sigma()
    n = a_dense.size
    beta_fs = 1.0 / (KT_CM * TWO_PI_C)
    res_cm = 1.0 / (n * DT_FS * 1e-15 * C_CM_S)

    print("=== quantum-corrected cumulant, full C(t), no frequency boundary ===")
    print(f"  record {n} x {DT_FS:.0f} fs = {n*DT_FS/1000:.2f} ps, resolution "
          f"{res_cm:.0f} cm^-1;  kT = {KT_CM:.1f} cm^-1")
    print(f"  experiment: peak {EXP_PEAK_CM:.0f}, FWHM {EXP_FWHM_CM:.0f}, "
          f"blue/red {EXP_ASYM:.2f}, Stokes {EXP_STOKES_CM:.0f} cm^-1\n")

    t = np.arange(0.0, 400.0, 0.5)
    grid = np.arange(15000.0, 25000.0, 2.0)
    out = {"n_frames": n, "dt_fs": DT_FS, "resolution_cm": res_cm, "kT_cm": KT_CM,
           "experiment": {"peak_cm": EXP_PEAK_CM, "fwhm_cm": EXP_FWHM_CM,
                          "asymmetry_blue_over_red": EXP_ASYM,
                          "stokes_shift_cm": EXP_STOKES_CM,
                          "emission_nm": EXP_EMISSION_NM},
           "sites": {}}
    curves = {"grid_cm": grid}

    for tag, x, xprod in (("site_A", a_dense, a_prod), ("site_B", b_dense, b_prod)):
        omega, s, var_tapered = power_spectrum(x)
        sigma_cm = float(np.sqrt(var_tapered) / TWO_PI_C)

        # ---- check 1: Parseval --------------------------------------------
        parseval = float(2.0 * np.trapezoid(s, omega) / (2.0 * np.pi))

        # ---- check 2: lambda from J equals the classical sigma^2 / 2kT ------
        nz = omega > 0
        lam_cm = float(np.trapezoid(0.5 * beta_fs * s[nz], omega[nz])
                       / np.pi * 2.0 / 2.0) / TWO_PI_C
        lam_classical_cm = sigma_cm ** 2 / (2.0 * KT_CM)

        # ---- check 3: spectral g == time-domain g on the SAME C(t) ---------
        g_cl = cumulant_g(t, omega, s, beta_fs, quantum=False)
        c_imp = implied_acf(s, n)
        lag = np.arange(n) * DT_FS
        tsub = t[t <= 60.0]
        g_td = np.array([np.trapezoid((ti - lag[lag <= ti]) * c_imp[lag <= ti],
                                      lag[lag <= ti]) for ti in tsub])
        closure = float(np.max(np.abs(g_td - g_cl[: tsub.size].real))
                        / max(np.max(np.abs(g_td)), 1e-30))

        # ---- check 4: the spectrum is vibrational, not per-frame noise ------
        nu = omega / TWO_PI_C
        w1 = s.copy(); w1[1:-1] *= 2.0
        hi = nu >= 2500.0
        frac_hi = float(np.trapezoid(w1[hi], omega[hi]) / np.trapezoid(w1, omega))
        white_hi = float((nu[-1] - 2500.0) / nu[-1])

        c = acf(x)
        tau_c, tau_cut = tau_correlation(c)
        bands, _ = band_table(omega, s)

        # ---- the only surviving partition, set by the record length --------
        sigma_prod = float(np.std(xprod, ddof=1))
        se_prod = sigma_prod / np.sqrt(2.0 * xprod.size)
        n_eff = max(n * DT_FS / (2.0 * tau_c), 1.0)
        se_dense = sigma_cm / np.sqrt(2.0 * n_eff)
        var_slow = sigma_prod ** 2 - sigma_cm ** 2
        se_var = float(np.hypot(2 * sigma_prod * se_prod, 2 * sigma_cm * se_dense))
        sigma_slow = float(np.sqrt(max(var_slow, 0.0)))
        slow_hi = float(np.sqrt(max(var_slow + 1.645 * se_var, 0.0)))

        g_q = cumulant_g(t, omega, s, beta_fs, quantum=True)
        i_q = lineshape(t, g_q, grid, EXP_PEAK_CM)
        i_cl = lineshape(t, g_cl, grid, EXP_PEAK_CM)
        i_qs = lineshape(t, g_q, grid, EXP_PEAK_CM, sigma_static_cm=sigma_slow)
        i_qhi = lineshape(t, g_q, grid, EXP_PEAK_CM, sigma_static_cm=slow_hi)
        curves[f"{tag}_quantum"] = i_q
        curves[f"{tag}_classical"] = i_cl
        curves[f"{tag}_quantum_static"] = i_qs
        curves[f"{tag}_psd_nu_cm"] = nu
        curves[f"{tag}_psd"] = w1

        m_q, m_cl = band_metrics(grid, i_q), band_metrics(grid, i_cl)
        m_qs, m_qhi = band_metrics(grid, i_qs), band_metrics(grid, i_qhi)
        # The one honest test of the band SHAPE: fix sigma_static by the width
        # alone -- and require its value to fall inside the range the 1 ns
        # ensemble independently permits -- then the asymmetry is a free
        # prediction that nothing in the model was tuned to.
        # predicted absorption/emission separation from the same g(t)
        lam_fs = lam_cm * TWO_PI_C
        wide = np.arange(13000.0, 26000.0, 2.0)
        abs_q = lineshape(t, g_q, wide, EXP_PEAK_CM)
        em_q = emission_lineshape(t, g_q, lam_fs, wide, EXP_PEAK_CM)
        stokes_peak = float(wide[abs_q.argmax()] - wide[em_q.argmax()])
        stokes_moment = float(first_moment(wide, abs_q) - first_moment(wide, em_q))

        sigma_needed = solve_static(t, g_q, grid, EXP_PEAK_CM, EXP_FWHM_CM)
        m_width = band_metrics(grid, lineshape(t, g_q, grid, EXP_PEAK_CM,
                                               sigma_needed))
        s_cl_need = solve_static(t, g_cl, grid, EXP_PEAK_CM, EXP_FWHM_CM)
        m_width_cl = band_metrics(grid, lineshape(t, g_cl, grid, EXP_PEAK_CM,
                                                  s_cl_need))
        y_w = lineshape(t, g_q, grid, EXP_PEAK_CM, sigma_needed)
        y_wcl = lineshape(t, g_cl, grid, EXP_PEAK_CM, s_cl_need)
        curves[f"{tag}_width_matched"] = y_w
        wing_q, wing_cl = wing_profile(grid, y_w), wing_profile(grid, y_wcl)
        decomp = band_decomposition(t, omega, s, beta_fs, grid, EXP_PEAK_CM)

        # calibrate the ONE amplitude against an independent observable -- the
        # emission maximum, not the band being predicted -- and re-predict.
        kappa = (EXP_STOKES_CM / 2.0) / lam_cm
        g_cal = cumulant_g(t, omega, kappa * s, beta_fs, quantum=True)
        m_cal = band_metrics(grid, lineshape(t, g_cal, grid, EXP_PEAK_CM))

        rec = dict(
            sigma_dense_cm=sigma_cm, sigma_dense_se_cm=float(se_dense),
            tau_c_fs=tau_c, tau_c_truncation_fs=tau_cut, n_eff=float(n_eff),
            kubo_sigma_tau=float(sigma_cm * TWO_PI_C * tau_c),
            sigma_production_1ns_cm=sigma_prod,
            sigma_slow_cm=sigma_slow, sigma_slow_upper95_cm=slow_hi,
            sigma_static_needed_for_width_cm=sigma_needed,
            lambda_cm=lam_cm, lambda_classical_check_cm=lam_classical_cm,
            stokes_two_lambda_cm=2.0 * lam_cm,
            stokes_peak_to_peak_predicted_cm=stokes_peak,
            stokes_first_moment_predicted_cm=stokes_moment,
            stokes_ratio_peak_over_exp=stokes_peak / EXP_STOKES_CM,
            variance_bands=bands,
            checks=dict(parseval_ratio=float(parseval / var_tapered),
                        lambda_ratio=float(lam_cm / lam_classical_cm),
                        spectral_vs_timedomain_relerr=closure,
                        variance_above_2500_cm=frac_hi,
                        white_noise_would_give=white_hi),
            quantum=m_q, classical=m_cl,
            quantum_plus_measured_static=m_qs,
            quantum_plus_static_upper95=m_qhi,
            stokes_calibrated=dict(kappa=float(kappa),
                                   sigma_cm=float(sigma_cm * np.sqrt(kappa)),
                                   **m_cal),
            band_decomposition=decomp,
            width_matched=dict(sigma_static_cm=float(sigma_needed),
                               permitted_upper95_cm=slow_hi,
                               inside_permitted=bool(sigma_needed <= slow_hi),
                               **m_width),
            width_matched_classical=dict(sigma_static_cm=float(s_cl_need),
                                         **m_width_cl),
            blue_wing=dict(quantum=wing_q, classical=wing_cl,
                           experiment=EXP_WING))
        out["sites"][tag] = rec

        print(f"  --- {tag} ---")
        print(f"    sigma (1 ps) {sigma_cm:6.1f} +/- {se_dense:4.1f}   "
              f"sigma (1 ns) {sigma_prod:6.1f} +/- {se_prod:4.1f}   "
              f"=> sigma_slow {sigma_slow:5.1f} (95% upper {slow_hi:.0f})")
        print(f"    tau_c {tau_c:4.1f} fs (truncated at the zero crossing, "
              f"{tau_cut:.0f} fs)   sigma*tau_c {rec['kubo_sigma_tau']:.2f}")
        print("    variance by band:  " + "  ".join(
            f"{b['lo_cm']:.0f}-{b['hi_cm']:.0f}:{100*b['variance_fraction']:.0f}%"
            for b in bands))
        print(f"    checks  Parseval {rec['checks']['parseval_ratio']:.6f}"
              f"   lambda/(s^2/2kT) {rec['checks']['lambda_ratio']:.4f}"
              f"   g closure {closure:.1e}"
              f"   var>2500 {100*frac_hi:.2f}% (white would be {100*white_hi:.0f}%)")
        print(f"    lambda {lam_cm:6.1f}   2*lambda (first moment) "
              f"{2*lam_cm:5.0f}  model moment {stokes_moment:5.0f}")
        print(f"    Stokes PEAK-TO-PEAK predicted {stokes_peak:5.0f} cm^-1  vs "
              f"measured {EXP_STOKES_CM:.0f}   "
              f"({rec['stokes_ratio_peak_over_exp']:.2f}x)  "
              f"<- the like-for-like comparison")
        print(f"    classical cumulant            FWHM {m_cl['fwhm_cm']:5.0f}  "
              f"blue/red {m_cl['asymmetry']:.2f}")
        print(f"    QUANTUM cumulant              FWHM {m_q['fwhm_cm']:5.0f}  "
              f"blue/red {m_q['asymmetry']:.2f}")
        print(f"    + measured static {sigma_slow:5.1f}       FWHM "
              f"{m_qs['fwhm_cm']:5.0f}  blue/red {m_qs['asymmetry']:.2f}")
        print(f"    + static at 95% upper {slow_hi:3.0f}     FWHM "
              f"{m_qhi['fwhm_cm']:5.0f}  blue/red {m_qhi['asymmetry']:.2f}")
        print(f"    EXPERIMENT                    FWHM {EXP_FWHM_CM:5.0f}  "
              f"blue/red {EXP_ASYM:.2f}"
              f"")
        print(f"    >> fix sigma_static by the WIDTH alone: {sigma_needed:.0f} cm^-1 "
              f"(permitted 0-{slow_hi:.0f}, "
              f"{'inside' if sigma_needed <= slow_hi else 'OUTSIDE'})")
        print(f"       asymmetry is then a FREE prediction: quantum "
              f"{m_width['asymmetry']:.2f}   classical {m_width_cl['asymmetry']:.2f}"
              f"   experiment {EXP_ASYM:.2f}")
        print("       blue wing (also free):  " + "   ".join(
            f"{k} q{wing_q[k]:.3f}/c{wing_cl[k]:.3f}/exp{EXP_WING[k]:.3f}"
            for k in wing_q))
        print(f"    Stokes-calibrated (kappa={kappa:.3f})  FWHM "
              f"{m_cal['fwhm_cm']:5.0f}  blue/red {m_cal['asymmetry']:.2f}"
              f"   <- matching lambda destroys the asymmetry")
        print("    which part of J(omega) sets the visible band:")
        for r in decomp:
            print(f"      {r['part']:<22s} sigma {r['sigma_cm']:5.0f}  "
                  f"lambda {r['lambda_cm']:5.0f}  FWHM {r['fwhm_cm']:5.0f}  "
                  f"blue/red {r['asymmetry']:.2f}")
        print()

    np.savez(ROOT / "results" / "v8dense_quantum_cumulant.npz", **curves)
    (ROOT / "results" / "v8dense_quantum_cumulant.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print(f"  wrote results/v8dense_quantum_cumulant.json and .npz")


if __name__ == "__main__":
    main()

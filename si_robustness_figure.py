#!/usr/bin/env python3
r"""
si_robustness_figure.py  --  Supporting-Information figure:
"Origin and robustness of the near-field excitonic coupling."

Composes, in one paper-styled 3-panel figure, the three sensitivity/decomposition
robustness results (each already produced piecemeal elsewhere):

  (a) Multipole convergence of J: partial sums of the primitive Cartesian
      multipole series (dip-dip -> +dip-quad -> +quad-quad -> +dip-oct) as a
      fraction of the full transition-density coupling J_TDC, showing the PDA
      and multipole-truncated models miss the bulk of the near-field coupling.
      Source: multipole_out/multipole_decomposition.csv (multipole_analysis.py).

  (b) Dephasing sensitivity: coherence decay |rho_12(t)| for T2* over 20-200 fs;
      every 1/e time is orders of magnitude below the ~8.3 ps Debye time, so the
      exciton-formation / dephasing timescale separation is preserved throughout.
      Physics reused from open_quantum_dynamics.py (solve_me).

  (c) Dielectric invariance: J(t) for static permittivity eps_s over 4-78. The
      t=0 coupling J(0) is set by eps_opt alone (curves coincide at t=0); only
      the fully relaxed long-time limit declines with eps_s.
      Physics reused from open_quantum_dynamics.py (J_of_t).

Output: si_out/Fig_SI_Robustness.pdf   (and the three panels separately).
"""
import argparse
import csv
from pathlib import Path

import numpy as np

import open_quantum_dynamics as oqd

# Palette consistent with the main figures.
_C_BAR = "#1f6aa5"
_C_CUM = "#c0392b"
_C_MEAN = "#222222"
FIGSIZE = (4.2, 3.4)


def _read_multipole(csv_path):
    terms, cum, tdc, pda = [], [], None, None
    with open(csv_path) as f:
        for row in csv.reader(f):
            if not row or row[0] == "term":
                continue
            name = row[0]
            if name == "full_TDC":
                tdc = float(row[2])
            elif name == "PDA_reference":
                pda = float(row[2])
            elif name in ("multipole_sum", "TDC_over_PDA"):
                continue
            else:
                terms.append(name)
                cum.append(float(row[2]))
    return terms, np.asarray(cum), tdc, pda


def _panel_multipole(plt, out, csv_path):
    terms, cum, tdc, pda = _read_multipole(csv_path)
    frac = 100.0 * cum / tdc
    labels = ["dip-dip\n(PDA)", "+dip-quad", "+quad-quad", "+dip-oct"]
    labels = labels[:len(terms)]
    x = np.arange(len(terms))
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(x, frac, color=_C_BAR, alpha=0.85, width=0.62, label="cumulative multipole")
    ax.axhline(100.0, color=_C_MEAN, ls="--", lw=1.4, label=r"full $J_{\mathrm{TDC}}$")
    for xi, fr in zip(x, frac):
        ax.text(xi, fr + 2, f"{fr:.0f}%", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(r"% of full $J_{\mathrm{TDC}}$")
    ax.set_ylim(0, 115)
    ax.grid(alpha=0.25, axis="y", lw=0.6)
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
              loc="upper left", fontsize=8)
    fig.tight_layout()
    p = out / "Fig_SI_Multipole.pdf"
    fig.savefig(p)
    plt.close(fig)
    return p, frac


def _panel_t2(plt, out, t2_list_fs, tf, dt):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    cmap = plt.cm.viridis(np.linspace(0.0, 0.85, len(t2_list_fs)))
    efolds = []
    for c, t2_fs in zip(cmap, t2_list_fs):
        p = oqd.make_params(t2_star=t2_fs / 1000.0)
        me = oqd.solve_me(p, tf, dt)
        coh = me["coh"] / (me["coh"][0] + 1e-30)
        ax.plot(me["t"] * 1000.0, coh, color=c, lw=1.7, label=fr"${t2_fs:.0f}$")
        below = np.where(me["coh"] <= me["coh"][0] / np.e)[0]
        efolds.append(me["t"][below[0]] * 1000.0 if below.size else np.nan)
    ax.set_xlim(0, 300)
    ax.set_xlabel("time (fs)")
    ax.set_ylabel(r"$|\rho_{12}(t)|/|\rho_{12}(0)|$")
    ax.grid(alpha=0.25, lw=0.6)
    leg = ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
                    loc="upper right", fontsize=8, ncol=2, title=r"$T_2^*$ (fs)")
    leg.get_title().set_fontsize(8)
    ax.text(0.97, 0.55,
            r"all $\tau_{1/e}\ll\tau_D=8.3$ ps",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    fig.tight_layout()
    p = out / "Fig_SI_T2.pdf"
    fig.savefig(p)
    plt.close(fig)
    return p, efolds


def _panel_eps(plt, out, eps_list, tf, dt):
    t = np.arange(0.0, tf + dt, dt)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    cmap = plt.cm.plasma(np.linspace(0.0, 0.82, len(eps_list)))
    j0s = []
    for c, es in zip(cmap, eps_list):
        p = oqd.make_params(eps_s=es)
        Jt = oqd.J_of_t(t, p)
        j0s.append(Jt[0])
        ax.plot(t, Jt, color=c, lw=1.7, label=fr"${es:g}$")
    ax.axhline(oqd.J_OPT, color=_C_MEAN, ls=":", lw=1.4,
               label=fr"$J(0)={oqd.J_OPT:.0f}$")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"$J(t)$ (cm$^{-1}$)")
    ax.grid(alpha=0.25, lw=0.6)
    leg = ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
                    loc="upper right", fontsize=8, ncol=2, title=r"$\varepsilon_s$")
    leg.get_title().set_fontsize(8)
    fig.tight_layout()
    p = out / "Fig_SI_Dielectric.pdf"
    fig.savefig(p)
    plt.close(fig)
    return p, (max(j0s) - min(j0s))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--multipole-csv", type=Path,
                    default=Path("multipole_out/multipole_decomposition.csv"))
    ap.add_argument("--t2-list", type=float, nargs="+", default=[20, 40, 60, 100, 200])
    ap.add_argument("--eps-list", type=float, nargs="+", default=[4, 10, 20, 40, 78])
    ap.add_argument("--tf", type=float, default=oqd.TF)
    ap.add_argument("--dt", type=float, default=oqd.DT)
    ap.add_argument("--out", type=Path, default=Path("si_out"))
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    plt = oqd._mpl()

    p_a, frac = _panel_multipole(plt, args.out, args.multipole_csv)
    p_b, efolds = _panel_t2(plt, args.out, args.t2_list, args.tf, args.dt)
    p_c, spread = _panel_eps(plt, args.out, args.eps_list, args.tf, args.dt)

    # Composed 3-panel figure.
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.4))
    for ax, lab in zip(axes, ("a", "b", "c")):
        ax.set_title(f"({lab})", loc="left", fontsize=12)

    # (a) multipole
    terms, cum, tdc, pda = _read_multipole(args.multipole_csv)
    fr = 100.0 * cum / tdc
    x = np.arange(len(terms))
    axes[0].bar(x, fr, color=_C_BAR, alpha=0.85, width=0.62)
    axes[0].axhline(100.0, color=_C_MEAN, ls="--", lw=1.4)
    for xi, f in zip(x, fr):
        axes[0].text(xi, f + 2, f"{f:.0f}%", ha="center", va="bottom", fontsize=8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["dip-dip", "+d-q", "+q-q", "+d-o"][:len(terms)], fontsize=8)
    axes[0].set_ylabel(r"% of full $J_{\mathrm{TDC}}$"); axes[0].set_ylim(0, 115)
    axes[0].grid(alpha=0.25, axis="y", lw=0.6)

    # (b) T2*
    cmap = plt.cm.viridis(np.linspace(0.0, 0.85, len(args.t2_list)))
    for c, t2_fs in zip(cmap, args.t2_list):
        p = oqd.make_params(t2_star=t2_fs / 1000.0)
        me = oqd.solve_me(p, args.tf, args.dt)
        axes[1].plot(me["t"] * 1000.0, me["coh"] / (me["coh"][0] + 1e-30),
                     color=c, lw=1.6, label=fr"${t2_fs:.0f}$")
    axes[1].set_xlim(0, 300); axes[1].set_xlabel("time (fs)")
    axes[1].set_ylabel(r"$|\rho_{12}|/|\rho_{12}(0)|$")
    axes[1].grid(alpha=0.25, lw=0.6)
    lg = axes[1].legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
                        loc="upper right", fontsize=7, ncol=2, title=r"$T_2^*$ (fs)")
    lg.get_title().set_fontsize(7)

    # (c) dielectric
    t = np.arange(0.0, args.tf + args.dt, args.dt)
    cmap = plt.cm.plasma(np.linspace(0.0, 0.82, len(args.eps_list)))
    for c, es in zip(cmap, args.eps_list):
        axes[2].plot(t, oqd.J_of_t(t, oqd.make_params(eps_s=es)), color=c, lw=1.6,
                     label=fr"${es:g}$")
    axes[2].axhline(oqd.J_OPT, color=_C_MEAN, ls=":", lw=1.4)
    axes[2].set_xlabel("time (ps)"); axes[2].set_ylabel(r"$J(t)$ (cm$^{-1}$)")
    axes[2].grid(alpha=0.25, lw=0.6)
    lg = axes[2].legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
                        loc="upper right", fontsize=7, ncol=2, title=r"$\varepsilon_s$")
    lg.get_title().set_fontsize(7)

    fig.tight_layout()
    composed = args.out / "Fig_SI_Robustness.pdf"
    fig.savefig(composed)
    plt.close(fig)

    print("=" * 60)
    print(f"  (a) multipole cumulative % of J_TDC: "
          f"{', '.join(f'{f:.0f}' for f in fr)}  (PDA={100*pda/tdc:.0f}%)")
    print(f"  (b) T2* 1/e coherence times (fs): "
          f"{', '.join(f'{t2:.0f}->{e:.0f}' for t2, e in zip(args.t2_list, efolds))}"
          f"   (all << 8300 fs Debye)")
    print(f"  (c) J(0) spread over eps_s {args.eps_list}: {spread:.2e} cm^-1 (invariant)")
    print("=" * 60)
    for k, v in {"multipole(a)": p_a, "T2(b)": p_b, "dielectric(c)": p_c,
                 "composed": composed}.items():
        print(f"  {k:13s}: {v}")


if __name__ == "__main__":
    main()

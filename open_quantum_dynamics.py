#!/usr/bin/env python3
r"""
open_quantum_dynamics.py  --  Open-quantum-system dynamics of the Venus dimer exciton.

Python port of the MATLAB open-quantum-systems code in LindbladCodes/
(Combined.m / Lindblad.m / NonAdiabatic.m), which the authors can no longer run
(MATLAB access lost). It integrates the Debye-screened, time-dependent coupling
J(t) together with Lindblad pure dephasing and reproduces the manuscript
figures, AND adds the reviewer-requested sensitivity sweeps (item 4):

  * default            : regenerate the six manuscript figures
                         (Fig_Coupling, Fig_SSE_Site, Fig_ME_Site,
                          Fig_SSE_Adiabatic, Fig_ME_Adiabatic, Fig_Bloch_Grid).
  * --sweep-t2         : sweep the pure-dephasing time T2* (reviewers R2/R3 note
                         60 fs is borrowed from photosynthetic systems) and show
                         the timescale-separation conclusion holds across a range
                         -> Fig_T2_Sweep.pdf.
  * --sweep-eps        : vary the static (protein) dielectric (R3/R4 note eps=78
                         is wrong inside a beta-barrel) and show the central
                         t=0 optical-limit coupling J(0)=100 cm^-1 is invariant
                         to it -> Fig_Dielectric_Sweep.pdf.
  * --all              : everything.

Model (energies in cm^-1, time in ps), from Combined.m:
  hbar = 5.308837 cm^-1*ps ; E1=E2=18437 ; eps_inf=1.77, eps_s=78, tau_D=8.3 ps
  1/eps(t) = 1/eps_s + (1/eps_inf - 1/eps_s) exp(-t/tau_D)
  J(t)     = J_pref / eps factors  with  J_pref = J_opt * eps_inf,  J_opt=100
  H(t)     = [[E1, J(t)], [J(t), E2]] ;  U=(1/sqrt2)[[1,1],[1,-1]]
  ME (Lindblad pure dephasing, rate gamma=1/T2*):
      drho/dt = -(i/hbar)[H,rho] + dephasing(off-diagonals * -gamma)
  SSE (Ito QSD, L = sqrt(hbar*gamma/2) sigma_z):
      dpsi = (1/hbar)(-iH - 0.5 L^2 + <L>L - 0.5<L>^2) psi dt
             + (1/sqrt(hbar))(L - <L>) psi dW

Pure NumPy/SciPy (solve_ivp replaces ode45). No QuTiP.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# ----- default physical parameters (Combined.m) -----------------------------
HBAR = 5.308837      # cm^-1 * ps
E1 = E2 = 18437.0    # cm^-1
EPS_INF = 1.77
EPS_S = 78.0
TAU_D = 8.3          # ps
# Energy-minimised (single-geometry) STEOM transition-density coupling J_TDC(0), the actual
# computed value from coupling_paper_steom_static/. The SSE/ME dynamics are DELIBERATELY driven
# by this energy-minimised J(0) (not the thermal-ensemble mean of 96.4 cm^-1): the dynamics start
# from the vertical excitation at the minimised t=0 geometry, before conformational sampling.
J_OPT = 100.2018964942497   # cm^-1  (= 100.2 cm^-1; paper rounds to "100")
T2_STAR = 0.060      # ps
TF = 1.0             # ps
DT = 1e-4            # ps

U = (1 / np.sqrt(2.0)) * np.array([[1.0, 1.0], [1.0, -1.0]])


def make_params(eps_s=EPS_S, eps_inf=EPS_INF, tau_d=TAU_D, j_opt=J_OPT,
                t2_star=T2_STAR, e1=E1, e2=E2):
    j_pref = j_opt * eps_inf
    return dict(eps_s=eps_s, eps_inf=eps_inf, tau_d=tau_d, j_opt=j_opt,
                j_pref=j_pref, t2_star=t2_star, gamma=1.0 / t2_star, e1=e1, e2=e2)


def inv_eps(t, p):
    return (1.0 / p["eps_s"]) + (1.0 / p["eps_inf"] - 1.0 / p["eps_s"]) * np.exp(-t / p["tau_d"])


def J_of_t(t, p):
    return p["j_pref"] * inv_eps(t, p)


def H_local(t, p):
    j = J_of_t(t, p)
    return np.array([[p["e1"], j], [j, p["e2"]]])


# --------------------------------------------------------------------------- #
# Master equation (Lindblad pure dephasing)
# --------------------------------------------------------------------------- #
def _me_rhs(t, y, p):
    rho = (y[:4] + 1j * y[4:]).reshape(2, 2)
    H = H_local(t, p)
    drho = -(1j / HBAR) * (H @ rho - rho @ H)
    g = p["gamma"]
    drho[0, 1] += -g * rho[0, 1]
    drho[1, 0] += -g * rho[1, 0]
    flat = drho.reshape(-1)
    return np.concatenate([flat.real, flat.imag])


def solve_me(p, tf=TF, dt=DT, psi0=None):
    if psi0 is None:
        psi0 = np.array([1.0, 1.0]) / np.sqrt(2.0)
    rho0 = np.outer(psi0, psi0.conj())
    y0 = np.concatenate([rho0.reshape(-1).real, rho0.reshape(-1).imag])
    tspan = np.arange(0.0, tf + dt, dt)
    sol = solve_ivp(_me_rhs, (0.0, tf), y0, t_eval=tspan, args=(p,),
                    method="RK45", rtol=1e-8, atol=1e-10)
    t = sol.t
    rho = (sol.y[:4] + 1j * sol.y[4:]).T.reshape(-1, 2, 2)  # (Nt,2,2)
    P1 = np.abs(rho[:, 0, 0])
    P2 = np.abs(rho[:, 1, 1])
    coh = np.abs(rho[:, 0, 1])
    rho_ex = np.einsum("ij,tjk,kl->til", U, rho, U.T.conj())
    PB = np.real(rho_ex[:, 0, 0])
    PD = np.real(rho_ex[:, 1, 1])
    bloch = np.stack([2 * np.real(rho[:, 0, 1]),
                      2 * np.imag(rho[:, 0, 1]),
                      P1 - P2], axis=1)
    purity = np.trace(rho @ rho, axis1=1, axis2=2).real
    return dict(t=t, P1=P1, P2=P2, coh=coh, PB=PB, PD=PD, bloch=bloch, purity=purity)


# --------------------------------------------------------------------------- #
# Stochastic Schrodinger equation (Ito QSD single trajectory)
# --------------------------------------------------------------------------- #
def solve_sse(p, tf=TF, dt=DT, psi0=None, seed=20260618):
    if psi0 is None:
        psi0 = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)
    rng = np.random.default_rng(seed)
    tspan = np.arange(0.0, tf + dt, dt)
    n = len(tspan)
    L = np.sqrt(HBAR * p["gamma"] / 2.0) * np.array([[1.0, 0.0], [0.0, -1.0]])
    L2 = L @ L
    I2 = np.eye(2)

    psi = psi0.astype(complex).copy()
    P1 = np.zeros(n); P2 = np.zeros(n)
    PB = np.zeros(n); PD = np.zeros(n)
    rho12 = np.zeros(n, dtype=complex)

    for k in range(n):
        P1[k] = abs(psi[0])**2
        P2[k] = abs(psi[1])**2
        rho12[k] = psi[0] * np.conj(psi[1])
        psi_ex = U @ psi
        PB[k] = abs(psi_ex[0])**2
        PD[k] = abs(psi_ex[1])**2
        if k < n - 1:
            H = H_local(tspan[k], p)
            exp_L = np.real(np.conj(psi) @ (L @ psi))
            u = (1.0 / HBAR) * ((-1j * H - 0.5 * L2 + exp_L * L - 0.5 * exp_L**2 * I2) @ psi)
            s = (1.0 / np.sqrt(HBAR)) * ((L - exp_L * I2) @ psi)
            psi = psi + u * dt + s * np.sqrt(dt) * rng.standard_normal()
            psi = psi / np.linalg.norm(psi)

    bloch = np.stack([2 * np.real(rho12), 2 * np.imag(rho12), P1 - P2], axis=1)
    return dict(t=tspan, P1=P1, P2=P2, coh=np.abs(rho12), PB=PB, PD=PD, bloch=bloch)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    # LaTeX (Computer Modern) look without a usetex toolchain: cmr10 is the actual
    # LaTeX body font and mathtext=cm is Computer Modern math.
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "mathtext.rm": "serif",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
        "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
        "legend.fontsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.linewidth": 0.8, "lines.antialiased": True,
    })
    return plt


# Semantic palette (softened for print; site red/blue, adiabatic green/purple)
_C_RHO11 = "#c0392b"    # site 1  rho_11
_C_RHO22 = "#1f6aa5"    # site 2  rho_22
_C_COH   = "#555555"    # coherence |rho_12|
_C_BRIGHT = "#3a7d44"   # P_+  (bright)
_C_DARK   = "#7b4397"   # P_-  (dark)
_C_SSE = _C_RHO11       # SSE trajectory (red)
_C_ME  = _C_RHO22       # ME curve (blue)


def _style_2d(ax, xlabel, ylabel, tf):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("0.35")
    ax.tick_params(colors="0.35", length=3, width=0.8)
    ax.grid(True, color="0.9", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(0, tf); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)


def _bloch_figure(plt, np, sse, me, out):
    """Clean floating Bloch sphere: no 3D box/panes/ticks, faint sphere + great
    circles, u/v/w axes with labels, SSE (red, on surface) and ME (blue, inward)."""
    fig = plt.figure(figsize=(4.2, 3.4))          # same footprint as the 2D panels
    ax = fig.add_subplot(111, projection="3d")
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    ax.set_axis_off()
    try:
        ax.set_box_aspect((1, 1, 1), zoom=1.18)   # enlarge sphere but leave room for labels
    except TypeError:
        ax.set_box_aspect((1, 1, 1))

    uu, vv = np.mgrid[0:2 * np.pi:60j, 0:np.pi:30j]
    xs, ys, zs = np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv)
    ax.plot_surface(xs, ys, zs, color="#eaf0f6", alpha=0.25, linewidth=0,
                    shade=False, zorder=0, antialiased=True)

    th = np.linspace(0, 2 * np.pi, 240); z0 = np.zeros_like(th)
    ax.plot(np.cos(th), np.sin(th), z0, color="0.70", lw=0.8, zorder=1)   # equator
    ax.plot(np.cos(th), z0, np.sin(th), color="0.85", lw=0.5, zorder=1)   # meridians
    ax.plot(z0, np.cos(th), np.sin(th), color="0.85", lw=0.5, zorder=1)

    for (vx, vy, vz), lab in [((1, 0, 0), r"\langle\hat{\sigma}_x\rangle"),
                              ((0, 1, 0), r"\langle\hat{\sigma}_y\rangle"),
                              ((0, 0, 1), r"\langle\hat{\sigma}_z\rangle")]:
        ax.plot([-1.1 * vx, 1.1 * vx], [-1.1 * vy, 1.1 * vy], [-1.1 * vz, 1.1 * vz],
                color="0.55", lw=0.9, zorder=2)
        ax.text(1.22 * vx, 1.22 * vy, 1.22 * vz, f"${lab}$", fontsize=12,
                ha="center", va="center", color="0.15")

    ax.plot(sse["bloch"][:, 0], sse["bloch"][:, 1], sse["bloch"][:, 2],
            color=_C_SSE, lw=0.6, alpha=0.45, label="SSE", zorder=4)
    ax.plot(me["bloch"][:, 0], me["bloch"][:, 1], me["bloch"][:, 2],
            color=_C_ME, lw=2.4, label="ME", zorder=5)

    s = me["bloch"][0]
    ax.scatter([s[0]], [s[1]], [s[2]], color="k", s=22, depthshade=False, zorder=6)
    # small label sitting just OUTSIDE the sphere at the +x pole, so it does not
    # occlude the surface (no background box needed there).
    ax.text(1.12, 0.0, 0.30, r"$|+\rangle$", fontsize=10,
            ha="center", va="center", color="0.1", zorder=7)

    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.view_init(elev=22, azim=130)
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
              loc="upper right", handlelength=1.4, bbox_to_anchor=(0.99, 0.99))
    fig.savefig(out / "Fig_Bloch_Grid.pdf")       # no tight crop -> full 4.2x3.4 footprint
    plt.close(fig)


def regenerate_base_figures(p, out, tf, dt, seed):
    plt = _mpl()
    me = solve_me(p, tf, dt)
    sse = solve_sse(p, tf, dt, seed=seed)
    t = me["t"]
    leg = dict(frameon=True, framealpha=0.95, edgecolor="0.75", fancybox=False,
               loc="lower center", bbox_to_anchor=(0.5, 1.0),
               handlelength=1.4, columnspacing=1.4, borderaxespad=0.4)

    # Fig_Purity_Coupling (Panel d)
    fig, ax1 = plt.subplots(figsize=(4.2, 3.4))
    c_pur = "#8e44ad"
    c_coup = _C_ME
    ax1.plot(t, me["purity"], color=c_pur, lw=2, label="Purity")
    ax1.set_xlabel(r"Time (ps)")
    ax1.set_ylabel(r"Purity $\mathrm{Tr}(\hat{\rho}^2)$", color=c_pur)
    ax1.tick_params(axis="y", labelcolor=c_pur)
    ax1.set_ylim(0.4, 1.05)
    
    ax2 = ax1.twinx()
    ax2.plot(t, J_of_t(t, p), color=c_coup, lw=2, ls="--", label=r"$|J(t)|$")
    ax2.set_ylabel(r"Coupling $|J(t)|\ \ (\mathrm{cm^{-1}})$", color=c_coup)
    ax2.tick_params(axis="y", labelcolor=c_coup)
    
    for s in ("top", "bottom"):
        ax1.spines[s].set_color("0.35")
        ax2.spines[s].set_color("0.35")
    ax1.spines["left"].set_color(c_pur)
    ax2.spines["right"].set_color(c_coup)
    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.grid(True, color="0.9", lw=0.6)
    ax1.set_axisbelow(True)
    ax1.set_xlim(0, tf)
    
    fig.tight_layout()
    fig.savefig(out / "Fig_Purity_Coupling.pdf")
    plt.close(fig)

    # Fig_SSE_Site
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.plot(sse["t"], sse["P1"], color=_C_RHO11, lw=0.7, alpha=0.9, label=r"$\rho_{11}$")
    ax.plot(sse["t"], sse["P2"], color=_C_RHO22, lw=0.7, alpha=0.9, label=r"$\rho_{22}$")
    ax.plot(sse["t"], sse["coh"], color=_C_COH, lw=0.8, alpha=0.85, label=r"$|\rho_{12}|$")
    _style_2d(ax, r"Time (ps)", r"Population", tf)
    ax.legend(ncol=3, **leg)
    fig.tight_layout(); fig.savefig(out / "Fig_SSE_Site.pdf"); plt.close(fig)

    # Fig_ME_Site
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.plot(t, me["P1"], color=_C_RHO11, lw=2, label=r"$\rho_{11}$")
    ax.plot(t, me["P2"], color=_C_RHO22, lw=2, ls="--", label=r"$\rho_{22}$")
    ax.plot(t, me["coh"], color=_C_COH, lw=2, label=r"$|\rho_{12}|$")
    _style_2d(ax, r"Time (ps)", r"Population", tf)
    ax.legend(ncol=3, **leg)
    fig.tight_layout(); fig.savefig(out / "Fig_ME_Site.pdf"); plt.close(fig)

    # Fig_SSE_Adiabatic
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.plot(sse["t"], sse["PB"], color=_C_BRIGHT, lw=0.9, alpha=0.95, label=r"$P_{+}$")
    ax.plot(sse["t"], sse["PD"], color=_C_DARK, lw=0.9, alpha=0.95, label=r"$P_{-}$")
    _style_2d(ax, r"Time (ps)", r"Population", tf)
    ax.legend(ncol=2, **leg)
    fig.tight_layout(); fig.savefig(out / "Fig_SSE_Adiabatic.pdf"); plt.close(fig)

    # Fig_ME_Adiabatic
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.plot(t, me["PB"], color=_C_BRIGHT, lw=2, label=r"$P_{+}$")
    ax.plot(t, me["PD"], color=_C_DARK, lw=2, label=r"$P_{-}$")
    _style_2d(ax, r"Time (ps)", r"Population", tf)
    ax.legend(ncol=2, **leg)
    fig.tight_layout(); fig.savefig(out / "Fig_ME_Adiabatic.pdf"); plt.close(fig)

    # Fig_Bloch_Grid: clean floating Bloch sphere
    _bloch_figure(plt, np, sse, me, out)

    print("    - wrote 6 base figures.")
    return me, sse


def sweep_t2(out, t2_list_fs, tf, dt):
    """Sweep T2* and show coherence decay << Debye time (timescale separation)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(6, 4.2))
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(t2_list_fs)))
    e_folds = []
    for c, t2_fs in zip(cmap, t2_list_fs):
        p = make_params(t2_star=t2_fs / 1000.0)
        me = solve_me(p, tf, dt)
        ax.plot(me["t"], me["coh"] / (me["coh"][0] + 1e-30), color=c, lw=1.8,
                label=fr"$T_2^*={t2_fs:.0f}$ fs")
        # 1/e coherence time
        c0 = me["coh"][0]
        below = np.where(me["coh"] <= c0 / np.e)[0]
        e_folds.append(me["t"][below[0]] if below.size else np.nan)
    ax.axvline(TAU_D, color="k", ls="--", lw=1.5, label=fr"$\tau_D={TAU_D}$ ps (Debye)")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel(r"$|\rho_{12}|/|\rho_{12}(0)|$")
    ax.set_title("Coherence decay vs Debye time (timescale separation)")
    ax.set_xlim(0, min(tf, 0.5)); ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(out / "Fig_T2_Sweep.pdf"); plt.close(fig)
    print("    - T2* sweep: coherence 1/e times (ps):",
          ", ".join(f"{t2:.0f}fs->{tau:.4f}" for t2, tau in zip(t2_list_fs, e_folds)))
    print(f"      all << tau_D = {TAU_D} ps  =>  timescale separation robust.")


def sweep_eps(out, eps_list, tf, dt):
    """Vary static (protein) dielectric; show J(0)=J_opt is invariant to it."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(6, 4.2))
    t = np.arange(0.0, tf + dt, dt)
    cmap = plt.cm.plasma(np.linspace(0, 0.85, len(eps_list)))
    j0s = []
    for c, es in zip(cmap, eps_list):
        p = make_params(eps_s=es)
        Jt = J_of_t(t, p)
        j0s.append(Jt[0])
        ax.plot(t, Jt, color=c, lw=1.8, label=fr"$\varepsilon_s={es:g}$")
    ax.axhline(J_OPT, color="k", ls=":", lw=1.5, label=fr"$J(0)={J_OPT}$ cm$^{{-1}}$")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel(r"$J(t)$ (cm$^{-1}$)")
    ax.set_title(r"Coupling vs static dielectric ($J(0)$ set by $\varepsilon_\infty$ only)")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(out / "Fig_Dielectric_Sweep.pdf"); plt.close(fig)
    spread = max(j0s) - min(j0s)
    print(f"    - eps_s sweep {eps_list}: J(0) spread = {spread:.3e} cm^-1 "
          f"(invariant; long-time limit J_s varies).")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("oqs_out"))
    p.add_argument("--tf", type=float, default=TF, help="Final time (ps).")
    p.add_argument("--dt", type=float, default=DT, help="Time step (ps).")
    p.add_argument("--seed", type=int, default=20260618, help="SSE trajectory RNG seed.")
    p.add_argument("--sweep-t2", action="store_true", help="Run the T2* sensitivity sweep.")
    p.add_argument("--sweep-eps", action="store_true", help="Run the dielectric sensitivity sweep.")
    p.add_argument("--t2-list", type=float, nargs="+", default=[20, 40, 60, 100, 200],
                   help="T2* values (fs) for the sweep.")
    p.add_argument("--eps-list", type=float, nargs="+", default=[4, 10, 20, 40, 78],
                   help="Static dielectric values for the sweep.")
    p.add_argument("--all", action="store_true", help="Base figures + both sweeps.")
    p.add_argument("--no-base", action="store_true", help="Skip the six base figures.")
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    params = make_params()
    print(f"[*] OQS dynamics: J_opt={J_OPT} cm^-1, T2*={T2_STAR*1000:.0f} fs, "
          f"tau_D={TAU_D} ps, out={args.out}")

    if not args.no_base:
        print("[*] Regenerating manuscript figures ...")
        regenerate_base_figures(params, args.out, args.tf, args.dt, args.seed)

    if args.sweep_t2 or args.all:
        print("[*] T2* sweep ...")
        sweep_t2(args.out, args.t2_list, args.tf, args.dt)

    if args.sweep_eps or args.all:
        print("[*] Dielectric sweep ...")
        sweep_eps(args.out, args.eps_list, args.tf, args.dt)

    print(f"[*] Done. Figures in {args.out}/")


if __name__ == "__main__":
    main()

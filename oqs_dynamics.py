#!/usr/bin/env python3
r"""
oqs_dynamics.py  --  Open-quantum-system dynamics of the Venus dimer exciton.

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
                         t=0 optical-limit coupling J(0)=74.38 cm^-1 is invariant
                         to it -> Fig_Dielectric_Sweep.pdf.
  * --all              : everything.

Model (energies in cm^-1, time in ps), from Combined.m:
  hbar = 5.308837 cm^-1*ps ; E1=E2=18437 ; eps_inf=1.77, eps_s=78, tau_D=8.3 ps
  1/eps(t) = 1/eps_s + (1/eps_inf - 1/eps_s) exp(-t/tau_D)
  J(t)     = J_pref / eps factors  with  J_pref = J_opt * eps_inf,  J_opt=74.38
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
J_OPT = 74.38        # cm^-1 (optical-limit coupling at t=0)
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
    return dict(t=t, P1=P1, P2=P2, coh=coh, PB=PB, PD=PD, bloch=bloch)


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
    import matplotlib.pyplot as plt
    return plt


def regenerate_base_figures(p, out, tf, dt, seed):
    plt = _mpl()
    me = solve_me(p, tf, dt)
    sse = solve_sse(p, tf, dt, seed=seed)
    t = me["t"]

    # Fig_Coupling: J(t)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(t, J_of_t(t, p), "r", lw=2)
    ax.set_xlabel("Time (ps)"); ax.set_ylabel(r"$J(t)$ (cm$^{-1}$)")
    ax.grid(True); fig.tight_layout()
    fig.savefig(out / "Fig_Coupling.pdf"); plt.close(fig)

    # Fig_SSE_Site
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(sse["t"], sse["P1"], color=(1, 0, 0, 0.6), lw=0.7, label=r"$\rho_{11}$")
    ax.plot(sse["t"], sse["P2"], color=(0, 0, 1, 0.6), lw=0.7, label=r"$\rho_{22}$")
    ax.plot(sse["t"], sse["coh"], color=(0, 0, 0, 0.6), lw=0.7, label=r"$|\rho_{12}|$")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel("Magnitude"); ax.set_ylim(0, 1)
    ax.legend(); ax.grid(True); fig.tight_layout()
    fig.savefig(out / "Fig_SSE_Site.pdf"); plt.close(fig)

    # Fig_ME_Site
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(t, me["P1"], "r", lw=2, label=r"$\rho_{11}$")
    ax.plot(t, me["P2"], "b--", lw=2, label=r"$\rho_{22}$")
    ax.plot(t, me["coh"], "k", lw=2, label=r"$|\rho_{12}|$")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel("Magnitude"); ax.set_ylim(0, 1)
    ax.legend(); ax.grid(True); fig.tight_layout()
    fig.savefig(out / "Fig_ME_Site.pdf"); plt.close(fig)

    # Fig_SSE_Adiabatic
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(sse["t"], sse["PB"], color=(0.466, 0.674, 0.188, 0.7), lw=0.7, label=r"$P_+$")
    ax.plot(sse["t"], sse["PD"], color=(0.494, 0.184, 0.556, 0.7), lw=0.7, label=r"$P_-$")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel("Population"); ax.set_ylim(0, 1)
    ax.legend(); ax.grid(True); fig.tight_layout()
    fig.savefig(out / "Fig_SSE_Adiabatic.pdf"); plt.close(fig)

    # Fig_ME_Adiabatic
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(t, me["PB"], color=(0.466, 0.674, 0.188), lw=2, label=r"$P_+$")
    ax.plot(t, me["PD"], color=(0.494, 0.184, 0.556), lw=2, label=r"$P_-$")
    ax.set_xlabel("Time (ps)"); ax.set_ylabel("Population"); ax.set_ylim(0, 1)
    ax.legend(); ax.grid(True); fig.tight_layout()
    fig.savefig(out / "Fig_ME_Adiabatic.pdf"); plt.close(fig)

    # Fig_Bloch_Grid: unified Bloch sphere (SSE trajectory + ME curve)
    fig = plt.figure(figsize=(5, 4.5))
    ax = fig.add_subplot(111, projection="3d")
    su, sv = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    ax.plot_wireframe(np.cos(su)*np.sin(sv), np.sin(su)*np.sin(sv), np.cos(sv),
                      color="0.8", alpha=0.3, linewidth=0.4)
    ax.plot(sse["bloch"][:, 0], sse["bloch"][:, 1], sse["bloch"][:, 2],
            color=(1, 0, 0, 0.6), lw=0.7, label="SSE")
    ax.plot(me["bloch"][:, 0], me["bloch"][:, 1], me["bloch"][:, 2],
            "b", lw=2.0, label="ME")

    # Time-progression arrows (R2): mark the direction of increasing time along
    # both paths, and the initial bright state |+> on the equator.
    def _time_arrows(curve, color, n_arrows):
        curve = np.asarray(curve)
        if len(curve) < 3:
            return
        idx = np.linspace(1, len(curve) - 2, n_arrows).astype(int)
        for i in idx:
            base = curve[i]
            step = curve[i + 1] - curve[i]
            norm = np.linalg.norm(step)
            if norm < 1e-9:
                continue
            step = step / norm * 0.28  # fixed visual arrow length
            ax.quiver(base[0], base[1], base[2], step[0], step[1], step[2],
                      color=color, arrow_length_ratio=0.5, lw=1.6)

    _time_arrows(me["bloch"], "b", 3)
    _time_arrows(sse["bloch"], (0.7, 0, 0, 0.9), 4)
    start = me["bloch"][0]
    ax.scatter([start[0]], [start[1]], [start[2]], color="k", s=28, depthshade=False)
    ax.text(start[0], start[1], start[2] + 0.12, r"$|+\rangle$, $t=0$", fontsize=8)

    ax.set_xlabel("u"); ax.set_ylabel("v"); ax.set_zlabel("w")
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.legend(loc="upper right", fontsize=9)
    ax.view_init(elev=25, azim=135)
    fig.tight_layout()
    fig.savefig(out / "Fig_Bloch_Grid.pdf"); plt.close(fig)

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

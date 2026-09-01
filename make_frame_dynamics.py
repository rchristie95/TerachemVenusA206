#!/usr/bin/env python3
"""Lindblad / stochastic-Schrodinger dynamics on the two extreme ensemble frames.

An ensemble-averaged mixing conceals how wide the detuning distribution is, so
rather than propagate a mean we take the two extreme members of the QM/MM
site-energy ensemble and propagate both:

  466 ps   Delta =    -0.17 cm^-1   near-degenerate, maximal delocalisation
  714 ps   Delta = -2495.52 cm^-1   maximal detuning, one chromophore favoured

Both are real frames of the v5 production trajectory (500 frames, def2-SVP
QM/MM with full-system embedding, transition densities recomputed per frame),
quoted by their time in the 1 ns production window.

H = [[Delta/2, J],[J, -Delta/2]] with each frame's OWN J rather than an ensemble
mean, which the per-frame densities now make available. Pure dephasing in the
SITE basis at T2* = 25 fs (measured; the campaign gives 20-31 fs across
configurations).

The Lindblad master equation gives the ensemble average; the diffusive
unraveling gives individual pure-state trajectories.

Superseded on 2026-08-31: this script previously used frames 465 and 1005 of the
1000-frame rigid-fit ensemble (Delta = +12.9 and -1825.0) with a single
J = 32.82 cm^-1. Both the ensemble and the coupling are replaced; see
tddft_ensemble_decoherence/V5_CAMPAIGN_RUNBOOK.md.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HBAR = 5.308837          # cm^-1 ps
T2 = 0.025               # ps, pure-dephasing time (v6t2 campaign, 2026-09-01)
# For L = sqrt(gamma/2) sigma_z the dissipator is (gamma/2)(sz rho sz - rho), and
# sigma_z flips the sign of an off-diagonal element, so d(rho_AB)/dt = -gamma
# rho_AB and the coherence decays as exp(-gamma t). The pure-dephasing rate is
# therefore 1/T2, not 2/T2.
#
# This was 2.0 / T2 until 2026-08-31, which delivered a 30 fs coherence time
# while the docstring and the adjacent panels from open_quantum_dynamics.py
# (gamma = 1/t2_star) both used 60 fs. Every panel this script produced dephased
# twice as fast as the figures beside it.
GAMMA = 1.0 / T2         # pure-dephasing rate (1/ps), sigma_z channel
# (label, Delta / cm^-1, J / cm^-1) -- J is now the frame's own recomputed value.
FRAMES = [("(a) 466 ps", -0.17, 24.55), ("(b) 714 ps", -2495.52, 24.65)]
SZ = np.array([[1, 0], [0, -1]], complex)

def hamiltonian(delta, J):
    return np.array([[delta / 2, J], [J, -delta / 2]], complex)

def lindblad(delta, J, t):
    H = hamiltonian(delta, J)
    rho = np.zeros((len(t), 2, 2), complex)
    # start in the bright eigenstate of this frame's Hamiltonian
    w, v = np.linalg.eigh(H)
    psi0 = v[:, np.argmax(np.abs(v.sum(axis=0)))]
    r = np.outer(psi0, psi0.conj())
    dt = t[1] - t[0]
    for i in range(len(t)):
        rho[i] = r
        drho = (-1j / HBAR) * (H @ r - r @ H) \
              + (GAMMA / 2) * (SZ @ r @ SZ - r)
        r = r + drho * dt
    return rho, psi0

def trajectory(delta, J, t, seed):
    """One diffusive (quantum-state-diffusion) pure-state path."""
    rng = np.random.default_rng(seed)
    H = hamiltonian(delta, J)
    w, v = np.linalg.eigh(H)
    psi = v[:, np.argmax(np.abs(v.sum(axis=0)))].astype(complex)
    dt = t[1] - t[0]
    out = np.zeros((len(t), 2), complex)
    L = np.sqrt(GAMMA / 2) * SZ
    for i in range(len(t)):
        out[i] = psi
        exp_L = (psi.conj() @ (L @ psi)).real
        dW = rng.normal(0.0, np.sqrt(dt))
        dpsi = (-1j / HBAR) * (H @ psi) * dt \
             + (exp_L * (L @ psi) - 0.5 * (L @ (L @ psi)) - 0.5 * exp_L**2 * psi) * dt \
             + (L @ psi - exp_L * psi) * dW
        psi = psi + dpsi
        psi = psi / np.linalg.norm(psi)
    return out

t = np.linspace(0, 0.5, 20001)      # ps
fig, axes = plt.subplots(2, 3, figsize=(11.5, 5.6))

for row, (label, delta, J) in enumerate(FRAMES):
    rho, psi0 = lindblad(delta, J, t)
    omega = np.hypot(delta, 2 * J)
    theta = 0.5 * np.arctan2(2 * J, delta)
    minor = min(np.sin(theta)**2, np.cos(theta)**2)

    # (i) site populations from the master equation
    ax = axes[row, 0]
    ax.plot(t * 1000, rho[:, 0, 0].real, color="tab:blue", lw=1.6, label="site A")
    ax.plot(t * 1000, rho[:, 1, 1].real, color="tab:red", lw=1.6, label="site B")
    ax.set_ylabel("site population")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{label}: $\\Delta$ = {delta:.2f} cm$^{{-1}}$, "
                 f"$J$ = {J:.2f}, $2|J|/\\Omega$ = {2*J/omega:.3f}", fontsize=9)
    if row == 0:
        ax.legend(fontsize=8, loc="center right")

    # (ii) coherence
    ax = axes[row, 1]
    ax.plot(t * 1000, np.abs(rho[:, 0, 1]), color="k", lw=1.6)
    ax.set_ylabel(r"$|\rho_{AB}|$")
    ax.set_ylim(-0.02, 0.55)
    ax.set_title(f"minor-site weight {100*minor:.2f}\\%", fontsize=9)

    # (iii) three stochastic trajectories, site-A population
    ax = axes[row, 2]
    for k, seed in enumerate((11, 12, 13)):
        psi = trajectory(delta, J, t, seed)
        ax.plot(t * 1000, np.abs(psi[:, 0])**2, lw=0.9, alpha=0.85)
    ax.plot(t * 1000, rho[:, 0, 0].real, color="k", lw=1.8, ls="--",
            label="ensemble")
    ax.set_ylabel("site-A population")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("individual trajectories", fontsize=9)
    if row == 0:
        ax.legend(fontsize=8, loc="center right")

for ax in axes[1, :]:
    ax.set_xlabel("time (fs)")
for ax in axes.ravel():
    ax.grid(alpha=0.25, lw=0.5)

fig.tight_layout()
fig.savefig("/home/robson/PetaChem/manuscript/Fig_FrameDynamics.pdf")
print("wrote Fig_FrameDynamics.pdf")
for label, delta, J in FRAMES:
    omega = np.hypot(delta, 2 * J)
    th = 0.5 * np.arctan2(2 * J, delta)
    print(f"  {label}: Delta={delta:9.2f}  J={J:.2f}  Omega={omega:7.1f}  "
          f"2J/Omega={2*J/omega:.4f}  minor weight={100*min(np.sin(th)**2,np.cos(th)**2):.2f}%")

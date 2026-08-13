#!/usr/bin/env python3
"""Lindblad / stochastic-Schrodinger dynamics on the two extreme ensemble frames.

An ensemble-averaged mixing conceals how wide the detuning distribution is, so
rather than propagate a mean we take the two extreme members of the QM/MM
site-energy ensemble and propagate both:

  frame 465   Delta =   +12.9 cm^-1   near-degenerate, maximal delocalisation
  frame 1005  Delta = -1825.0 cm^-1   maximal detuning, one chromophore favoured

Both are real frames of the production trajectory. H = [[Delta/2, J],[J, -Delta/2]]
with J = 32.82 cm^-1; pure dephasing in the SITE basis at T2* = 60 fs. The
Lindblad master equation gives the ensemble average; the diffusive unraveling
gives individual pure-state trajectories.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HBAR = 5.308837          # cm^-1 ps
J = 32.82                # cm^-1
T2 = 0.060               # ps
GAMMA = 2.0 / T2         # pure-dephasing rate (1/ps), sigma_z channel
FRAMES = [("(a) frame 465", 12.9), ("(b) frame 1005", -1825.0)]
SZ = np.array([[1, 0], [0, -1]], complex)

def hamiltonian(delta):
    return np.array([[delta / 2, J], [J, -delta / 2]], complex)

def lindblad(delta, t):
    H = hamiltonian(delta)
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

def trajectory(delta, t, seed):
    """One diffusive (quantum-state-diffusion) pure-state path."""
    rng = np.random.default_rng(seed)
    H = hamiltonian(delta)
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

for row, (label, delta) in enumerate(FRAMES):
    rho, psi0 = lindblad(delta, t)
    omega = np.hypot(delta, 2 * J)
    theta = 0.5 * np.arctan2(2 * J, delta)
    minor = min(np.sin(theta)**2, np.cos(theta)**2)

    # (i) site populations from the master equation
    ax = axes[row, 0]
    ax.plot(t * 1000, rho[:, 0, 0].real, color="tab:blue", lw=1.6, label="site A")
    ax.plot(t * 1000, rho[:, 1, 1].real, color="tab:red", lw=1.6, label="site B")
    ax.set_ylabel("site population")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{label}: $\\Delta$ = {delta:.1f} cm$^{{-1}}$, "
                 f"$2|J|/\\Omega$ = {2*J/omega:.3f}", fontsize=9)
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
        psi = trajectory(delta, t, seed)
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
for label, delta in FRAMES:
    omega = np.hypot(delta, 2 * J)
    th = 0.5 * np.arctan2(2 * J, delta)
    print(f"  {label}: Delta={delta:8.1f}  Omega={omega:7.1f}  "
          f"2J/Omega={2*J/omega:.4f}  minor weight={100*min(np.sin(th)**2,np.cos(th)**2):.2f}%")

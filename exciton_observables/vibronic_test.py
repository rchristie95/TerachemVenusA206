#!/usr/bin/env python3
"""Can a vibronic (Holstein) model escape the anisotropy/superradiance impossibility?

The purely electronic two-state model is dead: since cos^2 a <= |cos a|
pointwise, the anisotropy (<cos^2 a> = 0.436) and the superradiance
(<|cos a|> <= 0.33 at the measured detuning) cannot both hold for ANY geometry
or conformational ensemble. Vibronic structure is the one remaining mechanism
that changes BOTH relations rather than trading one against the other, because
it redistributes dipole strength across the vibronic manifold.

Two things are computed here on the REAL per-frame ensembles, at J = 30 cm^-1:

1. SUPERRADIANCE. The emitting manifold is the thermally populated set of
   vibronic eigenstates. Their dipole strengths need the interference term,
   |mu_k|^2 = |mu|^2 (a_k^2 + b_k^2 + 2 a_k b_k cos alpha), which
   `vibronic_exciton.frame_vibronic_observables` deliberately omits (it folds
   interference into the chirality factor instead), so it is added back here
   from the per-frame cos alpha.

2. CD APPARENT SPLITTING. The ensemble vibronic CD spectrum, fitted with
   Nguyen's constrained form (two couplet components pinned symmetrically at
   G1 -/+ delta), to see whether an apparent delta near their 130.79 cm^-1 can
   emerge from a true J of 30 cm^-1.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

REPO = Path("/home/robson/PetaChem")
sys.path.insert(0, str(REPO / "terachem_site_energy_cd"))
from vibronic_exciton import franck_condon_overlaps, vibronic_hamiltonian  # noqa: E402

KT = 208.509
J_CM = 30.0
OMEGA = 1450.0
S_HR = 0.35
N_MAX = 5
PHI, TAU_PS = 0.57, 3026.0
MEASURED_PS = 57.0


def load(name):
    d = np.load(REPO / f"terachem_site_energy_cd/results/ensembles/{name}")
    mu_a, mu_b = d["mu_a_au"], d["mu_b_au"]
    cos = np.sum(mu_a * mu_b, 1) / (np.linalg.norm(mu_a, axis=1) * np.linalg.norm(mu_b, axis=1))
    return d["e_a_cm"], d["e_b_cm"], cos, d["triple_product"]


def vibronic_frame(e_a, e_b, cos_alpha, j_cm=J_CM):
    """Eigenstates, dipole strengths (WITH interference), relative to monomer."""
    fc = franck_condon_overlaps(S_HR, N_MAX)
    H = vibronic_hamiltonian(e_a, e_b, j_cm, OMEGA, fc)
    energies, vectors = np.linalg.eigh(H)
    n = len(fc)
    a_k = fc @ vectors[:n, :]
    b_k = fc @ vectors[n:, :]
    # monomer reference: a single site carries sum_n f_n^2 = 1
    strength = a_k**2 + b_k**2 + 2.0 * a_k * b_k * cos_alpha
    chirality = np.einsum("n,nk,m,mk->k", fc, vectors[:n, :], fc, vectors[n:, :])
    return energies, strength, chirality


def main():
    print(f"Holstein vibronic model: J = {J_CM} cm^-1, omega = {OMEGA}, S = {S_HR}\n")
    fc = franck_condon_overlaps(S_HR, N_MAX)
    print(f"Franck-Condon: f_0^2 = {fc[0]**2:.4f}  -> the 0-0 inter-site coupling is")
    print(f"  J_eff = J f_0^2 = {J_CM*fc[0]**2:.2f} cm^-1, i.e. vibronic structure REDUCES")
    print(f"  the coupling between the two emitting states.\n")

    for name, label in (("ens_v2_all.npz", "crystal register (n=95)"),
                        ("ens_candidate2.npz", "candidate-2 register (n=23)")):
        e_a, e_b, cos, triple = load(name)
        dtau_vib, dtau_el = [], []
        for ea, eb, c in zip(e_a, e_b, cos):
            energies, strength, _ = vibronic_frame(ea, eb, c)
            # thermal population over the emitting manifold, referenced to its floor
            w = np.exp(-(energies - energies.min()) / KT)
            w /= w.sum()
            dimer = float(np.sum(w * strength))
            # Like-for-like monomer reference: the SAME vibronic manifold with
            # J = 0. Both dimer and monomer lose f_0^2 of their strength to the
            # sidebands, so that factor must cancel rather than appear as a
            # spurious 30% subradiance.
            e_m, s_m, _ = vibronic_frame(ea, eb, c, j_cm=0.0)
            w_m = np.exp(-(e_m - e_m.min()) / KT); w_m /= w_m.sum()
            monomer = float(np.sum(w_m * s_m))
            enh_vib = dimer / monomer - 1.0
            # electronic-only reference at the same frame
            delta = ea - eb
            omega = np.hypot(delta, 2 * J_CM)
            enh_el = -np.tanh(omega / (2 * KT)) * (2 * J_CM / omega) * c
            dtau_vib.append(PHI * enh_vib * TAU_PS)
            dtau_el.append(PHI * enh_el * TAU_PS)
        dtau_vib, dtau_el = np.array(dtau_vib), np.array(dtau_el)
        print(f"{label}:  <|cos a|> = {np.abs(cos).mean():.3f}, "
              f"<|Delta|> = {np.abs(e_a-e_b).mean():.0f} cm^-1")
        print(f"   predicted Dtau  electronic {dtau_el.mean():7.1f} ps   "
              f"vibronic {dtau_vib.mean():7.1f} ps   [measured {MEASURED_PS}]")
        print(f"   vibronic/electronic ratio = {dtau_vib.mean()/dtau_el.mean():.3f}\n")

    # ---- CD apparent splitting from the vibronic ensemble ----
    e_a, e_b, cos, triple = load("ens_v2_all.npz")
    centre = 0.5 * (e_a.mean() + e_b.mean())
    # Narrow window around the couplet only: the +/-1450 cm^-1 vibronic
    # sidebands are a separate feature (Nguyen fit those as a third band).
    grid = np.linspace(centre - 700, centre + 700, 4000)
    hwhm = 300.0
    cd = np.zeros_like(grid)
    for ea, eb, c, t in zip(e_a, e_b, cos, triple):
        energies, _, chir = vibronic_frame(ea, eb, c)
        rot = -np.pi * energies * chir * t
        for en, r in zip(energies, rot):
            cd += r * (hwhm / np.pi) / ((grid - en) ** 2 + hwhm**2)
    cd /= len(e_a)

    def couplet(x, g1, delta, a1, w1, a2, w2):
        return (a1 * (w1 / np.pi) / ((x - (g1 - delta))**2 + w1**2)
                + a2 * (w2 / np.pi) / ((x - (g1 + delta))**2 + w2**2))

    p0 = [centre, 150.0, cd.max() * 1e5, 300.0, -cd.max() * 1e5, 300.0]
    try:
        popt, _ = curve_fit(couplet, grid, cd, p0=p0, maxfev=40000)
        print(f"Nguyen-form constrained couplet fit of the VIBRONIC ensemble CD:")
        print(f"   apparent delta = {abs(popt[1]):.1f} cm^-1   "
              f"(apparent splitting 2delta = {2*abs(popt[1]):.1f})")
        print(f"   true 2|J|      = {2*J_CM:.1f} cm^-1")
        print(f"   Nguyen report  delta = 130.79 (2delta = 261.58)")
        print(f"   inflation over true J: {abs(popt[1])/J_CM:.2f}x")
    except Exception as exc:
        print(f"   couplet fit failed: {exc}")

    print("\nVERDICT on the impossibility:")
    c2_req = (2 * (2 * 0.30 / 0.52 - 1) + 1) / 3
    print(f"   anisotropy still requires <cos^2 a> = {c2_req:.3f}; vibronic structure")
    print(f"   does not change that relation, and it REDUCES the superradiance, so the")
    print(f"   bound <cos^2 a> <= <|cos a|> is violated by MORE, not less.")


if __name__ == "__main__":
    main()

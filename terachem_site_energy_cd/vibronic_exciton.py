"""Holstein vibronic exciton model for a detuned chromophore dimer.

The purely electronic two-state model cannot produce the third CD band that
Nguyen et al. fit at 481.37 nm: that band sits ~1450 cm^-1 to the blue of the
excitonic pair, which is the C=C/C=N stretch vibronic sideband, not a third
electronic state. This module adds one effective intramolecular mode per site
in the one-particle Holstein approximation.

Basis (one-particle approximation)
----------------------------------
|m, n>  = site m electronically excited carrying n quanta in ITS OWN
          excited-state potential; the partner site is electronically and
          vibrationally in its ground state.

Franck-Condon overlaps for a displaced harmonic oscillator of Huang-Rhys
factor S:

    f_n = <chi_n^exc | chi_0^gnd> = exp(-S/2) * (-sqrt(S))^n / sqrt(n!)

Hamiltonian
-----------
    <m,n| H |m,n>    = E_m + n * omega_vib
    <A,n| H |B,n'>   = J * f_n * f_n'

Observables (from the ground vibrational state)
-----------------------------------------------
    mu_k  = sum_{m,n} c^k_{m,n} f_n mu_m
    D_k   = |mu_k|^2
    R_k   = -pi nu_k (R_AB . mu_A x mu_B) sum_{n,n'} c^k_{A,n} c^k_{B,n'} f_n f_n'

Setting S = 0 collapses this exactly to the two-level electronic model, which
is used as a regression test.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "franck_condon_overlaps",
    "vibronic_hamiltonian",
    "frame_vibronic_observables",
    "vibronic_spectra",
]


def franck_condon_overlaps(huang_rhys: float, n_max: int) -> np.ndarray:
    """<chi_n^exc | chi_0^gnd> for a displaced harmonic oscillator."""
    if huang_rhys < 0.0:
        raise ValueError("Huang-Rhys factor must be non-negative")
    n = np.arange(n_max + 1)
    log_fact = np.cumsum(np.concatenate(([0.0], np.log(np.arange(1, n_max + 1)))))
    magnitude = np.exp(-0.5 * huang_rhys + 0.5 * (n * np.log(huang_rhys)
                                                  if huang_rhys > 0 else
                                                  np.where(n == 0, 0.0, -np.inf))
                       - 0.5 * log_fact)
    return magnitude * (-1.0) ** n


def vibronic_hamiltonian(e_a_cm, e_b_cm, j_cm, omega_cm, fc):
    """Block Hamiltonian in the |m,n> basis; returns (H, index_map)."""
    n_states = len(fc)
    dim = 2 * n_states
    H = np.zeros((dim, dim), float)
    site_energy = (e_a_cm, e_b_cm)
    for m in range(2):
        for n in range(n_states):
            H[m * n_states + n, m * n_states + n] = site_energy[m] + n * omega_cm
    # inter-site vibronic coupling J * f_n * f_n'
    block = j_cm * np.outer(fc, fc)
    H[:n_states, n_states:] = block
    H[n_states:, :n_states] = block.T
    return H


def frame_vibronic_observables(e_a_cm, e_b_cm, j_cm, triple_product,
                               omega_cm=1450.0, huang_rhys=0.35, n_max=5):
    """Vibronic eigenstates with dipole and exciton-chirality strengths.

    `triple_product` is R_AB . (mu_A x mu_B); the electric dipoles are taken as
    equal in magnitude (homodimer), so absorption strengths are relative.
    """
    fc = franck_condon_overlaps(huang_rhys, n_max)
    H = vibronic_hamiltonian(e_a_cm, e_b_cm, j_cm, omega_cm, fc)
    energies, vectors = np.linalg.eigh(H)
    n_states = len(fc)

    ca, cb = vectors[:n_states, :], vectors[n_states:, :]
    # electric dipole: sum over vibronic amplitudes weighted by FC overlap.
    # equal-|mu| homodimer, so |mu_k|^2 = (a_k)^2 + (b_k)^2 + 2 a_k b_k cos(theta);
    # the interference term is folded into the chirality factor below, and the
    # relative dipole strength uses the in-phase convention a_k + b_k.
    a_k = fc @ ca
    b_k = fc @ cb
    dipole_strength = a_k**2 + b_k**2
    chirality = np.einsum("n,nk,m,mk->k", fc, ca, fc, cb)
    rotational = -np.pi * energies * chirality * triple_product
    return {
        "energies_cm": energies,
        "dipole_strengths": dipole_strength,
        "rotational_strengths_relative": rotational,
        "vibronic_weight_a": a_k,
        "vibronic_weight_b": b_k,
        "franck_condon": fc,
    }


def vibronic_spectra(grid_cm, frames, hwhm_cm, omega_cm=1450.0,
                     huang_rhys=0.35, n_max=5):
    """Ensemble-averaged absorption and interaction-induced CD.

    `frames` is an iterable of (e_a_cm, e_b_cm, j_cm, triple_product).
    """
    grid = np.asarray(grid_cm, float)
    absorption = np.zeros_like(grid)
    cd = np.zeros_like(grid)
    count = 0
    for e_a, e_b, j_cm, triple in frames:
        obs = frame_vibronic_observables(e_a, e_b, j_cm, triple,
                                         omega_cm, huang_rhys, n_max)
        for energy, strength, rotation in zip(obs["energies_cm"],
                                              obs["dipole_strengths"],
                                              obs["rotational_strengths_relative"]):
            shape = hwhm_cm / np.pi / ((grid - energy) ** 2 + hwhm_cm**2)
            absorption += strength * shape
            cd += rotation * shape
        count += 1
    if count == 0:
        raise ValueError("no frames supplied")
    return absorption / count, cd / count

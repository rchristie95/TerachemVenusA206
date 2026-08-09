"""Nondegenerate two-site absorption and interaction-induced exciton CD."""

from __future__ import annotations

import numpy as np


def diagonalize_site_hamiltonian(e1_cm: float, e2_cm: float, j_cm: float):
    """Return ascending energies and site coefficients (columns are eigenstates)."""
    if j_cm == 0.0:
        if e1_cm <= e2_cm:
            return np.array([e1_cm, e2_cm], float), np.eye(2)
        return np.array([e2_cm, e1_cm], float), np.array([[0.0, 1.0], [1.0, 0.0]])
    return np.linalg.eigh(np.array([[e1_cm, j_cm], [j_cm, e2_cm]], float))


def frame_observables(e1_cm, e2_cm, j_cm, mu1, mu2, r1_A, r2_A):
    """Return eigenvectors, absorption strengths, and relative exciton-chirality strengths.

    Dipoles use one common Cartesian frame. Positions are Angstrom and energies
    are cm^-1. Rotational strengths omit the common unit-conversion prefactor,
    so they are interaction-induced relative signals, not absolute molar CD.
    """
    energies, coeff = diagonalize_site_hamiltonian(e1_cm, e2_cm, j_cm)
    mus = np.asarray([mu1, mu2], float)
    exciton_mu = coeff.T @ mus
    dipole_strength = np.einsum("ij,ij->i", exciton_mu, exciton_mu)
    triple = float(np.dot(np.asarray(r2_A) - np.asarray(r1_A), np.cross(mu1, mu2)))
    rotational = -np.pi * energies * coeff[0] * coeff[1] * triple
    omega = float(np.hypot(e1_cm - e2_cm, 2.0 * j_cm))
    mixing = 0.0 if omega == 0.0 else 2.0 * abs(j_cm) / omega
    return {
        "energies_cm": energies, "coefficients": coeff, "transition_dipoles": exciton_mu,
        "dipole_strengths": dipole_strength, "rotational_strengths_relative": rotational,
        "omega_cm": omega, "mixing": mixing, "localization_weights": coeff * coeff,
    }


def lorentzian(grid, center, hwhm):
    return hwhm / np.pi / ((np.asarray(grid) - center) ** 2 + hwhm**2)


def spectra(grid_cm, frames, hwhm_cm):
    absorption = np.zeros_like(np.asarray(grid_cm, float))
    cd = np.zeros_like(absorption)
    for frame in frames:
        obs = frame_observables(**frame)
        for energy, strength, rotation in zip(obs["energies_cm"], obs["dipole_strengths"], obs["rotational_strengths_relative"]):
            shape = lorentzian(grid_cm, energy, hwhm_cm)
            absorption += strength * shape
            cd += rotation * shape
    return absorption / len(frames), cd / len(frames)

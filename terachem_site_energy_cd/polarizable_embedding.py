"""Perturbative polarizable embedding (QM/MMpol) correction to site energies.

Our QM/MM lets the QM density relax in the static MM field, but the MM
environment never polarises in response. That missing back-reaction is the
electronic screening a fixed-charge force field cannot provide, and it is the
remaining candidate for the site-energy disorder being too large: rebuilding the
MD (neutral box, NPT density, seeded velocities) left the disorder statistically
unchanged (Mann-Whitney p = 0.45), so the defect is in the force-field model
rather than the setup.

Model
-----
Each MM atom b carries an isotropic polarizability alpha_b. The chromophore's
charge distribution creates a field F_b there; on excitation that field changes
by dF_b, driven by the STEOM excited-minus-ground atomic charges dq.

Polarization energy of a state is E_pol = -1/2 sum_b alpha_b |F_b|^2, so the
contribution to the EXCITATION energy is

    dE_pol = -sum_b alpha_b F_b^g . dF_b  -  1/2 sum_b alpha_b |dF_b|^2

The first term is the cross term between the ground-state field and its change;
the second is the self term. Both are negative-definite in the usual case, i.e.
polarization stabilises the excited state and REDUCES the magnitude of
environment-induced site-energy shifts - which is the screening we are missing.

This is first-order (non-iterative): induced dipoles do not mutually polarise.
That captures the leading effect and is a strict bound on it, since mutual
polarization in a Thole-damped treatment further reduces the net response.

Limitations, stated plainly: isotropic atom-centred polarizabilities, no Thole
damping at short range, no iteration to self-consistency, and the ground-state
chromophore charge distribution is taken from its MM (RESP) charges.
"""

from __future__ import annotations

import numpy as np

__all__ = ["ATOMIC_POLARIZABILITY_ANG3", "polarization_shift_cm", "BOHR_A"]

BOHR_A = 0.529177210903
HARTREE_CM = 219474.6314

# Isotropic atomic polarizabilities, Angstrom^3 (Thole/AMOEBA-style values).
ATOMIC_POLARIZABILITY_ANG3 = {
    "H": 0.496, "C": 1.334, "N": 1.073, "O": 0.837,
    "S": 2.926, "P": 1.828, "F": 0.507, "CL": 2.500,
    "NA": 0.120, "K": 0.780, "MG": 0.080, "ZN": 0.260,
}
DEFAULT_POLARIZABILITY_ANG3 = 1.0


def _polarizabilities(symbols) -> np.ndarray:
    return np.array([ATOMIC_POLARIZABILITY_ANG3.get(str(s).upper(),
                                                    DEFAULT_POLARIZABILITY_ANG3)
                     for s in symbols], float)


def _field(source_xyz_ang, source_q, target_xyz_ang, chunk=4096) -> np.ndarray:
    """Electrostatic field (atomic units) at each target from point charges."""
    src = np.asarray(source_xyz_ang, float) / BOHR_A
    tgt = np.asarray(target_xyz_ang, float) / BOHR_A
    q = np.asarray(source_q, float)
    out = np.zeros_like(tgt)
    for start in range(0, len(tgt), chunk):
        block = tgt[start:start + chunk]
        d = block[:, None, :] - src[None, :, :]
        r2 = np.einsum("ijk,ijk->ij", d, d)
        np.maximum(r2, 1.0e-6, out=r2)
        inv_r3 = r2 ** -1.5
        out[start:start + chunk] = np.einsum("ij,j,ijk->ik", inv_r3, q, d)
    return out


def polarization_shift_cm(qm_xyz_ang, qm_ground_q, qm_delta_q,
                          mm_xyz_ang, mm_symbols, polarizability_scale=1.0):
    """Return (dE_pol, cross_term, self_term) in cm^-1 for one site.

    A negative dE_pol lowers the excitation energy. What matters for the dimer is
    the DIFFERENCE between the two sites, since a common shift cancels in the
    detuning.
    """
    alpha = _polarizabilities(mm_symbols) * polarizability_scale / BOHR_A**3
    f_ground = _field(qm_xyz_ang, qm_ground_q, mm_xyz_ang)
    f_delta = _field(qm_xyz_ang, qm_delta_q, mm_xyz_ang)
    cross = -float(np.sum(alpha * np.einsum("ij,ij->i", f_ground, f_delta)))
    self_term = -0.5 * float(np.sum(alpha * np.einsum("ij,ij->i", f_delta, f_delta)))
    return ((cross + self_term) * HARTREE_CM,
            cross * HARTREE_CM,
            self_term * HARTREE_CM)

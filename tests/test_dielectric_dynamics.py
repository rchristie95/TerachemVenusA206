#!/usr/bin/env python3
"""Regression check for the first-picosecond dielectric-screening sensitivity."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import open_quantum_dynamics as oqd


def main() -> None:
    time_dependent = oqd.make_params()
    # eps_s == eps_inf makes inv_eps(t) constant, hence J(t) == J(0).
    frozen = oqd.make_params(eps_s=oqd.EPS_INF)

    inv_0 = float(oqd.inv_eps(0.0, time_dependent))
    inv_1 = float(oqd.inv_eps(1.0, time_dependent))
    drop_percent = 100.0 * (
        1.0 - oqd.J_of_t(1.0, time_dependent) / oqd.J_of_t(0.0, time_dependent)
    )
    assert np.isclose(inv_0, 0.5649717514124294, rtol=0.0, atol=1.0e-15)
    assert np.isclose(inv_1, 0.5022987804944742, rtol=0.0, atol=1.0e-15)
    assert np.isclose(drop_percent, 11.093115852478075, rtol=0.0, atol=1.0e-12)

    # A LOCALISED initial state is required for this comparison to mean
    # anything. With the default symmetric psi0 = (1,1)/sqrt(2) and E1 == E2,
    # rho0 commutes with H(t) for every value of J, so -(i/hbar)[H,rho]
    # vanishes identically at every timestep and the evolving and frozen runs
    # are bit-identical no matter what J_of_t returns -- the old
    # `max(differences) < 1e-12` assertion passed even with J_of_t
    # deliberately corrupted, providing zero regression protection.
    localised = np.array([1.0, 0.0])
    evolving_result = oqd.solve_me(time_dependent, tf=1.0, dt=1.0e-3, psi0=localised)
    frozen_result = oqd.solve_me(frozen, tf=1.0, dt=1.0e-3, psi0=localised)
    keys = ("P1", "P2", "coh", "PB", "PD", "bloch", "purity")
    differences = {
        key: float(np.max(np.abs(evolving_result[key] - frozen_result[key])))
        for key in keys
    }
    # J drops 11.1% over the first ps, so a site-localised population MUST
    # transfer differently under evolving and frozen screening. Observed
    # separation is ~1.3e-2 in the Bloch vector and ~6.5e-3 in P1.
    assert max(differences.values()) > 1.0e-3, differences
    assert differences["P1"] > 1.0e-3, differences
    print(
        f"PASS: 1/eps 0->1 ps = {inv_0:.10f}->{inv_1:.10f}; "
        f"J drop = {drop_percent:.6f}%; max ME difference = "
        f"{max(differences.values()):.3e}"
    )


def test_dielectric_dynamics() -> None:
    main()


if __name__ == "__main__":
    main()

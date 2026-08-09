import numpy as np

from terachem_site_energy_cd.nondegenerate_spectra import frame_observables


BASE = dict(mu1=np.array([1.0, 0.2, 0.0]), mu2=np.array([0.1, 1.0, 0.3]),
            r1_A=np.array([0.0, 0.0, 0.0]), r2_A=np.array([0.0, 0.0, 20.0]))


def test_degenerate_limit():
    result = frame_observables(19000.0, 19000.0, 33.0, **BASE)
    assert np.allclose(result["energies_cm"], [18967.0, 19033.0])
    assert np.isclose(result["mixing"], 1.0)


def test_zero_coupling_is_localized_and_has_no_interaction_cd():
    result = frame_observables(18800.0, 19200.0, 0.0, **BASE)
    assert np.allclose(result["localization_weights"], np.eye(2))
    assert np.allclose(result["rotational_strengths_relative"], 0.0)


def test_site_swap_invariance():
    a = frame_observables(18800.0, 19200.0, 33.0, **BASE)
    swapped = dict(mu1=BASE["mu2"], mu2=BASE["mu1"], r1_A=BASE["r2_A"], r2_A=BASE["r1_A"])
    b = frame_observables(19200.0, 18800.0, 33.0, **swapped)
    for key in ("energies_cm", "dipole_strengths", "rotational_strengths_relative"):
        assert np.allclose(a[key], b[key])


def test_handedness_reversal_flips_cd_only():
    a = frame_observables(18800.0, 19200.0, 33.0, **BASE)
    mirrored = {key: np.array(value) * np.array([1.0, 1.0, -1.0]) for key, value in BASE.items() if key.startswith("mu") or key.startswith("r")}
    b = frame_observables(18800.0, 19200.0, 33.0, **mirrored)
    assert np.allclose(a["dipole_strengths"], b["dipole_strengths"])
    assert np.allclose(a["rotational_strengths_relative"], -b["rotational_strengths_relative"])


def test_large_detuning_reduces_mixing_formula():
    result = frame_observables(18000.0, 19000.0, 25.0, **BASE)
    assert np.isclose(result["mixing"], 50.0 / np.sqrt(1000.0**2 + 50.0**2))

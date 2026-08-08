"""Regression tests for rejecting pre-correction coupling caches."""

import copy

import pytest

import reproduce_paper


NEW_METADATA = {
    "tdc_units": {
        "status": "corrected",
        "pair_distance_unit": "angstrom",
        "reciprocal_distance_to_atomic_units": 0.529177210903,
    }
}

MIGRATED_METADATA = {
    "unit_correction": {
        "status": "corrected",
        "new_reciprocal_distance_factor": 0.529177210903,
    }
}


@pytest.mark.parametrize("metadata", [NEW_METADATA, MIGRATED_METADATA])
def test_corrected_metadata_is_accepted(metadata):
    assert reproduce_paper._has_corrected_tdc_units(copy.deepcopy(metadata))


@pytest.mark.parametrize(
    "metadata",
    [
        {},
        {"tdc_units": {"status": "corrected"}},
        {
            "tdc_units": {
                "status": "corrected",
                "pair_distance_unit": "angstrom",
                "reciprocal_distance_to_atomic_units": 1.889726124626,
            }
        },
        {
            "unit_correction": {
                "status": "corrected",
                "new_reciprocal_distance_factor": 1.889726124626,
            }
        },
    ],
)
def test_missing_or_old_conversion_is_rejected(metadata):
    assert not reproduce_paper._has_corrected_tdc_units(metadata)


def test_coupling_summary_rejects_unmarked_cache(tmp_path, monkeypatch):
    out = tmp_path / "cache"
    out.mkdir()
    (out / "coupling_distribution.json").write_text(
        '{"mean": 117.2, "std": 5.5, "n": 1000, "epsilon": 1.77}'
    )
    monkeypatch.setattr(reproduce_paper, "REPO", tmp_path)

    with pytest.raises(ValueError, match="Refusing unverified coupling cache"):
        reproduce_paper.coupling_summary("cache")

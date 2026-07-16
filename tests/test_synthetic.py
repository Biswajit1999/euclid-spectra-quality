from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.masks import CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED
from euclid_nisp_contamination_audit.synthetic import build_synthetic_catalogue, build_synthetic_spectrum


def test_build_synthetic_catalogue_deterministic() -> None:
    a = build_synthetic_catalogue(n_sources=50, seed=1)
    b = build_synthetic_catalogue(n_sources=50, seed=1)
    assert np.array_equal(a.spe_cont_snr, b.spe_cont_snr)


def test_build_synthetic_catalogue_rejects_small_n() -> None:
    with pytest.raises(ValueError):
        build_synthetic_catalogue(n_sources=2)


def test_build_synthetic_catalogue_has_both_groups() -> None:
    cat = build_synthetic_catalogue(n_sources=200, seed=7)
    assert CONTAMINATION_GROUP_CLEAN in cat.contamination_group
    assert CONTAMINATION_GROUP_FLAGGED in cat.contamination_group


def test_build_synthetic_spectrum_shapes_match() -> None:
    spec = build_synthetic_spectrum(object_id=1, seed=1, contaminated=False)
    n = spec.wavelength_angstrom.size
    assert spec.signal.shape == (n,)
    assert spec.uncertainty.shape == (n,)
    assert spec.mask.shape == (n,)
    assert spec.quality.shape == (n,)


def test_contaminated_spectrum_has_higher_noise() -> None:
    clean = build_synthetic_spectrum(object_id=1, seed=1, contaminated=False)
    flagged = build_synthetic_spectrum(object_id=2, seed=1, contaminated=True)
    assert np.mean(flagged.uncertainty) > np.mean(clean.uncertainty)

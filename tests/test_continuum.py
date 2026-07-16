from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.continuum import estimate_continuum
from euclid_nisp_contamination_audit.exceptions import InsufficientDataError


def test_estimate_continuum_recovers_flat_level() -> None:
    rng = np.random.default_rng(0)
    wavelength = np.linspace(12000, 18000, 200)
    signal = np.full(200, 3.0) + rng.normal(0, 0.1, size=200)
    uncertainty = np.full(200, 0.1)
    good = np.ones(200, dtype=bool)

    result = estimate_continuum(wavelength, signal, uncertainty, good)
    assert np.allclose(result.continuum, 3.0, atol=0.5)
    assert result.continuum_snr > 5
    assert result.n_pixels == 200


def test_estimate_continuum_insufficient_pixels() -> None:
    wavelength = np.linspace(12000, 13000, 5)
    signal = np.ones(5)
    uncertainty = np.ones(5) * 0.1
    good = np.ones(5, dtype=bool)
    with pytest.raises(InsufficientDataError):
        estimate_continuum(wavelength, signal, uncertainty, good)


def test_estimate_continuum_ignores_bad_pixels() -> None:
    wavelength = np.linspace(12000, 18000, 100)
    signal = np.full(100, 3.0)
    signal[:20] = 1000.0  # gross outliers, masked out
    uncertainty = np.full(100, 0.1)
    good = np.ones(100, dtype=bool)
    good[:20] = False

    result = estimate_continuum(wavelength, signal, uncertainty, good)
    assert np.median(result.continuum) < 10

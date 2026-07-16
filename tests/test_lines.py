from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.exceptions import InsufficientDataError
from euclid_nisp_contamination_audit.lines import (
    inject_gaussian_line,
    line_recovery_rate,
    recover_injected_line,
)


def test_inject_and_recover_strong_line() -> None:
    rng = np.random.default_rng(42)
    wavelength = np.linspace(15000, 16000, 500)
    continuum = np.full(500, 2.0)
    uncertainty = np.full(500, 0.05)
    center = 15500.0
    sigma = 8.0
    amplitude = 5.0
    injected_flux = amplitude * sigma * np.sqrt(2 * np.pi)

    signal = inject_gaussian_line(wavelength, continuum, uncertainty, center, amplitude, sigma, rng)
    result = recover_injected_line(wavelength, signal, uncertainty, center, window_angstrom=100.0, injected_flux=injected_flux)

    assert result.detected
    assert 0.7 < result.flux_ratio < 1.3


def test_recover_injected_line_insufficient_window() -> None:
    wavelength = np.linspace(15000, 15010, 3)
    signal = np.array([1.0, 2.0, 1.0])
    uncertainty = np.array([0.1, 0.1, 0.1])
    with pytest.raises(InsufficientDataError):
        recover_injected_line(wavelength, signal, uncertainty, 15005.0, window_angstrom=5.0, injected_flux=1.0)


def test_line_recovery_rate() -> None:
    snr = np.array([10.0, 2.0, 6.0, np.nan, 4.9])
    rate = line_recovery_rate(snr, snr_threshold=5.0)
    assert rate == pytest.approx(2 / 4)


def test_line_recovery_rate_no_finite_values() -> None:
    with pytest.raises(InsufficientDataError):
        line_recovery_rate(np.array([np.nan, np.nan]))

"""Robust local continuum estimation and continuum S/N for NISP 1-D spectra."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from euclid_nisp_contamination_audit.exceptions import InsufficientDataError

_MIN_GOOD_PIXELS = 10


@dataclass(frozen=True)
class ContinuumEstimate:
    continuum: np.ndarray
    residual: np.ndarray
    continuum_snr: float
    n_pixels: int


def estimate_continuum(
    wavelength: np.ndarray,
    signal: np.ndarray,
    uncertainty: np.ndarray,
    good_pixel_mask: np.ndarray,
    poly_degree: int = 2,
) -> ContinuumEstimate:
    """Fit a low-order polynomial continuum to the good pixels of one spectrum
    and compute a continuum S/N as median(continuum / local rms residual).

    Requires at least `_MIN_GOOD_PIXELS` good, finite pixels; raises
    InsufficientDataError otherwise rather than returning a degenerate fit.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    signal = np.asarray(signal, dtype=float)
    uncertainty = np.asarray(uncertainty, dtype=float)
    good_pixel_mask = np.asarray(good_pixel_mask, dtype=bool)

    finite = np.isfinite(wavelength) & np.isfinite(signal) & np.isfinite(uncertainty)
    usable = good_pixel_mask & finite & (uncertainty > 0)
    n_good = int(np.sum(usable))
    if n_good < _MIN_GOOD_PIXELS:
        raise InsufficientDataError(
            f"only {n_good} usable pixels (< {_MIN_GOOD_PIXELS} required) for continuum estimation"
        )

    degree = min(poly_degree, max(0, n_good - 2))
    coeffs = np.polyfit(wavelength[usable], signal[usable], deg=degree)
    continuum_full = np.polyval(coeffs, wavelength)
    residual_full = signal - continuum_full

    residual_good = residual_full[usable]
    rms = float(np.std(residual_good))
    continuum_level = float(np.median(np.abs(continuum_full[usable])))
    continuum_snr = continuum_level / rms if rms > 0 else float("inf")

    return ContinuumEstimate(
        continuum=continuum_full,
        residual=residual_full,
        continuum_snr=continuum_snr,
        n_pixels=n_good,
    )


__all__ = ["ContinuumEstimate", "estimate_continuum"]

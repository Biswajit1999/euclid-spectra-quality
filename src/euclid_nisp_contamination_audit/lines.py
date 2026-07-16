"""Synthetic emission-line injection/recovery, and a wrapper around the
released per-object best line S/N used as the real-catalogue line-recovery
indicator.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from euclid_nisp_contamination_audit.exceptions import ConvergenceError, InsufficientDataError

DEFAULT_LINE_RECOVERY_SNR_THRESHOLD = 5.0


def _gaussian(x: np.ndarray, amplitude: float, center: float, sigma: float, offset: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2) + offset


@dataclass(frozen=True)
class LineInjectionResult:
    injected_flux: float
    recovered_flux: float
    recovered_flux_err: float
    recovered_snr: float
    detected: bool
    flux_ratio: float


def inject_gaussian_line(
    wavelength: np.ndarray,
    continuum_signal: np.ndarray,
    uncertainty: np.ndarray,
    center_angstrom: float,
    amplitude: float,
    sigma_angstrom: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return a new signal array with a Gaussian emission line of known
    `amplitude` (same flux units as `continuum_signal`) added at
    `center_angstrom`, plus one new Gaussian-noise realisation drawn from
    `uncertainty` (so repeated injections are independent noise draws on the
    same real continuum, not the same noisy pixels re-used).
    """
    wavelength = np.asarray(wavelength, dtype=float)
    continuum_signal = np.asarray(continuum_signal, dtype=float)
    uncertainty = np.asarray(uncertainty, dtype=float)
    if wavelength.shape != continuum_signal.shape or wavelength.shape != uncertainty.shape:
        raise InsufficientDataError("wavelength, continuum_signal and uncertainty must share a shape")

    line = _gaussian(wavelength, amplitude, center_angstrom, sigma_angstrom, offset=0.0)
    noise = rng.normal(loc=0.0, scale=np.clip(uncertainty, a_min=1e-6, a_max=None))
    return continuum_signal + line + noise


def recover_injected_line(
    wavelength: np.ndarray,
    signal_with_line: np.ndarray,
    uncertainty: np.ndarray,
    expected_center_angstrom: float,
    window_angstrom: float,
    injected_flux: float,
    snr_threshold: float = DEFAULT_LINE_RECOVERY_SNR_THRESHOLD,
) -> LineInjectionResult:
    """Fit a Gaussian-plus-offset within `window_angstrom` of the expected
    line centre and report whether the injected line was recovered above
    `snr_threshold`. Flux here is reported as the fitted Gaussian integral
    (amplitude * sigma * sqrt(2*pi)), matching `injected_flux`'s convention
    of amplitude * sigma * sqrt(2*pi) for a fair recovered/injected ratio.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    signal_with_line = np.asarray(signal_with_line, dtype=float)
    uncertainty = np.asarray(uncertainty, dtype=float)

    in_window = np.abs(wavelength - expected_center_angstrom) <= window_angstrom
    finite = np.isfinite(wavelength) & np.isfinite(signal_with_line) & np.isfinite(uncertainty)
    usable = in_window & finite & (uncertainty > 0)
    n_usable = int(np.sum(usable))
    if n_usable < 5:
        raise InsufficientDataError(
            f"only {n_usable} usable pixels in the {window_angstrom} Angstrom fit window"
        )

    x = wavelength[usable]
    y = signal_with_line[usable]
    yerr = uncertainty[usable]

    amp_guess = float(np.max(y) - np.median(y))
    p0 = [max(amp_guess, 1e-6), expected_center_angstrom, window_angstrom / 4.0, float(np.median(y))]
    try:
        popt, pcov = curve_fit(
            _gaussian, x, y, p0=p0, sigma=yerr, absolute_sigma=True, maxfev=5000
        )
    except (RuntimeError, ValueError) as exc:
        raise ConvergenceError(f"Gaussian line fit failed to converge: {exc}") from exc

    amplitude, _center, sigma, _offset = popt
    perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(4, np.nan)

    recovered_flux = float(amplitude * abs(sigma) * np.sqrt(2.0 * np.pi))
    # propagate amplitude/sigma uncertainty into the flux uncertainty (partial derivatives
    # of flux = amplitude * sigma * sqrt(2*pi), ignoring the amplitude-sigma covariance term)

    amp_err, sigma_err = perr[0], perr[2]
    flux_err = float(
        np.sqrt(2.0 * np.pi)
        * np.sqrt((sigma * amp_err) ** 2 + (amplitude * sigma_err) ** 2)
    ) if np.isfinite(amp_err) and np.isfinite(sigma_err) else float("nan")

    recovered_snr = recovered_flux / flux_err if flux_err and np.isfinite(flux_err) and flux_err > 0 else 0.0
    detected = bool(np.isfinite(recovered_snr) and recovered_snr >= snr_threshold)
    flux_ratio = recovered_flux / injected_flux if injected_flux != 0 else float("nan")

    return LineInjectionResult(
        injected_flux=float(injected_flux),
        recovered_flux=recovered_flux,
        recovered_flux_err=flux_err,
        recovered_snr=float(recovered_snr) if np.isfinite(recovered_snr) else 0.0,
        detected=detected,
        flux_ratio=flux_ratio,
    )


def line_recovery_rate(
    best_line_snr: np.ndarray, snr_threshold: float = DEFAULT_LINE_RECOVERY_SNR_THRESHOLD
) -> float:
    """Fraction of objects whose best released line S/N (spe_line_snr_gf)
    clears `snr_threshold` — the real-catalogue line-recovery indicator.
    """
    best_line_snr = np.asarray(best_line_snr, dtype=float)
    finite = best_line_snr[np.isfinite(best_line_snr)]
    if finite.size == 0:
        raise InsufficientDataError("no finite best_line_snr values to compute a recovery rate")
    return float(np.mean(finite >= snr_threshold))


__all__ = [
    "DEFAULT_LINE_RECOVERY_SNR_THRESHOLD",
    "LineInjectionResult",
    "inject_gaussian_line",
    "line_recovery_rate",
    "recover_injected_line",
]

"""Synthetic NISP-like spectra and catalogue generation, clearly labelled and
never presented as real data. Shared by the test suite (injection-recovery
validation, negative control) and `scripts/run_analysis.py --demo` /
`scripts/make_figures.py --demo`, so demo and test fixtures never duplicate
this logic.

The synthetic contamination effect injected here (higher noise + a spurious
second-trace continuum bump for the "flagged" group) is a stand-in used only
to validate that the pipeline can detect a *known* injected effect; it is not
a model of the real Euclid SIR contamination physics.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from euclid_nisp_contamination_audit.masks import CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED

NISP_WAVELENGTH_MIN_ANGSTROM = 12000.0  # verified: NISP red grism range ~1.25-1.85 micron (2503.15307)
NISP_WAVELENGTH_MAX_ANGSTROM = 18500.0
DEFAULT_N_PIXELS = 400
DEFAULT_CONTAMINATION_FRACTION = 0.18  # matches the live-observed ~16% flagged fraction on tile 102160061


@dataclass(frozen=True)
class SyntheticSpectrum:
    object_id: int
    wavelength_angstrom: np.ndarray
    signal: np.ndarray
    uncertainty: np.ndarray
    mask: np.ndarray
    quality: np.ndarray
    contamination_group: str
    true_z: float


def build_synthetic_spectrum(
    *,
    object_id: int,
    seed: int,
    contaminated: bool,
    true_z: float = 1.0,
    continuum_level: float = 2.0,
    line_center_rest_angstrom: float = 6564.6,  # H-alpha, rest-frame
    line_amplitude: float = 0.0,
    n_pixels: int = DEFAULT_N_PIXELS,
) -> SyntheticSpectrum:
    """One synthetic NISP-like combined 1-D spectrum with the same column
    layout as the real IRSA VOTable product (WAVELENGTH/SIGNAL/UNCERTAINTY/
    MASK/QUALITY). Contaminated spectra get ~3x noise and a broad spurious
    continuum bump (a stand-in for an overlapping second-order trace).
    """
    rng = np.random.default_rng(seed)
    wavelength = np.linspace(NISP_WAVELENGTH_MIN_ANGSTROM, NISP_WAVELENGTH_MAX_ANGSTROM, n_pixels)

    base_uncertainty = np.full(n_pixels, 0.15 if not contaminated else 0.45)
    continuum = np.full(n_pixels, continuum_level)

    observed_line_center = line_center_rest_angstrom * (1.0 + true_z)
    if line_amplitude > 0 and NISP_WAVELENGTH_MIN_ANGSTROM < observed_line_center < NISP_WAVELENGTH_MAX_ANGSTROM:
        continuum = continuum + line_amplitude * np.exp(
            -0.5 * ((wavelength - observed_line_center) / 25.0) ** 2
        )

    if contaminated:
        bump_center = rng.uniform(NISP_WAVELENGTH_MIN_ANGSTROM, NISP_WAVELENGTH_MAX_ANGSTROM)
        continuum = continuum + 0.6 * np.exp(-0.5 * ((wavelength - bump_center) / 900.0) ** 2)

    noise = rng.normal(loc=0.0, scale=base_uncertainty)
    signal = continuum + noise

    mask = np.zeros(n_pixels, dtype=int)
    quality = np.zeros(n_pixels, dtype=int)
    edge = n_pixels // 40
    if edge > 0:
        mask[:edge] = 1
        mask[-edge:] = 1

    group = CONTAMINATION_GROUP_FLAGGED if contaminated else CONTAMINATION_GROUP_CLEAN
    return SyntheticSpectrum(
        object_id=object_id,
        wavelength_angstrom=wavelength,
        signal=signal,
        uncertainty=base_uncertainty,
        mask=mask,
        quality=quality,
        contamination_group=group,
        true_z=true_z,
    )


@dataclass(frozen=True)
class SyntheticCatalogue:
    object_id: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    context_warning_flag: np.ndarray
    context_error_flag: np.ndarray
    contamination_group: np.ndarray
    spe_z: np.ndarray
    spe_z_err: np.ndarray
    spe_z_true: np.ndarray
    spe_z_rel: np.ndarray
    spe_cont_snr: np.ndarray
    best_line_snr: np.ndarray


def build_synthetic_catalogue(
    *,
    n_sources: int = 300,
    seed: int = 20260713,
    contamination_fraction: float = DEFAULT_CONTAMINATION_FRACTION,
    ra0: float = 60.0,
    dec0: float = -48.0,
) -> SyntheticCatalogue:
    """Deterministic synthetic Euclid-Q1-NISP-catalogue-like table with a
    known, injected contamination effect: the flagged group has a lower
    continuum S/N, lower line-recovery S/N, lower redshift reliability, and a
    higher redshift outlier rate than the clean group, by construction — so
    the pipeline's job in validation is to *recover* this known signal.
    """
    if n_sources < 8:
        raise ValueError("n_sources must be at least 8")
    if not (0.0 < contamination_fraction < 1.0):
        raise ValueError("contamination_fraction must be in (0, 1)")

    rng = np.random.default_rng(seed)
    object_id = np.arange(1, n_sources + 1, dtype=np.int64)
    ra = ra0 + rng.uniform(-0.05, 0.05, size=n_sources)
    dec = dec0 + rng.uniform(-0.05, 0.05, size=n_sources)

    flagged = rng.uniform(size=n_sources) < contamination_fraction
    context_warning_flag = np.where(flagged, 32768, 0)
    context_error_flag = np.where(rng.uniform(size=n_sources) < contamination_fraction * 0.05, 64, 0)
    contamination_group = np.where(flagged, CONTAMINATION_GROUP_FLAGGED, CONTAMINATION_GROUP_CLEAN)

    z_true = rng.uniform(0.5, 1.8, size=n_sources)
    z_scatter = np.where(flagged, 0.06, 0.01)
    spe_z = z_true + rng.normal(0.0, z_scatter)
    spe_z_err = np.abs(rng.normal(loc=z_scatter, scale=z_scatter * 0.2))

    z_rel_center = np.where(flagged, 0.45, 0.85)
    spe_z_rel = np.clip(rng.normal(loc=z_rel_center, scale=0.12), 0.0, 1.0)

    cont_snr_center = np.where(flagged, 4.0, 12.0)
    spe_cont_snr = np.clip(rng.normal(loc=cont_snr_center, scale=cont_snr_center * 0.25), 0.1, None)

    line_snr_center = np.where(flagged, 3.0, 9.0)
    best_line_snr = np.clip(rng.normal(loc=line_snr_center, scale=line_snr_center * 0.3), 0.0, None)

    return SyntheticCatalogue(
        object_id=object_id,
        ra=ra,
        dec=dec,
        context_warning_flag=context_warning_flag,
        context_error_flag=context_error_flag,
        contamination_group=contamination_group,
        spe_z=spe_z,
        spe_z_err=spe_z_err,
        spe_z_true=z_true,
        spe_z_rel=spe_z_rel,
        spe_cont_snr=spe_cont_snr,
        best_line_snr=best_line_snr,
    )


__all__ = [
    "DEFAULT_CONTAMINATION_FRACTION",
    "NISP_WAVELENGTH_MAX_ANGSTROM",
    "NISP_WAVELENGTH_MIN_ANGSTROM",
    "SyntheticCatalogue",
    "SyntheticSpectrum",
    "build_synthetic_catalogue",
    "build_synthetic_spectrum",
]

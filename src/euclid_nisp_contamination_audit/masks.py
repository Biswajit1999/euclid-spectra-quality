"""Contamination-group labelling and finite/quality masking.

The released Euclid Q1 NISP `spe_context_warning_flag` / `spe_context_error_flag`
columns are the operational "contamination indicator" for this audit: the SIR
processing paper (arXiv:2503.15307) states the pipeline "subtracts
cross-contaminations, minimizes self-contamination" during spectral
extraction, and the "context" flags are the only released per-object quality
bits plausibly tied to that step. The exact bit semantics are not published
in the queryable TAP column metadata (confirmed live — the `description`
field for these columns is a one-line label, not a bit dictionary), so this
mapping is treated as a documented assumption, not an internal Euclid
definition (see docs/ASSUMPTIONS_AND_LIMITATIONS.md).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from euclid_nisp_contamination_audit.exceptions import DataSchemaError

CONTAMINATION_GROUP_CLEAN = "clean"
CONTAMINATION_GROUP_FLAGGED = "flagged"


@dataclass(frozen=True)
class ThresholdDefinition:
    """One candidate definition of 'contaminated' used for threshold-sensitivity testing."""

    name: str
    description: str


THRESHOLD_DEFINITIONS: tuple[ThresholdDefinition, ...] = (
    ThresholdDefinition("context_warning_only", "spe_context_warning_flag != 0"),
    ThresholdDefinition(
        "context_warning_or_error",
        "spe_context_warning_flag != 0 OR spe_context_error_flag != 0",
    ),
    ThresholdDefinition(
        "context_error_only_strict",
        "spe_context_error_flag != 0 (most conservative: severe contamination only)",
    ),
)


def contamination_group_labels(
    context_warning_flag: np.ndarray,
    context_error_flag: np.ndarray,
    threshold: str = "context_warning_only",
) -> np.ndarray:
    """Return an array of 'clean'/'flagged' labels for one of THRESHOLD_DEFINITIONS."""
    warning = np.asarray(context_warning_flag)
    error = np.asarray(context_error_flag)
    if warning.shape != error.shape:
        raise DataSchemaError("context_warning_flag and context_error_flag must have the same shape")

    valid_names = {t.name for t in THRESHOLD_DEFINITIONS}
    if threshold not in valid_names:
        raise DataSchemaError(f"unknown threshold '{threshold}', expected one of {valid_names}")

    if threshold == "context_warning_only":
        flagged = warning != 0
    elif threshold == "context_warning_or_error":
        flagged = (warning != 0) | (error != 0)
    else:  # context_error_only_strict
        flagged = error != 0

    return np.where(flagged, CONTAMINATION_GROUP_FLAGGED, CONTAMINATION_GROUP_CLEAN)


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    """Elementwise AND of `np.isfinite` across all input arrays (broadcastable, same shape)."""
    if not arrays:
        raise DataSchemaError("finite_mask requires at least one array")
    shape = arrays[0].shape
    mask = np.ones(shape, dtype=bool)
    for arr in arrays:
        arr = np.asarray(arr, dtype=float)
        if arr.shape != shape:
            raise DataSchemaError(f"array shape mismatch in finite_mask: {arr.shape} vs {shape}")
        mask &= np.isfinite(arr)
    return mask


def spectral_pixel_mask(mask_column: np.ndarray, quality_column: np.ndarray) -> np.ndarray:
    """Boolean 'good pixel' mask from the VOTable MASK/QUALITY columns.

    Convention (consistent with the released Euclid combined-spectrum
    product): MASK == 0 and QUALITY == 0 marks an unflagged science pixel.
    Any nonzero value in either column excludes the pixel.
    """
    mask_column = np.asarray(mask_column)
    quality_column = np.asarray(quality_column)
    if mask_column.shape != quality_column.shape:
        raise DataSchemaError("MASK and QUALITY columns must have the same shape")
    return (mask_column == 0) & (quality_column == 0)


__all__ = [
    "CONTAMINATION_GROUP_CLEAN",
    "CONTAMINATION_GROUP_FLAGGED",
    "THRESHOLD_DEFINITIONS",
    "ThresholdDefinition",
    "contamination_group_labels",
    "finite_mask",
    "spectral_pixel_mask",
]

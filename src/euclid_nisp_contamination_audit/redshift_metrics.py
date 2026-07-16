"""Redshift-reliability metrics: normalised error, outlier fraction, and
released reliability-flag summaries, split by contamination group.
"""
from __future__ import annotations

import numpy as np

from euclid_nisp_contamination_audit.exceptions import InsufficientDataError

DEFAULT_OUTLIER_THRESHOLD = 0.15  # standard |delta_z| / (1+z_true) spec-z outlier convention


def normalized_redshift_error(z_estimate: np.ndarray, z_reference: np.ndarray) -> np.ndarray:
    """Delta_z / (1 + z_reference), the standard normalised redshift error."""
    z_estimate = np.asarray(z_estimate, dtype=float)
    z_reference = np.asarray(z_reference, dtype=float)
    if z_estimate.shape != z_reference.shape:
        raise InsufficientDataError("z_estimate and z_reference must have the same shape")
    return (z_estimate - z_reference) / (1.0 + z_reference)


def outlier_fraction(
    delta_z_norm: np.ndarray, threshold: float = DEFAULT_OUTLIER_THRESHOLD
) -> float:
    """Fraction of |delta_z / (1+z)| exceeding `threshold`."""
    delta_z_norm = np.asarray(delta_z_norm, dtype=float)
    finite = delta_z_norm[np.isfinite(delta_z_norm)]
    if finite.size == 0:
        raise InsufficientDataError("no finite normalised redshift errors to compute an outlier fraction")
    return float(np.mean(np.abs(finite) > threshold))


def reliability_summary_by_group(
    z_rel: np.ndarray, group_labels: np.ndarray
) -> dict[str, dict[str, float]]:
    """Per-contamination-group summary of the released `spe_z_rel` reliability
    score: mean, median, and fraction below the conventional 'unreliable'
    threshold of 0.5 (Euclid SPE redshift-reliability convention: spe_z_rel is
    a probability-like score in [0, 1]).
    """
    z_rel = np.asarray(z_rel, dtype=float)
    group_labels = np.asarray(group_labels)
    if z_rel.shape != group_labels.shape:
        raise InsufficientDataError("z_rel and group_labels must have the same shape")

    summary: dict[str, dict[str, float]] = {}
    for group in sorted(set(group_labels.tolist())):
        in_group = group_labels == group
        values = z_rel[in_group]
        values = values[np.isfinite(values)]
        if values.size == 0:
            continue
        summary[str(group)] = {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "fraction_unreliable": float(np.mean(values < 0.5)),
            "n": int(values.size),
        }
    return summary


__all__ = [
    "DEFAULT_OUTLIER_THRESHOLD",
    "normalized_redshift_error",
    "outlier_fraction",
    "reliability_summary_by_group",
]

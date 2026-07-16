from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.exceptions import DataSchemaError
from euclid_nisp_contamination_audit.masks import (
    CONTAMINATION_GROUP_CLEAN,
    CONTAMINATION_GROUP_FLAGGED,
    contamination_group_labels,
    finite_mask,
    spectral_pixel_mask,
)


def test_contamination_group_labels_basic() -> None:
    warning = np.array([0, 32768, 0, 32768])
    error = np.array([0, 0, 0, 0])
    labels = contamination_group_labels(warning, error)
    assert list(labels) == [
        CONTAMINATION_GROUP_CLEAN,
        CONTAMINATION_GROUP_FLAGGED,
        CONTAMINATION_GROUP_CLEAN,
        CONTAMINATION_GROUP_FLAGGED,
    ]


def test_contamination_group_labels_unknown_threshold() -> None:
    with pytest.raises(DataSchemaError):
        contamination_group_labels(np.array([0]), np.array([0]), threshold="nonexistent")


def test_contamination_group_labels_shape_mismatch() -> None:
    with pytest.raises(DataSchemaError):
        contamination_group_labels(np.array([0, 1]), np.array([0]))


def test_finite_mask() -> None:
    a = np.array([1.0, np.nan, 3.0])
    b = np.array([1.0, 2.0, np.inf])
    mask = finite_mask(a, b)
    assert list(mask) == [True, False, False]


def test_finite_mask_requires_arrays() -> None:
    with pytest.raises(DataSchemaError):
        finite_mask()


def test_spectral_pixel_mask() -> None:
    mask_col = np.array([0, 1, 0])
    quality_col = np.array([0, 0, 1])
    good = spectral_pixel_mask(mask_col, quality_col)
    assert list(good) == [True, False, False]

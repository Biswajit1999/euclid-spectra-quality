from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.exceptions import InsufficientDataError
from euclid_nisp_contamination_audit.redshift_metrics import (
    normalized_redshift_error,
    outlier_fraction,
    reliability_summary_by_group,
)


def test_normalized_redshift_error() -> None:
    z_est = np.array([1.1, 2.0])
    z_ref = np.array([1.0, 2.0])
    delta = normalized_redshift_error(z_est, z_ref)
    assert delta[0] == pytest.approx(0.1 / 2.0)
    assert delta[1] == pytest.approx(0.0)


def test_outlier_fraction() -> None:
    delta = np.array([0.01, 0.5, -0.3, 0.02, np.nan])
    frac = outlier_fraction(delta, threshold=0.15)
    assert frac == pytest.approx(2 / 4)


def test_outlier_fraction_no_finite() -> None:
    with pytest.raises(InsufficientDataError):
        outlier_fraction(np.array([np.nan]))


def test_reliability_summary_by_group() -> None:
    z_rel = np.array([0.9, 0.8, 0.2, 0.1])
    groups = np.array(["clean", "clean", "flagged", "flagged"])
    summary = reliability_summary_by_group(z_rel, groups)
    assert summary["clean"]["n"] == 2
    assert summary["flagged"]["fraction_unreliable"] == 1.0
    assert summary["clean"]["fraction_unreliable"] == 0.0

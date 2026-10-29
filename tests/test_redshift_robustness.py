"""Synthetic truth-labelled regression tests (not Euclid Q1 measurements)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))
from redshift_robustness import (  # noqa: E402
    matched_snr_outlier_contrast, redshift_coverage, reliability_calibration,
)


def synthetic(seed: int = 23):
    rng = np.random.default_rng(seed)
    n = 600
    groups = np.array(["clean"] * 300 + ["flagged"] * 300)
    truth = rng.uniform(.7, 1.7, n)
    snr = rng.uniform(4, 16, n)  # same SNR support in both groups
    sigma = .02 * (1 + truth)
    estimate = truth + rng.normal(0, sigma)
    contaminated = (groups == "flagged") & (rng.random(n) < .3)
    estimate[contaminated] += .5 * (1 + truth[contaminated])
    probability = np.where(groups == "clean", .9, .7)
    return estimate, truth, sigma, groups, snr, probability


def test_gaussian_pull_coverage_and_outlier_degradation():
    z, t, sigma, groups, _, _ = synthetic()
    result = redshift_coverage(z, t, sigma, groups)
    assert .55 < result["clean"]["coverage_1sigma"] < .8
    assert result["flagged"]["coverage_2sigma"] < result["clean"]["coverage_2sigma"]
    assert result["flagged"]["mean_squared_pull"] > result["clean"]["mean_squared_pull"]


def test_exact_coverage_boundary_and_validation():
    z = np.array([1., 1.1, 1.2, 1.3]); t = np.ones(4); s = np.full(4, .1)
    result = redshift_coverage(z, t, s, np.array(["x"] * 4))
    assert result["x"]["coverage_1sigma"] == pytest.approx(.5)
    with pytest.raises(ValueError):
        redshift_coverage(z, t, np.zeros(4), np.array(["x"] * 4))


def test_reliability_calibration_perfect_probabilities():
    t = np.ones(40); z = t.copy(); z[20:] = 2.
    score = np.r_[np.ones(20), np.zeros(20)]
    result = reliability_calibration(z, t, score)
    assert result["brier_score"] == pytest.approx(0.)
    assert result["expected_calibration_error"] == pytest.approx(0.)


def test_reliability_rejects_invalid_probability():
    with pytest.raises(ValueError):
        reliability_calibration(np.ones(4), np.ones(4), np.array([0, .5, 1, 1.1]))


def test_snr_matched_recovers_injected_excess():
    z, t, _, groups, snr, _ = synthetic()
    result = matched_snr_outlier_contrast(z, t, groups, snr, n_bootstrap=300)
    assert result["n_strata"] == 5
    assert result["n_matched"] > 550
    assert .18 < result["outlier_rate_difference_b_minus_a"] < .42
    assert result["ci_95_percentile"][0] > 0


def test_bootstrap_is_deterministic():
    z, t, _, groups, snr, _ = synthetic()
    a = matched_snr_outlier_contrast(z, t, groups, snr, n_bootstrap=40, seed=11)
    b = matched_snr_outlier_contrast(z, t, groups, snr, n_bootstrap=40, seed=11)
    assert a == b


def test_degenerate_real_tile_fails_closed():
    z, t, _, groups, snr, _ = synthetic()
    with pytest.raises(ValueError, match="degenerate"):
        matched_snr_outlier_contrast(z, t, np.array(["clean"] * len(groups)), snr)


def test_no_snr_overlap_fails_closed():
    z, t, _, groups, snr, _ = synthetic()
    snr[groups == "flagged"] = 100 + snr[groups == "flagged"]
    with pytest.raises(ValueError, match="no SNR stratum"):
        matched_snr_outlier_contrast(z, t, groups, snr)

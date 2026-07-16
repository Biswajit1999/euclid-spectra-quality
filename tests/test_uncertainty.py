from __future__ import annotations

import numpy as np
import pytest

from euclid_nisp_contamination_audit.exceptions import ConvergenceError, InsufficientDataError
from euclid_nisp_contamination_audit.uncertainty import (
    bootstrap_statistic,
    check_fit_convergence,
    permuted_label_null_distribution,
)


def test_bootstrap_statistic_ci_contains_estimate() -> None:
    rng = np.random.default_rng(1)
    data = rng.normal(5.0, 1.0, size=200)
    result = bootstrap_statistic(data, np.mean, n_resamples=500, seed=20260713)
    assert result.ci_low <= result.estimate <= result.ci_high
    assert result.n_resamples == 500


def test_bootstrap_statistic_insufficient_data() -> None:
    with pytest.raises(InsufficientDataError):
        bootstrap_statistic(np.array([1.0]), np.mean, n_resamples=10, seed=1)


def test_check_fit_convergence_well_conditioned() -> None:
    pcov = np.eye(3) * 1e-4
    result = check_fit_convergence(pcov, residuals=np.array([0.1, -0.1, 0.05]), dof=3)
    assert result.converged
    assert result.reduced_chi_square is not None


def test_check_fit_convergence_ill_conditioned_raises() -> None:
    pcov = np.array([[1e10, 0.0], [0.0, 1e-10]])
    with pytest.raises(ConvergenceError):
        check_fit_convergence(pcov, max_condition_number=1e6)


def test_check_fit_convergence_nonfinite_raises() -> None:
    pcov = np.array([[np.nan, 0.0], [0.0, 1.0]])
    with pytest.raises(ConvergenceError):
        check_fit_convergence(pcov)


def test_permuted_label_null_distribution_detects_real_effect() -> None:
    rng = np.random.default_rng(2)
    clean = rng.normal(10.0, 1.0, size=100)
    flagged = rng.normal(4.0, 1.0, size=100)
    values = np.concatenate([clean, flagged])
    labels = np.array(["clean"] * 100 + ["flagged"] * 100)

    result = permuted_label_null_distribution(values, labels, "flagged", "clean", n_permutations=500, seed=3)
    assert result.p_value_two_sided < 0.01
    assert abs(result.observed_difference) > 3 * result.null_std


def test_permuted_label_null_distribution_no_real_effect() -> None:
    # Labels are independent of `values` by construction, so the observed
    # group difference should sit within the bulk of its own null
    # distribution; a >=5% two-sided test has an intrinsic ~1-in-20 false
    # positive rate, so this asserts against a generous 3-sigma-equivalent
    # bound rather than a fixed alpha to avoid a rare, non-representative seed.
    rng = np.random.default_rng(4)
    values = rng.normal(5.0, 1.0, size=200)
    labels = rng.permutation(np.array(["clean"] * 100 + ["flagged"] * 100))

    result = permuted_label_null_distribution(values, labels, "flagged", "clean", n_permutations=2000, seed=5)
    assert abs(result.observed_difference - result.null_mean) < 3 * result.null_std


def test_permuted_label_null_distribution_insufficient_data() -> None:
    with pytest.raises(InsufficientDataError):
        permuted_label_null_distribution(
            np.array([1.0]), np.array(["clean"]), "flagged", "clean", n_permutations=10, seed=1
        )

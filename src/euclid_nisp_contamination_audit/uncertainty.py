"""Uncertainty quantification, split into two explicit categories.

- Observational/statistical uncertainty: bootstrap resampling over sources
  (`bootstrap_statistic`).
- Numerical/convergence uncertainty: fit covariance conditioning and reduced
  chi-square from the nonlinear least-squares solver
  (`check_fit_convergence`).

Reported separately per docs/VALIDATION_CONTRACT.md's requirement to not conflate
observational uncertainty with numerical convergence uncertainty.

`permuted_label_null_distribution` implements the permuted-label negative
control required by docs/VALIDATION_CONTRACT.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from euclid_nisp_contamination_audit.exceptions import ConvergenceError, InsufficientDataError


@dataclass(frozen=True)
class BootstrapResult:
    estimate: float
    ci_low: float
    ci_high: float
    n_resamples: int
    confidence_level: float


def bootstrap_statistic(
    data: np.ndarray,
    statistic: Callable[[np.ndarray], float],
    n_resamples: int,
    seed: int,
    confidence_level: float = 0.95,
) -> BootstrapResult:
    """Bootstrap a scalar statistic over `data`, resampling with replacement."""
    arr = np.asarray(data, dtype=float)
    if arr.size < 2:
        raise InsufficientDataError(f"need at least 2 samples to bootstrap, got {arr.size}")
    if not (0.0 < confidence_level < 1.0):
        raise ValueError("confidence_level must be in (0, 1)")

    rng = np.random.default_rng(seed)
    n = arr.size
    resampled = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        resampled[i] = statistic(arr[idx])

    alpha = 1 - confidence_level
    lo, hi = np.percentile(resampled, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return BootstrapResult(
        estimate=float(statistic(arr)),
        ci_low=float(lo),
        ci_high=float(hi),
        n_resamples=n_resamples,
        confidence_level=confidence_level,
    )


@dataclass(frozen=True)
class FitConvergence:
    converged: bool
    covariance_condition_number: float
    reduced_chi_square: float | None


def check_fit_convergence(
    pcov: np.ndarray,
    residuals: np.ndarray | None = None,
    dof: int | None = None,
    max_condition_number: float = 1e10,
) -> FitConvergence:
    """Assess numerical convergence of a `scipy.optimize.curve_fit` result.

    Raises ConvergenceError (not a silently-returned False) if the
    covariance matrix is non-finite or ill-conditioned beyond
    `max_condition_number`, per the docs/VALIDATION_CONTRACT.md stop condition for failed
    numerical convergence.
    """
    pcov = np.asarray(pcov, dtype=float)
    if pcov.ndim != 2 or pcov.shape[0] != pcov.shape[1]:
        raise ConvergenceError(f"covariance matrix has invalid shape {pcov.shape}")
    if not np.all(np.isfinite(pcov)):
        raise ConvergenceError("fit covariance matrix contains non-finite values")

    condition_number = float(np.linalg.cond(pcov))
    if condition_number > max_condition_number or not np.isfinite(condition_number):
        raise ConvergenceError(
            f"fit covariance condition number {condition_number:.3e} exceeds "
            f"threshold {max_condition_number:.3e}"
        )

    reduced_chi_square = None
    if residuals is not None and dof:
        if dof <= 0:
            raise ConvergenceError(f"degrees of freedom must be positive, got {dof}")
        reduced_chi_square = float(np.sum(np.asarray(residuals, dtype=float) ** 2) / dof)

    return FitConvergence(
        converged=True,
        covariance_condition_number=condition_number,
        reduced_chi_square=reduced_chi_square,
    )


@dataclass(frozen=True)
class PermutationTestResult:
    observed_difference: float
    null_mean: float
    null_std: float
    p_value_two_sided: float
    n_permutations: int


def permuted_label_null_distribution(
    values: np.ndarray,
    group_labels: np.ndarray,
    group_a: str,
    group_b: str,
    n_permutations: int,
    seed: int,
) -> PermutationTestResult:
    """Negative control: shuffle `group_labels` `n_permutations` times and
    recompute mean(group_a) - mean(group_b) each time, to build a null
    distribution for the observed group difference. A real contamination
    effect should place the observed (unpermuted) difference far in the tail
    of this null; if it does not, the group difference is statistically
    indistinguishable from a random label split.
    """
    values = np.asarray(values, dtype=float)
    group_labels = np.asarray(group_labels)
    if values.shape != group_labels.shape:
        raise InsufficientDataError("values and group_labels must have the same shape")

    mask_a = group_labels == group_a
    mask_b = group_labels == group_b
    if int(np.sum(mask_a)) < 2 or int(np.sum(mask_b)) < 2:
        raise InsufficientDataError(
            f"need at least 2 samples in each of '{group_a}' and '{group_b}' for a permutation test"
        )

    observed_difference = float(np.mean(values[mask_a]) - np.mean(values[mask_b]))

    rng = np.random.default_rng(seed)
    in_either = mask_a | mask_b
    sub_values = values[in_either]
    sub_labels = group_labels[in_either].copy()
    n_a = int(np.sum(mask_a))

    null_diffs = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        permuted = rng.permutation(sub_labels.size)
        perm_a = sub_values[permuted[:n_a]]
        perm_b = sub_values[permuted[n_a:]]
        null_diffs[i] = np.mean(perm_a) - np.mean(perm_b)

    null_mean = float(np.mean(null_diffs))
    null_std = float(np.std(null_diffs))
    p_value = float(np.mean(np.abs(null_diffs) >= np.abs(observed_difference)))

    return PermutationTestResult(
        observed_difference=observed_difference,
        null_mean=null_mean,
        null_std=null_std,
        p_value_two_sided=p_value,
        n_permutations=n_permutations,
    )


__all__ = [
    "BootstrapResult",
    "FitConvergence",
    "PermutationTestResult",
    "bootstrap_statistic",
    "check_fit_convergence",
    "permuted_label_null_distribution",
]

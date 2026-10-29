"""Truth-labelled, synthetic-first Euclid NISP reliability diagnostics.

The selected Q1 real tile has no flagged objects and no external spectroscopic
truth in this repository: these functions MUST NOT be used to manufacture a
real contamination contrast or truth-calibration claim.
"""
from __future__ import annotations

import numpy as np


def _vectors(*arrays: np.ndarray) -> list[np.ndarray]:
    result = [np.asarray(a) for a in arrays]
    if any(a.ndim != 1 or a.shape != result[0].shape for a in result):
        raise ValueError("inputs must be one-dimensional and have matching shapes")
    if len(result[0]) < 2:
        raise ValueError("at least two objects are required")
    return result


def redshift_coverage(z_est: np.ndarray, z_true: np.ndarray, sigma_z: np.ndarray,
                      groups: np.ndarray) -> dict:
    """Empirical 1σ/2σ coverage and mean squared standardized residual by group.

    Requires *external* reference redshifts; never use model z as its own truth.
    """
    z, truth, sigma, label = _vectors(z_est, z_true, sigma_z, groups)
    if not (np.all(np.isfinite(z)) and np.all(np.isfinite(truth)) and
            np.all(np.isfinite(sigma)) and np.all(sigma > 0) and np.all(truth > -1)):
        raise ValueError("finite redshifts, positive errors and z_true > -1 required")
    residual = (z - truth) / sigma
    return {str(g): {"n": int(np.sum(label == g)),
                     "coverage_1sigma": float(np.mean(np.abs(residual[label == g]) <= 1 + 1e-12)),
                     "coverage_2sigma": float(np.mean(np.abs(residual[label == g]) <= 2 + 1e-12)),
                     "mean_squared_pull": float(np.mean(residual[label == g] ** 2)),
                     "mean_pull": float(np.mean(residual[label == g]))}
            for g in sorted(set(label.tolist())) if np.sum(label == g) >= 2}


def reliability_calibration(z_est: np.ndarray, z_true: np.ndarray,
                            reliability: np.ndarray, *, threshold: float = 0.15,
                            n_bins: int = 5) -> dict:
    """Brier score and expected calibration error if score is P(non-outlier).

    These are empirical diagnostics of the declared probability interpretation,
    not a claim that Euclid spe_z_rel is already calibrated this way.
    """
    z, truth, score = _vectors(z_est, z_true, reliability)
    if (not np.all(np.isfinite(z)) or not np.all(np.isfinite(truth)) or
        not np.all(np.isfinite(score)) or np.any(truth <= -1) or
        np.any((score < 0) | (score > 1)) or not np.isfinite(threshold) or
        threshold <= 0 or n_bins < 2):
        raise ValueError("invalid reference, probability, threshold or bins")
    correct = (np.abs((z - truth) / (1 + truth)) <= threshold).astype(float)
    edges = np.linspace(0, 1, n_bins + 1)
    assigned = np.minimum(np.searchsorted(edges, score, side="right") - 1, n_bins - 1)
    summary = []
    ece = 0.
    for j in range(n_bins):
        mask = assigned == j
        if not np.any(mask):
            continue
        confidence, observed = float(np.mean(score[mask])), float(np.mean(correct[mask]))
        count = int(np.sum(mask))
        ece += count / len(score) * abs(confidence - observed)
        summary.append({"bin": j, "n": count, "mean_score": confidence,
                        "observed_non_outlier_fraction": observed})
    return {"n": len(score), "threshold": threshold,
            "brier_score": float(np.mean((score - correct) ** 2)),
            "expected_calibration_error": float(ece), "bins": summary}


def matched_snr_outlier_contrast(z_est: np.ndarray, z_true: np.ndarray,
                                 groups: np.ndarray, snr: np.ndarray, *,
                                 group_a: str = "clean", group_b: str = "flagged",
                                 threshold: float = 0.15, n_bins: int = 5,
                                 n_bootstrap: int = 400, seed: int = 20261029) -> dict:
    """SNR-stratified (B - A) outlier-rate contrast, within-cell bootstrap.

    Quantile strata without >=3 members of *both* groups are excluded. The
    target is the pooled *matched* subset, not the full catalogue population.
    """
    z, truth, label, quality = _vectors(z_est, z_true, groups, snr)
    if (not np.all(np.isfinite(z)) or not np.all(np.isfinite(truth)) or
        not np.all(np.isfinite(quality)) or np.any(truth <= -1) or
        np.any(quality <= 0) or threshold <= 0 or not np.isfinite(threshold) or
        n_bins < 2 or n_bootstrap < 20 or group_a == group_b):
        raise ValueError("invalid data or bootstrap settings")
    chosen = (label == group_a) | (label == group_b)
    if np.sum(label == group_a) < 3 or np.sum(label == group_b) < 3:
        raise ValueError("degenerate contamination groups")
    # Exclude all objects outside *common support* before quantile matching.
    # Otherwise a quantile stratum spanning an empty SNR gap could falsely
    # present two disjoint populations as matched.
    support_low = max(np.min(quality[label == group_a]), np.min(quality[label == group_b]))
    support_high = min(np.max(quality[label == group_a]), np.max(quality[label == group_b]))
    if support_low >= support_high:
        raise ValueError("no SNR stratum has adequate overlap between groups")
    chosen &= (quality >= support_low) & (quality <= support_high)
    boundaries = np.quantile(quality[chosen], np.linspace(0, 1, n_bins + 1))
    if np.any(np.diff(boundaries) <= 0):
        raise ValueError("non-distinct SNR quantile boundaries")
    bins = np.minimum(np.searchsorted(boundaries, quality, side="right") - 1, n_bins - 1)
    is_outlier = (np.abs((z - truth) / (1 + truth)) > threshold).astype(float)
    cells = []
    for j in range(n_bins):
        aa = is_outlier[chosen & (bins == j) & (label == group_a)]
        bb = is_outlier[chosen & (bins == j) & (label == group_b)]
        if len(aa) >= 3 and len(bb) >= 3:
            cells.append((j, aa, bb))
    if not cells:
        raise ValueError("no SNR stratum has adequate overlap between groups")
    weights = np.array([len(a) + len(b) for _, a, b in cells], dtype=float)
    weights /= np.sum(weights)
    contrasts = np.array([np.mean(b) - np.mean(a) for _, a, b in cells])
    rng = np.random.default_rng(seed)
    trials = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        draws = [np.mean(b[rng.integers(len(b), size=len(b))]) -
                 np.mean(a[rng.integers(len(a), size=len(a))]) for _, a, b in cells]
        trials[i] = float(weights @ draws)
    return {
        "outlier_rate_difference_b_minus_a": float(weights @ contrasts),
        "ci_95_percentile": list(map(float, np.quantile(trials, [.025, .975]))),
        "n_matched": int(sum(len(a) + len(b) for _, a, b in cells)),
        "n_total": int(np.sum((label == group_a) | (label == group_b))),
        "n_common_support": int(np.sum(chosen)), "n_strata": len(cells),
        "strata": [{"bin": j, "n_a": len(a), "n_b": len(b),
                    "difference": float(np.mean(b) - np.mean(a))} for j, a, b in cells],
        "seed": seed, "n_bootstrap": n_bootstrap,
    }

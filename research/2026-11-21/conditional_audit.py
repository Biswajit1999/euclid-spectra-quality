"""Truth-labelled conditional redshift calibration; synthetic fixtures only."""
import numpy as np


def normalized_redshift_error(z_est, z_ref):
    a = np.asarray(z_est, dtype=float)
    b = np.asarray(z_ref, dtype=float)
    if a.shape != b.shape or a.size == 0 or not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(b <= -1):
        raise ValueError("finite matched redshifts with z_ref > -1 required")
    return np.abs(a-b)/(1+b)


def stratified_outlier_difference(errors, flagged, strata, threshold=.15, n_boot=1000, seed=0):
    """Flagged minus clean outlier rate, averaged over strata by source count.

    Bootstrap resamples independently within each flag-by-stratum cell.
    This is descriptive, not a causal estimate.
    """
    e = np.asarray(errors, dtype=float)
    f = np.asarray(flagged, dtype=bool)
    g = np.asarray(strata)
    if e.ndim != 1 or e.size == 0 or f.shape != e.shape or g.shape != e.shape or not np.isfinite(e).all() or np.any(e < 0) or threshold < 0 or n_boot < 2:
        raise ValueError("invalid paired inputs")
    cells = []
    for label in np.unique(g):
        mask = g == label
        yes = (e[mask & f] > threshold).astype(float)
        no = (e[mask & ~f] > threshold).astype(float)
        if not len(yes) or not len(no):
            raise ValueError("both groups required in every stratum")
        cells.append((yes, no, np.count_nonzero(mask)/len(e)))
    estimate = sum(w*(a.mean()-b.mean()) for a,b,w in cells)
    rng = np.random.default_rng(seed)
    draws = np.zeros(n_boot)
    for a,b,w in cells:
        ia = rng.integers(len(a),size=(n_boot,len(a)))
        ib = rng.integers(len(b),size=(n_boot,len(b)))
        draws += w*(a[ia].mean(axis=1)-b[ib].mean(axis=1))
    return float(estimate), np.quantile(draws,[.025,.975])


def interval_coverage(z_est, z_ref, sigma, multiplier=1.):
    """Conditional empirical coverage, assuming sigma is externally calibrated."""
    a = np.asarray(z_est,dtype=float)
    b = np.asarray(z_ref,dtype=float)
    s = np.asarray(sigma,dtype=float)
    if a.shape != b.shape or a.shape != s.shape or a.size == 0 or not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(s).all() or np.any(s <= 0) or multiplier <= 0:
        raise ValueError("invalid estimates, truth, sigma or multiplier")
    return float(np.mean(np.abs(a-b) <= multiplier*s))

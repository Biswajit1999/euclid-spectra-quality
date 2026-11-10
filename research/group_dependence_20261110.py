"""Synthetic-only dependence and weighted redshift residual diagnostics.

A duplicated source is not an independent redshift validation target.
"""
from collections import defaultdict
from math import isfinite, sqrt


def source_level_pull_summary(source_ids, z_fit, z_truth, z_sigma):
    """One inverse-variance combined residual per source; report RMS pulls.

    This is a descriptive sensitivity statistic, not a calibrated coverage test:
    within-source covariance is unknown, so combined sigma is NOT estimated.
    """
    n = len(source_ids)
    if not n or any(len(v) != n for v in (z_fit, z_truth, z_sigma)):
        raise ValueError("input lengths must agree and be nonempty")
    grouped = defaultdict(list)
    for sid, fit, truth, sigma in zip(source_ids, z_fit, z_truth, z_sigma):
        if not all(isfinite(float(v)) for v in (fit, truth, sigma)) or sigma <= 0:
            raise ValueError("invalid redshift or uncertainty")
        grouped[str(sid)].append(((fit-truth)/sigma, 1/sigma**2))
    per_source = []
    for sid in sorted(grouped):
        rows = grouped[sid]
        per_source.append(sum(p*w for p,w in rows)/sum(w for _,w in rows))
    return {
        "rows": n,
        "unique_sources": len(grouped),
        "source_mean_pull": sum(per_source)/len(per_source),
        "source_rms_pull": sqrt(sum(p*p for p in per_source)/len(per_source)),
    }


def design_effect(n, cluster_size, intraclass_correlation):
    """Equal-sized clusters: n_eff = n/[1+(m-1)rho]."""
    if not isinstance(n, int) or not isinstance(cluster_size, int):
        raise ValueError("n and cluster_size must be integers")
    if n < 1 or cluster_size < 1 or n % cluster_size:
        raise ValueError("n must be positive and divisible by cluster_size")
    if not isfinite(intraclass_correlation):
        raise ValueError("rho must be finite")
    lower = -1/(cluster_size-1) if cluster_size > 1 else -1
    if not lower < intraclass_correlation <= 1:
        raise ValueError("rho outside strictly positive-definite range")
    deff = 1+(cluster_size-1)*intraclass_correlation
    return {"design_effect":deff, "effective_n":n/deff}

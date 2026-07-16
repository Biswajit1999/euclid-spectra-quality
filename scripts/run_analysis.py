"""Run the Euclid NISP spectral-contamination reliability audit: either the
synthetic --demo smoke path, or the real-data pipeline over
data/source_catalog.csv (written by scripts/fetch_data.py).

Peak memory is measured with the stdlib `tracemalloc` (Python-level
allocations) rather than a full process-RSS profiler such as psutil, which is
not part of this project's pinned dependency set.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import tracemalloc
from pathlib import Path

from euclid_nisp_contamination_audit import __version__
from euclid_nisp_contamination_audit.config import load_config
from euclid_nisp_contamination_audit.core import CatalogueRow, read_catalogue_csv, run_pipeline
from euclid_nisp_contamination_audit.exceptions import ProjectError
from euclid_nisp_contamination_audit.logging_utils import get_logger
from euclid_nisp_contamination_audit.masks import CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED
from euclid_nisp_contamination_audit.provenance import get_git_commit, sha256_config
from euclid_nisp_contamination_audit.results_io import Metric, write_summary
from euclid_nisp_contamination_audit.synthetic import build_synthetic_catalogue

LOGGER = get_logger(__name__)
REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_benchmark(path: Path, label: str, wall_time_s: float, peak_memory_mib: float, dataset_size: int) -> None:
    payload = {
        "label": label,
        "wall_time_seconds": wall_time_s,
        "peak_memory_mib": peak_memory_mib,
        "peak_memory_method": "tracemalloc (Python-level allocations, not full process RSS)",
        "dataset_size": dataset_size,
        "python_version": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "package_version": __version__,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    existing = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else []
    existing.append(payload)
    path.write_text(json.dumps(existing, indent=2), encoding="utf-8")


def _build_metrics(result, project_title: str) -> list[Metric]:
    metrics = [
        Metric(name="n_total", estimate=float(result.n_total), units="count", sample_size=result.n_total),
        Metric(name="n_clean", estimate=float(result.n_clean), units="count", sample_size=result.n_total),
        Metric(name="n_flagged", estimate=float(result.n_flagged), units="count", sample_size=result.n_total),
    ]
    for gm in result.continuum_snr_by_group:
        if gm.bootstrap is None:
            continue
        metrics.append(
            Metric(
                name=f"continuum_snr_{gm.group}",
                estimate=gm.bootstrap.estimate,
                uncertainty_low=gm.bootstrap.ci_low,
                uncertainty_high=gm.bootstrap.ci_high,
                units="dimensionless (S/N)",
                sample_size=gm.n,
            )
        )
    for group in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        if group in result.line_recovery_by_group:
            metrics.append(
                Metric(
                    name=f"line_recovery_rate_{group}",
                    estimate=result.line_recovery_by_group[group],
                    units="fraction",
                    sample_size=result.n_clean if group == CONTAMINATION_GROUP_CLEAN else result.n_flagged,
                )
            )
        if group in result.redshift_reliability_by_group:
            rel = result.redshift_reliability_by_group[group]
            metrics.append(
                Metric(
                    name=f"redshift_reliability_mean_{group}",
                    estimate=rel["mean"], units="dimensionless (spe_z_rel)", sample_size=rel["n"],
                )
            )
            metrics.append(
                Metric(
                    name=f"redshift_fraction_unreliable_{group}",
                    estimate=rel["fraction_unreliable"], units="fraction", sample_size=rel["n"],
                )
            )
    for row in result.threshold_sensitivity:
        metrics.append(
            Metric(
                name=f"threshold_sensitivity_{row.threshold_name}_flagged_cont_snr",
                estimate=row.flagged_mean_cont_snr, units="dimensionless (S/N)", sample_size=row.n_flagged,
            )
        )
    metrics.append(
        Metric(
            name="negative_control_observed_minus_null_mean_diff",
            estimate=result.negative_control.observed_difference - result.negative_control.null_mean,
            units="dimensionless (S/N)",
            sample_size=result.n_total,
        )
    )
    metrics.append(
        Metric(
            name="negative_control_p_value_two_sided",
            estimate=result.negative_control.p_value_two_sided,
            units="probability",
            sample_size=result.n_total,
        )
    )
    return metrics


def run_demo() -> None:
    config = load_config(REPO_ROOT / "config" / "analysis.yml")
    tracemalloc.start()
    start = time.perf_counter()

    catalogue = build_synthetic_catalogue(n_sources=300, seed=config.execution.seed)
    rows = [
        CatalogueRow(
            object_id=int(catalogue.object_id[i]), ra=float(catalogue.ra[i]), dec=float(catalogue.dec[i]),
            context_warning_flag=int(catalogue.context_warning_flag[i]),
            context_error_flag=int(catalogue.context_error_flag[i]),
            spe_cont_snr=float(catalogue.spe_cont_snr[i]), spe_z=float(catalogue.spe_z[i]),
            spe_z_err=float(catalogue.spe_z_err[i]), spe_z_rel=float(catalogue.spe_z_rel[i]),
            best_line_snr=float(catalogue.best_line_snr[i]),
        )
        for i in range(catalogue.object_id.size)
    ]
    result = run_pipeline(rows, config)

    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    out = REPO_ROOT / "results"
    out.mkdir(exist_ok=True)
    metrics = _build_metrics(result, "Euclid NISP Spectral-Contamination Reliability Audit (demo)")
    payload = write_summary(
        out / "summary.json",
        project="Euclid NISP Spectral-Contamination Reliability Audit (demo smoke test)",
        data_kind="synthetic_demo",
        metrics=metrics,
        provenance={
            "config_sha256": sha256_config(REPO_ROOT / "config" / "analysis.yml"),
            "git_commit": get_git_commit(REPO_ROOT),
            "package_version": __version__,
        },
        warnings=result.warnings,
    )
    (out / "warnings.json").write_text(json.dumps(result.warnings, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    _write_benchmark(out / "benchmarks.json", "demo", elapsed, peak / (1024 * 1024), len(rows))


def run_real_data(config_path: Path, catalogue_csv: Path, results_dir: Path) -> None:
    config = load_config(config_path)
    try:
        rows = read_catalogue_csv(catalogue_csv)
    except ProjectError as exc:
        raise SystemExit(
            f"Cannot run the real-data pipeline: {exc}. Run scripts/fetch_data.py "
            "(with explicit operator authorization) first."
        ) from exc

    tracemalloc.start()
    start = time.perf_counter()
    result = run_pipeline(rows, config)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    metrics = _build_metrics(result, config.project.title)
    provenance = {
        "config_sha256": sha256_config(config_path),
        "git_commit": get_git_commit(REPO_ROOT),
        "package_version": __version__,
        "n_sources": len(rows),
    }

    results_dir.mkdir(exist_ok=True)
    write_summary(
        results_dir / "summary.json",
        project=config.project.title,
        data_kind=config.input.data_mode,
        metrics=metrics,
        provenance=provenance,
        warnings=result.warnings,
    )
    (results_dir / "warnings.json").write_text(json.dumps(result.warnings, indent=2), encoding="utf-8")
    _write_benchmark(results_dir / "benchmarks.json", "real_data", elapsed, peak / (1024 * 1024), len(rows))
    print(f"Wrote {results_dir / 'summary.json'} ({len(metrics)} metrics, {len(result.warnings)} warnings)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", action="store_true", help="Run synthetic smoke data only")
    parser.add_argument("--config", type=Path, default=Path("config/analysis.yml"))
    parser.add_argument("--catalogue-csv", type=Path, default=Path("data/source_catalog.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    if args.demo:
        run_demo()
        return

    run_real_data(args.config, args.catalogue_csv, args.results_dir)


if __name__ == "__main__":
    main()

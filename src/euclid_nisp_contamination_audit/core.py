"""Pipeline orchestration composing the reusable scientific modules.

`run_pipeline` is the single non-notebook entry point for the catalogue-level
analysis (contamination group vs continuum S/N, line recovery, redshift
reliability). It is network-free — real/synthetic rows are passed in — so it
is exercised in tests against the synthetic catalogue generator.

The starter functions below (`Summary`, `validate_numeric`, `robust_summary`,
`demo_series`) are kept because `tests/test_starter_core.py` exercises them
and they remain useful small utilities used elsewhere in the pipeline.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from euclid_nisp_contamination_audit.config import AnalysisConfig
from euclid_nisp_contamination_audit.exceptions import (
    DataSchemaError,
    InsufficientDataError,
)
from euclid_nisp_contamination_audit.lines import DEFAULT_LINE_RECOVERY_SNR_THRESHOLD, line_recovery_rate
from euclid_nisp_contamination_audit.logging_utils import get_logger
from euclid_nisp_contamination_audit.masks import (
    CONTAMINATION_GROUP_CLEAN,
    CONTAMINATION_GROUP_FLAGGED,
    THRESHOLD_DEFINITIONS,
    contamination_group_labels,
)
from euclid_nisp_contamination_audit.redshift_metrics import reliability_summary_by_group
from euclid_nisp_contamination_audit.uncertainty import (
    BootstrapResult,
    PermutationTestResult,
    bootstrap_statistic,
    permuted_label_null_distribution,
)

LOGGER = get_logger(__name__)


@dataclass(frozen=True)
class Summary:
    count: int
    median: float
    mad: float


def validate_numeric(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError("values must be one-dimensional")
    if arr.size == 0:
        raise ValueError("values must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError("values contain non-finite entries")
    return arr


def robust_summary(values: np.ndarray) -> Summary:
    arr = validate_numeric(values)
    median = float(np.median(arr))
    mad = float(np.median(np.abs(arr - median)))
    return Summary(count=int(arr.size), median=median, mad=mad)


def demo_series(seed: int = 20260713, size: int = 128) -> np.ndarray:
    """Return deterministic synthetic data labelled only for smoke testing."""
    if size < 8:
        raise ValueError("size must be at least 8")
    rng = np.random.default_rng(seed)
    return rng.normal(loc=0.0, scale=1.0, size=size)


@dataclass(frozen=True)
class CatalogueRow:
    object_id: int
    ra: float
    dec: float
    context_warning_flag: int
    context_error_flag: int
    spe_cont_snr: float
    spe_z: float
    spe_z_err: float
    spe_z_rel: float
    best_line_snr: float
    spe_z_true: float | None = None  # only populated for synthetic validation rows


def read_catalogue_csv(path: str | Path) -> list[CatalogueRow]:
    """Read the joined catalogue CSV written by scripts/fetch_data.py.

    Contamination-indicator column note (real-data finding, see
    IMPLEMENTATION_PLAN.md / LOCAL_COMPLETION_REPORT.md): live TAP diagnostics
    on tile 102160061 showed `spe_context_warning_flag` is degenerate within
    the population that actually has an extracted spectrum -- it marks
    "spectrum was extracted" (bit 32768 set for exactly the ~15% of
    quality-table rows that have any galaxy-candidate/continuum-features row
    at all), not per-object contamination, so every row joined through
    spe_galaxy_candidates comes out 100% flagged (a degenerate split; see
    `masks.contamination_group_labels`'s warning path). The released
    `spe_error_flag` column, which retains real within-sample variance
    (~4.5% nonzero among objects with an extracted spectrum in a 5000-row
    diagnostic sample), is used instead as the operational "context warning"
    input here. `spe_context_error_flag` (the rarer, more severe tier) is
    kept as read for the alternate threshold-sensitivity definitions.
    """
    csv_path = Path(path)
    if not csv_path.is_file():
        raise InsufficientDataError(f"catalogue CSV not found: {csv_path}")
    rows: list[CatalogueRow] = []
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                CatalogueRow(
                    object_id=int(raw["object_id"]),
                    ra=float(raw["ra"]),
                    dec=float(raw["dec"]),
                    context_warning_flag=int(raw["spe_error_flag"]),
                    context_error_flag=int(raw["spe_context_error_flag"]),
                    spe_cont_snr=float(raw["spe_cont_snr"]),
                    spe_z=float(raw["spe_z"]),
                    spe_z_err=float(raw["spe_z_err"]),
                    spe_z_rel=float(raw["spe_z_rel"]),
                    best_line_snr=float(raw["best_line_snr"]) if raw["best_line_snr"] not in ("", None) else float("nan"),
                )
            )
    if not rows:
        raise InsufficientDataError(f"catalogue CSV {csv_path} has no rows")
    return rows


@dataclass(frozen=True)
class GroupMetric:
    group: str
    n: int
    bootstrap: BootstrapResult | None


@dataclass(frozen=True)
class ThresholdSensitivityRow:
    threshold_name: str
    threshold_description: str
    n_flagged: int
    n_clean: int
    flagged_mean_cont_snr: float
    clean_mean_cont_snr: float


@dataclass(frozen=True)
class PipelineResult:
    n_total: int
    n_clean: int
    n_flagged: int
    contamination_group: np.ndarray
    continuum_snr_by_group: list[GroupMetric]
    line_recovery_by_group: dict[str, float]
    redshift_reliability_by_group: dict[str, dict[str, float]]
    threshold_sensitivity: list[ThresholdSensitivityRow]
    negative_control: PermutationTestResult
    warnings: list[str] = field(default_factory=list)


def run_pipeline(
    rows: list[CatalogueRow],
    config: AnalysisConfig,
    line_snr_threshold: float = DEFAULT_LINE_RECOVERY_SNR_THRESHOLD,
) -> PipelineResult:
    """Run the contamination-vs-continuum-S/N / line-recovery / redshift-
    reliability audit over a joined catalogue.

    Raises InsufficientDataError immediately if `rows` is empty. Per-row
    failures are not expected at this stage (rows are already validated
    scalars from the catalogue CSV), but any InsufficientDataError,
    ConvergenceError or DataSchemaError raised by a per-group computation is
    caught and converted into a warning rather than aborting the whole run.
    """
    if not rows:
        raise InsufficientDataError("run_pipeline received an empty catalogue row list")

    all_warnings: list[str] = []
    n_total = len(rows)

    context_warning = np.array([r.context_warning_flag for r in rows])
    context_error = np.array([r.context_error_flag for r in rows])
    cont_snr = np.array([r.spe_cont_snr for r in rows], dtype=float)
    z_rel = np.array([r.spe_z_rel for r in rows], dtype=float)
    best_line_snr = np.array([r.best_line_snr for r in rows], dtype=float)

    group = contamination_group_labels(context_warning, context_error, threshold="context_warning_only")
    n_clean = int(np.sum(group == CONTAMINATION_GROUP_CLEAN))
    n_flagged = int(np.sum(group == CONTAMINATION_GROUP_FLAGGED))
    if n_clean == 0 or n_flagged == 0:
        all_warnings.append(
            f"contamination split is degenerate (n_clean={n_clean}, n_flagged={n_flagged}); "
            "group comparisons below may be unreliable or empty"
        )

    # 1. continuum S/N by group, with bootstrap CI
    continuum_snr_by_group: list[GroupMetric] = []
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        values = cont_snr[(group == grp) & np.isfinite(cont_snr)]
        if values.size == 0:
            all_warnings.append(f"no finite continuum S/N values in group '{grp}'")
            continuum_snr_by_group.append(GroupMetric(group=grp, n=0, bootstrap=None))
            continue
        try:
            boot = bootstrap_statistic(
                values, np.mean, config.validation.bootstrap_resamples, config.execution.seed,
                config.validation.confidence_level,
            )
        except InsufficientDataError as exc:
            all_warnings.append(f"continuum S/N bootstrap skipped for group '{grp}': {exc}")
            boot = None
        continuum_snr_by_group.append(GroupMetric(group=grp, n=int(values.size), bootstrap=boot))

    # 2. line recovery rate by group
    line_recovery_by_group: dict[str, float] = {}
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        values = best_line_snr[group == grp]
        try:
            line_recovery_by_group[grp] = line_recovery_rate(values, line_snr_threshold)
        except InsufficientDataError as exc:
            all_warnings.append(f"line recovery rate skipped for group '{grp}': {exc}")

    # 3. redshift reliability by group (released spe_z_rel; no external truth for real data)
    try:
        redshift_reliability_by_group = reliability_summary_by_group(z_rel, group)
    except InsufficientDataError as exc:
        all_warnings.append(f"redshift reliability summary skipped: {exc}")
        redshift_reliability_by_group = {}

    # 4. threshold sensitivity
    threshold_sensitivity: list[ThresholdSensitivityRow] = []
    for definition in THRESHOLD_DEFINITIONS:
        try:
            alt_group = contamination_group_labels(context_warning, context_error, threshold=definition.name)
        except DataSchemaError as exc:
            all_warnings.append(f"threshold '{definition.name}' skipped: {exc}")
            continue
        alt_flagged = alt_group == CONTAMINATION_GROUP_FLAGGED
        alt_clean = alt_group == CONTAMINATION_GROUP_CLEAN
        flagged_vals = cont_snr[alt_flagged & np.isfinite(cont_snr)]
        clean_vals = cont_snr[alt_clean & np.isfinite(cont_snr)]
        threshold_sensitivity.append(
            ThresholdSensitivityRow(
                threshold_name=definition.name,
                threshold_description=definition.description,
                n_flagged=int(np.sum(alt_flagged)),
                n_clean=int(np.sum(alt_clean)),
                flagged_mean_cont_snr=float(np.mean(flagged_vals)) if flagged_vals.size else float("nan"),
                clean_mean_cont_snr=float(np.mean(clean_vals)) if clean_vals.size else float("nan"),
            )
        )

    # 5. permuted-label negative control on continuum S/N
    finite_snr_mask = np.isfinite(cont_snr)
    try:
        negative_control = permuted_label_null_distribution(
            cont_snr[finite_snr_mask],
            group[finite_snr_mask],
            CONTAMINATION_GROUP_FLAGGED,
            CONTAMINATION_GROUP_CLEAN,
            n_permutations=config.validation.bootstrap_resamples,
            seed=config.execution.seed,
        )
    except InsufficientDataError as exc:
        all_warnings.append(f"permuted-label negative control skipped: {exc}")
        negative_control = PermutationTestResult(
            observed_difference=float("nan"), null_mean=float("nan"), null_std=float("nan"),
            p_value_two_sided=float("nan"), n_permutations=0,
        )

    return PipelineResult(
        n_total=n_total,
        n_clean=n_clean,
        n_flagged=n_flagged,
        contamination_group=group,
        continuum_snr_by_group=continuum_snr_by_group,
        line_recovery_by_group=line_recovery_by_group,
        redshift_reliability_by_group=redshift_reliability_by_group,
        threshold_sensitivity=threshold_sensitivity,
        negative_control=negative_control,
        warnings=all_warnings,
    )

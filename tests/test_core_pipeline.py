from __future__ import annotations

import pytest

from euclid_nisp_contamination_audit.core import CatalogueRow, run_pipeline
from euclid_nisp_contamination_audit.exceptions import InsufficientDataError
from euclid_nisp_contamination_audit.masks import CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED


def _rows_from_catalogue(catalogue) -> list[CatalogueRow]:
    rows = []
    for i in range(catalogue.object_id.size):
        rows.append(
            CatalogueRow(
                object_id=int(catalogue.object_id[i]),
                ra=float(catalogue.ra[i]),
                dec=float(catalogue.dec[i]),
                context_warning_flag=int(catalogue.context_warning_flag[i]),
                context_error_flag=int(catalogue.context_error_flag[i]),
                spe_cont_snr=float(catalogue.spe_cont_snr[i]),
                spe_z=float(catalogue.spe_z[i]),
                spe_z_err=float(catalogue.spe_z_err[i]),
                spe_z_rel=float(catalogue.spe_z_rel[i]),
                best_line_snr=float(catalogue.best_line_snr[i]),
            )
        )
    return rows


def test_run_pipeline_empty_raises(analysis_config) -> None:
    with pytest.raises(InsufficientDataError):
        run_pipeline([], analysis_config)


def test_run_pipeline_recovers_known_synthetic_effect(analysis_config, synthetic_catalogue) -> None:
    rows = _rows_from_catalogue(synthetic_catalogue)
    result = run_pipeline(rows, analysis_config)

    assert result.n_total == len(rows)
    assert result.n_clean > 0 and result.n_flagged > 0

    clean_metric = next(m for m in result.continuum_snr_by_group if m.group == CONTAMINATION_GROUP_CLEAN)
    flagged_metric = next(m for m in result.continuum_snr_by_group if m.group == CONTAMINATION_GROUP_FLAGGED)
    assert clean_metric.bootstrap is not None and flagged_metric.bootstrap is not None
    # the synthetic generator injects a known, strong effect: clean group has much higher continuum S/N
    assert clean_metric.bootstrap.estimate > flagged_metric.bootstrap.estimate

    assert result.line_recovery_by_group[CONTAMINATION_GROUP_CLEAN] > result.line_recovery_by_group[CONTAMINATION_GROUP_FLAGGED]

    assert result.redshift_reliability_by_group[CONTAMINATION_GROUP_CLEAN]["mean"] > (
        result.redshift_reliability_by_group[CONTAMINATION_GROUP_FLAGGED]["mean"]
    )

    assert len(result.threshold_sensitivity) == 3

    # the negative control must show the observed clean-vs-flagged effect is NOT consistent with noise
    assert result.negative_control.p_value_two_sided < 0.05


def test_run_pipeline_warns_on_degenerate_split(analysis_config) -> None:
    rows = [
        CatalogueRow(
            object_id=i, ra=60.0, dec=-48.0, context_warning_flag=0, context_error_flag=0,
            spe_cont_snr=10.0, spe_z=1.0, spe_z_err=0.01, spe_z_rel=0.9, best_line_snr=8.0,
        )
        for i in range(1, 21)
    ]
    result = run_pipeline(rows, analysis_config)
    assert result.n_flagged == 0
    assert any("degenerate" in w for w in result.warnings)

from __future__ import annotations

import json

import pytest

from euclid_nisp_contamination_audit.exceptions import DataSchemaError
from euclid_nisp_contamination_audit.results_io import Metric, validate_summary, write_summary


def test_write_summary_roundtrip(tmp_path) -> None:
    metrics = [Metric(name="n_sources", estimate=100.0, units="count", sample_size=100)]
    out_path = tmp_path / "summary.json"
    payload = write_summary(
        out_path, project="test", data_kind="synthetic_demo", metrics=metrics,
        provenance={"git_commit": "abc"}, warnings=["w1"],
    )
    on_disk = json.loads(out_path.read_text(encoding="utf-8"))
    assert on_disk == payload
    assert on_disk["metrics"][0]["name"] == "n_sources"


def test_validate_summary_missing_key() -> None:
    with pytest.raises(DataSchemaError):
        validate_summary({"project": "x"})


def test_validate_summary_bad_metric() -> None:
    payload = {
        "project": "x", "data_kind": "synthetic_demo", "metrics": [{"name": "a"}],
        "provenance": {}, "warnings": [],
    }
    with pytest.raises(DataSchemaError):
        validate_summary(payload)

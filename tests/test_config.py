from __future__ import annotations

import pytest

from euclid_nisp_contamination_audit.exceptions import DataSchemaError
from euclid_nisp_contamination_audit.config import load_config


def test_load_config_real_file() -> None:
    cfg = load_config("config/analysis.yml")
    assert cfg.project.repository == "euclid-nisp-contamination-audit"
    assert cfg.execution.seed == 20260713
    assert 0 < cfg.validation.confidence_level < 1


def test_load_config_missing_file() -> None:
    with pytest.raises(DataSchemaError):
        load_config("config/does_not_exist.yml")


def test_load_config_missing_section(tmp_path) -> None:
    bad = tmp_path / "bad.yml"
    bad.write_text("project:\n  title: x\n", encoding="utf-8")
    with pytest.raises(DataSchemaError):
        load_config(bad)

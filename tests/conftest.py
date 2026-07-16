from __future__ import annotations

import pytest

from euclid_nisp_contamination_audit.config import load_config
from euclid_nisp_contamination_audit.synthetic import build_synthetic_catalogue, build_synthetic_spectrum

CONFIG_PATH = "config/analysis.yml"


@pytest.fixture()
def analysis_config():
    return load_config(CONFIG_PATH)


@pytest.fixture()
def synthetic_catalogue():
    return build_synthetic_catalogue(n_sources=300, seed=20260713)


@pytest.fixture()
def synthetic_clean_spectrum():
    return build_synthetic_spectrum(object_id=1, seed=1, contaminated=False, true_z=1.0)


@pytest.fixture()
def synthetic_contaminated_spectrum():
    return build_synthetic_spectrum(object_id=2, seed=2, contaminated=True, true_z=1.0)

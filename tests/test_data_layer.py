from __future__ import annotations

import numpy as np
import pytest
from astropy.io.votable import from_table, writeto
from astropy.table import Table

from euclid_nisp_contamination_audit.core import read_catalogue_csv
from euclid_nisp_contamination_audit.exceptions import DataSchemaError, InsufficientDataError
from euclid_nisp_contamination_audit.io import read_spectrum_votable


def _write_votable(path, colnames, n=50) -> None:
    data = {name: np.linspace(1, 2, n) for name in colnames}
    table = Table(data)
    vot = from_table(table)
    writeto(vot, str(path))


def test_read_spectrum_votable_success(tmp_path) -> None:
    path = tmp_path / "spec.vot"
    _write_votable(path, ["WAVELENGTH", "SIGNAL", "UNCERTAINTY", "MASK", "QUALITY"])
    spectrum = read_spectrum_votable(path, object_id=42)
    assert spectrum.object_id == 42
    assert spectrum.wavelength_angstrom.size == 50


def test_read_spectrum_votable_missing_file(tmp_path) -> None:
    with pytest.raises(DataSchemaError):
        read_spectrum_votable(tmp_path / "nope.vot", object_id=1)


def test_read_spectrum_votable_missing_columns(tmp_path) -> None:
    path = tmp_path / "bad.vot"
    _write_votable(path, ["WAVELENGTH", "SIGNAL"])
    with pytest.raises(DataSchemaError):
        read_spectrum_votable(path, object_id=1)


def test_read_catalogue_csv_missing_file(tmp_path) -> None:
    with pytest.raises(InsufficientDataError):
        read_catalogue_csv(tmp_path / "nope.csv")


def test_read_catalogue_csv_roundtrip(tmp_path) -> None:
    path = tmp_path / "cat.csv"
    path.write_text(
        "object_id,ra,dec,spe_error_flag,spe_context_error_flag,spe_cont_snr,spe_z,spe_z_err,spe_z_rel,best_line_snr\n"
        "1,60.0,-48.0,0,0,10.0,1.0,0.01,0.9,8.0\n",
        encoding="utf-8",
    )
    rows = read_catalogue_csv(path)
    assert len(rows) == 1
    assert rows[0].object_id == 1

from __future__ import annotations

import pytest

from euclid_nisp_contamination_audit.exceptions import ProvenanceError
from euclid_nisp_contamination_audit.provenance import (
    ManifestRow,
    append_manifest_row,
    get_git_commit,
    read_manifest,
    sha256_bytes,
    sha256_file,
)


def _row(product_id: str = "obj_1") -> ManifestRow:
    return ManifestRow(
        product_id=product_id, source="IRSA/Euclid-Q1", source_url="https://irsa.ipac.caltech.edu",
        retrieved_utc="2026-07-14T00:00:00+00:00", sha256="deadbeef", file_size_bytes=123,
        selection_reason="test", licence_or_terms="test terms",
    )


def test_manifest_roundtrip(tmp_path) -> None:
    manifest_path = tmp_path / "manifest.csv"
    append_manifest_row(manifest_path, _row("a"))
    append_manifest_row(manifest_path, _row("b"))
    rows = read_manifest(manifest_path)
    assert len(rows) == 2
    assert rows[0]["product_id"] == "a"


def test_read_manifest_missing_file(tmp_path) -> None:
    with pytest.raises(ProvenanceError):
        read_manifest(tmp_path / "nope.csv")


def test_sha256_file_matches_bytes(tmp_path) -> None:
    p = tmp_path / "f.bin"
    p.write_bytes(b"hello world")
    assert sha256_file(p) == sha256_bytes(b"hello world")


def test_get_git_commit_never_raises(tmp_path) -> None:
    result = get_git_commit(tmp_path)
    assert isinstance(result, str)
    assert result  # either a commit hash or the LOCAL_UNCOMMITTED sentinel

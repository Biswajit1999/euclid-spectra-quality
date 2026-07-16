"""IRSA/Euclid Q1 NISP catalogue query wrappers and VOTable spectrum I/O.

Real network access (`query_joined_catalogue`, `fetch_spectrum_votable`) is
isolated in this module so the rest of the package stays testable offline.
All TAP queries here are restricted by `tileid` and/or an explicit
`object_id` list, per docs/DATASET_PLAN.md and the hang observed in the
sibling astrometry project from an unrestricted query.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from euclid_nisp_contamination_audit.exceptions import ArchiveAccessError, DataSchemaError

CATALOGUE_COLUMNS = (
    "object_id",
    "ra",
    "dec",
    "tileid",
    "spe_context_error_flag",
    "spe_context_warning_flag",
    "spe_error_flag",
    "spe_warning_flag",
    "spe_grism",
    "spe_npix",
    "spe_z",
    "spe_z_err",
    "spe_z_rel",
    "spe_cont_snr",
    "spe_subclass",
    "best_line_snr",
)


def _masked_scalar(value):
    """Return None for a masked/NULL TAP result cell, else the raw value."""
    return None if np.ma.is_masked(value) else value


def query_joined_catalogue(tileid: int, n_sources: int, grism: str = "RGS"):
    """Deterministic (ORDER BY object_id), tile-restricted join of the three
    Euclid Q1 NISP spectroscopy catalogue tables needed for this audit:
    quality flags (contamination proxy), galaxy redshift candidates, and MER
    RA/Dec. Returns an astropy Table. Line-recovery SNR is queried separately
    per selected object_id (see `query_best_line_snr`) because joining the
    line-features table (many rows per object) inside the same query is slow
    without an object_id restriction.
    """
    try:
        from astroquery.ipac.irsa import Irsa
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ArchiveAccessError("astroquery is not installed in this environment") from exc

    query = f"""
        SELECT TOP {n_sources} q.object_id, m.ra, m.dec, m.tileid,
               q.spe_context_error_flag, q.spe_context_warning_flag,
               q.spe_error_flag, q.spe_warning_flag, q.spe_grism, q.spe_npix,
               g.spe_z, g.spe_z_err, g.spe_z_rel, g.spe_cont_snr, g.spe_subclass
        FROM euclid_q1_spectro_zcatalog_spe_quality q
        JOIN euclid_q1_spectro_zcatalog_spe_galaxy_candidates g
            ON q.object_id = g.object_id AND g.spe_rank = 0
        JOIN euclid_q1_mer_catalogue m ON q.object_id = m.object_id
        WHERE q.spe_grism = '{grism}' AND m.tileid = {int(tileid)}
        ORDER BY q.object_id
    """
    try:
        result = Irsa.query_tap(query)
        table = result.to_table()
    except Exception as exc:  # noqa: BLE001
        raise ArchiveAccessError(f"IRSA TAP joined-catalogue query failed: {exc}") from exc
    if len(table) == 0:
        raise ArchiveAccessError(
            f"IRSA TAP query returned zero rows for tileid={tileid}, grism={grism}"
        )
    return table


_MAX_IN_LIST_SIZE = 900  # IRSA's Oracle TAP backend rejects IN-lists over 1000 expressions (confirmed live)


def query_best_line_snr(object_ids: list[int]) -> dict[int, float]:
    """Best (max) Gaussian-fit line S/N per object_id, restricted to an
    explicit object_id list (never an unrestricted scan of the line-features
    table, which has many rows per object across the whole survey).

    Batches the IN-list at `_MAX_IN_LIST_SIZE`: a live query with 6000
    object_ids in a single IN-clause failed with
    "ORA-01795: maximum number of expressions in a list is 1000" from IRSA's
    Oracle TAP backend.
    """
    try:
        from astroquery.ipac.irsa import Irsa
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ArchiveAccessError("astroquery is not installed in this environment") from exc
    if not object_ids:
        raise DataSchemaError("object_ids must not be empty")

    best: dict[int, float] = {}
    for start in range(0, len(object_ids), _MAX_IN_LIST_SIZE):
        batch = object_ids[start : start + _MAX_IN_LIST_SIZE]
        id_list = ",".join(str(int(o)) for o in batch)
        query = f"""
            SELECT object_id, spe_line_snr_gf
            FROM euclid_q1_spe_lines_line_features
            WHERE object_id IN ({id_list}) AND spe_line_snr_gf IS NOT NULL
        """
        try:
            result = Irsa.query_tap(query)
            table = result.to_table()
        except Exception as exc:  # noqa: BLE001
            raise ArchiveAccessError(f"IRSA TAP line-features query failed: {exc}") from exc
        best.update(_reduce_best_line_snr(table))
    return best


def _reduce_best_line_snr(table) -> dict[int, float]:
    best: dict[int, float] = {}
    for row in table:
        oid = int(row["object_id"])
        snr = float(row["spe_line_snr_gf"])
        if not np.isfinite(snr):
            continue
        if oid not in best or snr > best[oid]:
            best[oid] = snr
    return best


def find_spectrum_download_url(tileid: int, object_id: int) -> tuple[str, str]:
    """Resolve the IRSA `spectrumdm/convert` URL for one object's combined
    1-D SIR spectrum via the object_id-to-spectral-file association table.
    Returns (full_url, obs_publisher_did).
    """
    try:
        from astroquery.ipac.irsa import Irsa
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ArchiveAccessError("astroquery is not installed in this environment") from exc

    query = f"""
        SELECT path, obs_publisher_did
        FROM euclid.objectid_spectrafile_association_q1
        WHERE objectid = {int(object_id)} AND tileid = {int(tileid)} AND path IS NOT NULL
    """
    try:
        result = Irsa.query_tap(query)
        table = result.to_table()
    except Exception as exc:  # noqa: BLE001
        raise ArchiveAccessError(f"IRSA TAP spectrum-association query failed: {exc}") from exc
    if len(table) == 0:
        raise ArchiveAccessError(f"no spectral file association found for object_id={object_id}")

    path = str(table[0]["path"])
    did = str(table[0]["obs_publisher_did"])
    return f"https://irsa.ipac.caltech.edu/{path}", did


@dataclass(frozen=True)
class Spectrum1D:
    object_id: int
    wavelength_angstrom: np.ndarray
    signal: np.ndarray
    uncertainty: np.ndarray
    mask: np.ndarray
    quality: np.ndarray


def download_spectrum_votable(url: str, out_path: str | Path, timeout_seconds: int = 60) -> bytes:
    """Download the raw VOTable bytes for one spectrum. Real network I/O,
    isolated here so the parser (`read_spectrum_votable`) stays offline-testable.
    """
    import requests

    try:
        response = requests.get(url, timeout=timeout_seconds)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise ArchiveAccessError(f"spectrum download failed for {url}: {exc}") from exc
    content = response.content
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(content)
    return content


def read_spectrum_votable(path: str | Path, object_id: int) -> Spectrum1D:
    """Parse a locally-saved IRSA `spectrumdm/convert` VOTable spectrum into a
    Spectrum1D. Raises DataSchemaError if expected columns are missing.
    """
    from astropy.io.votable import parse

    path = Path(path)
    if not path.is_file():
        raise DataSchemaError(f"spectrum file not found: {path}")
    try:
        votable = parse(str(path))
        table = votable.get_first_table().to_table()
    except Exception as exc:  # noqa: BLE001
        raise DataSchemaError(f"failed to parse spectrum VOTable {path}: {exc}") from exc

    required = ("WAVELENGTH", "SIGNAL", "UNCERTAINTY", "MASK", "QUALITY")
    missing = [c for c in required if c not in table.colnames]
    if missing:
        raise DataSchemaError(f"spectrum VOTable {path} missing columns: {missing}")

    wavelength = np.asarray(table["WAVELENGTH"], dtype=float)
    signal = np.asarray(table["SIGNAL"], dtype=float)
    uncertainty = np.asarray(table["UNCERTAINTY"], dtype=float)
    mask = np.asarray(table["MASK"], dtype=int)
    quality = np.asarray(table["QUALITY"], dtype=int)

    if not (wavelength.shape == signal.shape == uncertainty.shape == mask.shape == quality.shape):
        raise DataSchemaError(f"inconsistent array lengths in spectrum VOTable {path}")

    return Spectrum1D(
        object_id=object_id,
        wavelength_angstrom=wavelength,
        signal=signal,
        uncertainty=uncertainty,
        mask=mask,
        quality=quality,
    )


__all__ = [
    "CATALOGUE_COLUMNS",
    "Spectrum1D",
    "download_spectrum_votable",
    "find_spectrum_download_url",
    "query_best_line_snr",
    "query_joined_catalogue",
    "read_spectrum_votable",
]

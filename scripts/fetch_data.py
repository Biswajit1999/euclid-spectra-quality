"""Deterministic, provenance-recording fetch of real public Euclid Q1 NISP
spectroscopy products from IRSA's Euclid Q1 mirror.

`astroquery.esa.euclid` is unavailable in pinned astroquery==0.4.7 (confirmed
ImportError in the euclid-q1-vis-psf-astrometry-audit sibling project), so
this goes through `astroquery.ipac.irsa.Irsa`'s TAP service, exactly as that
sibling. All TAP queries are restricted to a single Q1 tile (see
IMPLEMENTATION_PLAN.md Sec 2) to avoid the multi-minute unrestricted-query
hang observed there.

This script performs real network requests and must only be invoked with
explicit user authorization for the session, per docs/DATASET_PLAN.md.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from euclid_nisp_contamination_audit.exceptions import ArchiveAccessError
from euclid_nisp_contamination_audit.io import (
    CATALOGUE_COLUMNS,
    download_spectrum_votable,
    find_spectrum_download_url,
    query_best_line_snr,
    query_joined_catalogue,
)
from euclid_nisp_contamination_audit.logging_utils import get_logger
from euclid_nisp_contamination_audit.provenance import ManifestRow, append_manifest_row, sha256_file

LOGGER = get_logger(__name__)

# Single Q1 tile, RGS grism, confirmed live to hold 17414 combined spectra
# with a non-null spectral-file path (see IMPLEMENTATION_PLAN.md Sec 2).
TILE_ID = 102160061
GRISM = "RGS"
SOURCE_URL = "https://irsa.ipac.caltech.edu"
LICENCE_TERMS = (
    "Euclid Q1 public data release (ESA Euclid Consortium), mirrored by IRSA/IPAC; "
    "standard IRSA archive usage terms apply, https://irsa.ipac.caltech.edu/holdings/copyright.html"
)


def _masked_scalar(value, default=None):
    return default if np.ma.is_masked(value) else value


def _write_catalogue_csv(table, best_line_snr: dict[int, float], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CATALOGUE_COLUMNS))
        writer.writeheader()
        for row in table:
            object_id = int(row["object_id"])
            writer.writerow(
                {
                    "object_id": object_id,
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "tileid": int(row["tileid"]),
                    "spe_context_error_flag": int(row["spe_context_error_flag"]),
                    "spe_context_warning_flag": int(row["spe_context_warning_flag"]),
                    "spe_error_flag": int(row["spe_error_flag"]),
                    "spe_warning_flag": int(row["spe_warning_flag"]),
                    "spe_grism": str(row["spe_grism"]),
                    "spe_npix": int(_masked_scalar(row["spe_npix"], 0)),
                    "spe_z": float(_masked_scalar(row["spe_z"], np.nan)),
                    "spe_z_err": float(_masked_scalar(row["spe_z_err"], np.nan)),
                    "spe_z_rel": float(_masked_scalar(row["spe_z_rel"], np.nan)),
                    "spe_cont_snr": float(_masked_scalar(row["spe_cont_snr"], np.nan)),
                    "spe_subclass": str(_masked_scalar(row["spe_subclass"], "")),
                    "best_line_snr": best_line_snr.get(object_id, np.nan),
                }
            )
            n_written += 1
    return n_written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-sources", type=int, default=800, help="Deterministic catalogue sample size.")
    parser.add_argument("--n-spectra", type=int, default=40, help="Number of real 1-D spectra to download.")
    parser.add_argument("--tileid", type=int, default=TILE_ID)
    parser.add_argument("--out-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifest.csv"))
    parser.add_argument("--catalogue-csv", type=Path, default=Path("data/source_catalog.csv"))
    parser.add_argument(
        "--i-have-authorization",
        action="store_true",
        help=(
            "Required flag confirming the operator has explicitly authorized this "
            "real network download in the current session."
        ),
    )
    args = parser.parse_args()

    if not args.i_have_authorization:
        raise SystemExit(
            "Refusing to download real archive data without --i-have-authorization. "
            "This flag exists so the download only runs after the operator has "
            "explicitly confirmed it in the current session (see docs/DATASET_PLAN.md)."
        )

    LOGGER.info("Querying IRSA TAP joined catalogue for tileid=%d, grism=%s", args.tileid, GRISM)
    table = query_joined_catalogue(args.tileid, args.n_sources, grism=GRISM)
    object_ids = [int(row["object_id"]) for row in table]
    LOGGER.info("Selected %d catalogue rows; querying best line S/N per object", len(object_ids))
    best_line_snr = query_best_line_snr(object_ids)

    n_catalogue_rows = _write_catalogue_csv(table, best_line_snr, args.catalogue_csv)
    LOGGER.info("Wrote %d catalogue rows to %s", n_catalogue_rows, args.catalogue_csv)

    retrieved_utc = datetime.now(timezone.utc).isoformat()
    catalogue_sha = sha256_file(args.catalogue_csv)
    catalogue_size = args.catalogue_csv.stat().st_size
    append_manifest_row(
        args.manifest,
        ManifestRow(
            product_id=f"euclid_q1_nisp_catalogue_tile{args.tileid}",
            source="IRSA/Euclid-Q1",
            source_url=SOURCE_URL,
            retrieved_utc=retrieved_utc,
            sha256=catalogue_sha,
            file_size_bytes=catalogue_size,
            selection_reason=(
                f"deterministic TOP {args.n_sources} join of spe_quality + spe_galaxy_candidates + "
                f"mer_catalogue, ORDER BY object_id, grism={GRISM}, tileid={args.tileid}"
            ),
            licence_or_terms=LICENCE_TERMS,
        ),
    )

    # Deterministic (sorted) subset of real 1-D spectrum downloads, needed for
    # the low/high-contamination-spectrum figure and the line-injection
    # validation run alongside real data.
    n_downloaded = 0
    for object_id in sorted(object_ids)[: args.n_spectra]:
        product_id = f"euclid_q1_nisp_spectrum_{object_id}"
        out_path = args.out_dir / f"{product_id}.vot"
        try:
            url, did = find_spectrum_download_url(args.tileid, object_id)
            download_spectrum_votable(url, out_path)
        except ArchiveAccessError as exc:
            LOGGER.warning("Skipping object_id=%d spectrum download: %s", object_id, exc)
            continue

        digest = sha256_file(out_path)
        size = out_path.stat().st_size
        append_manifest_row(
            args.manifest,
            ManifestRow(
                product_id=product_id,
                source="IRSA/Euclid-Q1",
                source_url=url,
                retrieved_utc=retrieved_utc,
                sha256=digest,
                file_size_bytes=size,
                selection_reason=(
                    f"deterministic sorted subset (first {args.n_spectra} object_id values) of tile "
                    f"{args.tileid} catalogue sample; obs_publisher_did={did}"
                ),
                licence_or_terms=LICENCE_TERMS,
            ),
        )
        n_downloaded += 1
        LOGGER.info("Downloaded spectrum for object_id=%d (%d bytes)", object_id, size)

    if n_downloaded == 0:
        raise ArchiveAccessError("no spectra were successfully downloaded")

    print(
        f"Wrote {n_catalogue_rows} catalogue rows to {args.catalogue_csv}; "
        f"downloaded {n_downloaded} real spectra to {args.out_dir}"
    )


if __name__ == "__main__":
    main()

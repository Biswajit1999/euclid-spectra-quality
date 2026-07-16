# Implementation Plan — Euclid NISP Spectral-Contamination Reliability Audit

## 0. Prior-work check

No prior partial implementation existed (`src/` held only stub `NotImplementedError` modules,
no `IMPLEMENTATION_PLAN.md`). Starting fresh. A dedicated conda env
`euclid-nisp-contamination-audit` (python 3.11) already existed with all pinned
dependencies installed — reused as-is.

## 1. Literature verification (docs/LITERATURE_SEEDS.md)

Verified live via WebFetch against arXiv abstract pages:

| Seed | Verified title | arXiv | Status |
|---|---|---|---|
| Euclid Q1 overview | "Euclid Quick Data Release (Q1) — Data release overview", Euclid Collaboration: H. Aussel et al. | 2503.15302 | VERIFIED |
| Euclid SIR spectroscopic processing | "Euclid Quick Data Release (Q1): From spectrograms to spectra: the SIR spectroscopic Processing Function", Euclid Collaboration: Y. Copin et al. | 2503.15307 | VERIFIED — abstract explicitly states the SIR function "subtracts cross-contaminations, minimizes self-contamination" |
| Euclid NIR processing | "Euclid Quick Data Release (Q1). NIR processing and data products", Euclid Collaboration: G. Polenta et al. | 2503.15304 | VERIFIED |
| "Euclid Q1 spectroscopy science papers" | Generic seed, not a single citable paper | — | TODO_VERIFY — left as a general pointer in references.bib, not cited as a concrete work |
| Astropy/specutils official packages | astropy.org, specutils.readthedocs.io | — | VERIFIED via project software citation (not a paper) |

## 2. Real-data access plan (verified live against IRSA TAP, not assumed)

`astroquery.esa.euclid` is unavailable in pinned astroquery==0.4.7 (confirmed by ImportError in
sibling project). Used `astroquery.ipac.irsa.Irsa` instead, exactly as the
euclid-q1-vis-psf-astrometry-audit sibling. Discovered live (do not assume names):

- `Irsa.list_catalogs()` lists Euclid Q1 NISP spectroscopy tables. Relevant ones used here:
  - `euclid_q1_spectro_zcatalog_spe_quality` — per-object processing quality flags, incl.
    `spe_context_error_flag` / `spe_context_warning_flag` (the SIR pipeline's "context" flags,
    the released proxy for spectral-trace contamination — the SIR paper (2503.15307) confirms the
    pipeline explicitly performs cross-contamination subtraction; the exact bit semantics of the
    context flag are not published in the queryable TAP column metadata, so this project's use of
    "context flag != 0" as the contamination indicator is documented as an explicit, named
    assumption, not an internal Euclid definition — see docs/ASSUMPTIONS_AND_LIMITATIONS.md).
  - `euclid_q1_spectro_zcatalog_spe_galaxy_candidates` — `spe_z`, `spe_z_err`, `spe_z_rel`,
    `spe_cont_snr` (continuum S/N), `spe_subclass`.
  - `euclid_q1_spe_lines_line_features` — per-line `spe_line_snr_gf`, `spe_line_flux_gf`,
    `spe_line_qual_gf`, `spe_line_aon` (amplitude-over-noise).
  - `euclid.objectid_spectrafile_association_q1` — `path`/`hdu` needed to build a download URL
    for the actual combined 1-D NISP spectrum (`WAVELENGTH, SIGNAL, MASK, QUALITY, VAR, NDITH,
    UNCERTAINTY` columns), returned as a VOTable by IRSA's `spectrumdm/convert` API, confirmed
    live with a real 200-OK download (108 KB, 531-pixel spectrum) at
    `https://irsa.ipac.caltech.edu/{path}` for `objectid=2721023473674756138`, `tileid=102160061`.
- All TAP queries are restricted with `WHERE tileid = 102160061` (single Q1 tile, RGS grism,
  17414 spectra with a non-null file path) to avoid the 5+ minute unrestricted-query hang observed
  in the sibling project. A tile-count query confirmed this live before any bulk query was run.
- Real per-object spectrum downloads (needed for line-injection-into-real-continua and the
  low/high-contamination spectrum figure) are restricted to a small, deterministic
  (sorted-by-object_id) subset of `--n-spectra` objects (default 40), not the full 17k-row tile.
- Licence/terms: standard IRSA archive usage terms
  (https://irsa.ipac.caltech.edu/holdings/copyright.html); Euclid Q1 is a public ESA data release,
  confirmed PUBLIC (no proprietary/embargo restriction) via the IRSA Euclid Q1 holdings page.

## 3. Module plan (mirrors docs/RESEARCH_BLUEPRINT.md "Reusable scientific modules")

- `config.py`, `exceptions.py`, `logging_utils.py`, `provenance.py`, `results_io.py` — ported
  near-verbatim from euclid-q1-vis-psf-astrometry-audit, renamed to this package.
- `io.py` — TAP catalog query wrappers + VOTable spectrum reader/writer (local cache under
  `data/raw/`).
- `masks.py` — contamination-group labelling (`context_flag != 0`), finite-value masking,
  quality-flag masking.
- `continuum.py` — robust local continuum estimate (median/polynomial in a rest-frame-agnostic
  observed-wavelength window) and continuum S/N.
- `lines.py` — synthetic Gaussian emission-line injection into a real continuum + recovery via
  simple matched-filter/Gaussian fit; also exposes a wrapper around the released
  `spe_line_snr_gf` recovery flag for the real-catalogue line-recovery metric.
- `redshift_metrics.py` — `delta_z/(1+z)`, outlier fraction (|Δz|/(1+z) > 0.15, the standard
  photo-z/spec-z outlier convention), reliability-vs-contamination grouping.
- `uncertainty.py` — `bootstrap_statistic` (1000 resamples, seed 20260713) and
  `check_fit_convergence` (covariance condition number + reduced chi-square), kept separate.
- `synthetic.py` — shared synthetic spectrum + catalogue generator for tests/`--demo`.
- `core.py` — `run_pipeline` orchestrator over the joined catalogue rows, per-row
  try/except(InsufficientDataError, ConvergenceError, DataSchemaError) → warning, never a hard
  abort; raises `InsufficientDataError` immediately on an empty input.
- `plotting.py` — shared low-level plotting helpers used by `scripts/make_figures.py`.

## 4. Validation gate (docs/VALIDATION_CONTRACT.md) — must pass before real data is touched

1. Line injection into real continua (inject known Gaussian flux into an observed continuum +
   noise, recover with a fit, report recovered/injected flux ratio and detection completeness by
   S/N bin).
2. Bootstrap contamination groups (bootstrap CI on group-mean continuum S/N, line-recovery rate,
   redshift outlier fraction).
3. Redshift outlier fraction with bootstrap CI, split by contamination group.
4. Threshold sensitivity (recompute the clean/flagged split and headline metrics for several
   candidate flag thresholds; report stability).
5. Permuted-label negative control (shuffle contamination-group labels and show the group
   difference collapses to within its null distribution).

Run first on the synthetic generator output; only proceed to `scripts/fetch_data.py` once this
passes.

## 5. Figures (docs/FIGURE_AND_UI_SPEC.md)

1. `low_high_contamination_spectra` — example spectra, one clean / one flagged.
2. `snr_vs_contamination` — continuum S/N distribution by contamination group.
3. `redshift_error` — Δz/(1+z) distribution by contamination group.
4. `outlier_fraction` — outlier fraction with bootstrap CI by contamination group and threshold.
5. `line_recovery` — injected-vs-recovered line flux / completeness curve.
6. `field_map` — RA/Dec scatter of the selected sample coloured by contamination group.

## 6. Order of work

Foundation → data layer (fetch_data.py + synthetic.py) → scientific modules → validation gate on
synthetic data → benchmarks → demo figures/report skeleton → React dashboard (eslint/recharts fix
first) → real data fetch (`--i-have-authorization`) → real pipeline run → regenerate figures on
real data → update report Results/Limitations with real numbers → final pytest/ruff/mypy/npm
build → `LOCAL_COMPLETION_REPORT.md` / `_PROJECT_LOG.md`.

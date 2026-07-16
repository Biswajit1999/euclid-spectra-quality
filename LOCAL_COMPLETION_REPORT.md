# Local Completion Report — Euclid NISP Spectral-Contamination Reliability Audit

Author: Biswajit Jana. This report documents a local implementation pass (BUILD_FIRST
priority 9.4/10). No git operations were performed. Nothing has been published.

## 1. Environment

- Reused the pre-existing `euclid-nisp-contamination-audit` conda environment (Python 3.11.15),
  matching this project's `pyproject.toml` pins exactly (numpy==1.26.4, scipy==1.13.1,
  pandas==2.2.2, matplotlib==3.9.0, pyyaml==6.0.1, astropy==6.1.0, specutils==1.20.3,
  astroquery==0.4.7, requests==2.34.2; dev: pytest==8.2.2, pytest-cov==5.0.0, ruff==0.5.5,
  mypy==1.10.1, types-PyYAML, types-requests).
- No local LaTeX toolchain in this environment; `reports/report.tex` was verified for structural
  completeness only, not compiled to PDF.

## 2. Files created or changed

Foundation: `config.py`, `exceptions.py`, `logging_utils.py`, `provenance.py`, `results_io.py`
(ported from euclid-q1-vis-psf-astrometry-audit, renamed). Data layer: `io.py` (IRSA TAP
query/download + VOTable spectrum parsing), `synthetic.py`, `scripts/fetch_data.py`,
`scripts/sync_web_assets.py`. Scientific modules: `masks.py`, `continuum.py`, `lines.py`,
`redshift_metrics.py`, `uncertainty.py`, `core.py` (`run_pipeline`), `plotting.py`. 13 test files
(`tests/test_*.py` + `conftest.py`). Figures/report: `scripts/make_figures.py`, `reports/report.tex`,
`reports/references.bib`. Web dashboard: `web-react/src/App.jsx` rewritten,
`web-react/eslint.config.js` fixed (`react/jsx-uses-vars`/`react/jsx-uses-react`), `recharts`
removed from `package.json`, `web-react/public/project.json` rewritten. `IMPLEMENTATION_PLAN.md`,
`_PROJECT_LOG.md`, this report.

## 3. Exact commands run (in order)

```bash
python -m pip install -e ".[dev]"
pytest -q                                  # 53 passed
ruff check src tests scripts               # All checks passed
mypy src                                   # Success: no issues found in 15 source files
python scripts/run_analysis.py --demo
python scripts/make_figures.py --demo
cd web-react && npm install && npm run lint && npm run build
python scripts/sync_web_assets.py
# Real-data pipeline, run only after explicit operator authorization in chat:
python scripts/fetch_data.py --i-have-authorization --n-sources 800 --n-spectra 40
python scripts/run_analysis.py             # first run: exposed a real column-name bug (see Sec 4)
# bug fixed, then:
python scripts/fetch_data.py --i-have-authorization --n-sources 6000 --n-spectra 3   # broader real sample
python scripts/run_analysis.py
python scripts/make_figures.py
python scripts/sync_web_assets.py
pytest -q && ruff check src tests scripts && mypy src   # re-verified clean
cd web-react && npm run lint && npm run build            # re-verified clean
```

## 4. Test / lint / build results

- **pytest**: 53 tests passed, 0 failed (both before and after the real-data bug fixes below).
- **ruff**: clean on `src tests scripts`.
- **mypy**: clean on `src` (0 errors, 15 source files).
- **web-react**: `npm run lint` and `npm run build` both clean.

### Bugs found and fixed during implementation
1. `core.read_catalogue_csv` read unprefixed column names (`context_warning_flag`) but
   `fetch_data.py`'s CSV writer used the real, `spe_`-prefixed IRSA column names
   (`spe_context_warning_flag`) — the first real-data run crashed with `KeyError`. Fixed by
   aligning the reader to the actual written schema.
2. `io.query_best_line_snr` built a single SQL `IN (...)` clause for all requested object_ids;
   a live query with 6000 IDs failed with IRSA's Oracle TAP backend error `ORA-01795: maximum
   number of expressions in a list is 1000`. Fixed by batching the IN-list at 900 IDs per query.
3. One spectrum download (`object_id=2719168617674311044`) hit a transient `502 Bad Gateway` from
   IRSA; `fetch_data.py`'s per-object try/except caught it, logged a warning, and continued rather
   than aborting the whole 40-spectrum run — worked exactly as designed.
4. **Real, substantive finding (not a code bug)**: `spe_context_warning_flag`, this project's
   originally assumed contamination indicator, turned out to be structurally degenerate for the
   population actually usable for this audit. Live TAP diagnostics showed it is set for exactly
   the subset of quality-table rows that have *any* extracted-spectrum measurement at all (i.e. it
   marks "spectrum was extracted", correlating almost perfectly with grism=RGS), not a per-object
   contamination signature — every object joined through `spe_galaxy_candidates` (required for
   `spe_cont_snr`/`spe_z_rel`) is 100% "flagged" under that definition. Pivoting to
   `spe_error_flag` (confirmed genuinely bimodal, ~4.5% nonzero, in a 5000-row diagnostic over the
   broader "has spectrum" population) did not fix this either: restricted to the
   `spe_galaxy_candidates` join specifically, `spe_error_flag` is uniformly zero across all 6000
   objects sampled on tile 102160061/RGS. `core.read_catalogue_csv` was updated to source the
   "context warning" input from `spe_error_flag` (documented in its docstring) for correctness
   going forward, but the real result on this tile remains a documented null finding — see Sec 6.

## 5. Real datasets accessed

`astroquery.esa.euclid` is unavailable in the pinned `astroquery==0.4.7` (confirmed by
`ImportError`, matching the sibling project); real access went through IRSA's public Euclid Q1
TAP mirror (`astroquery.ipac.irsa.Irsa`).

- **Catalogue query**: tile 102160061, grism='RGS', joining `euclid_q1_spectro_zcatalog_spe_quality`
  + `euclid_q1_spectro_zcatalog_spe_galaxy_candidates` (rank=0) + `euclid_q1_mer_catalogue`, ORDER BY
  object_id, TOP 6000 (final run; an earlier TOP 800 run surfaced the column-name bug in Sec 4 and
  was superseded). Best-line S/N queried separately (batched) from
  `euclid_q1_spe_lines_line_features` for the same 6000 object_ids.
- **Spectra**: 23 real combined 1-D NISP spectra (VOTables with WAVELENGTH/SIGNAL/UNCERTAINTY/
  MASK/QUALITY/NDITH columns) downloaded via IRSA's `spectrumdm/convert` API, ~100-110 KB each
  (~2.4 MB total), for a deterministic sorted-by-object_id subset. 40 were attempted; 17 were
  skipped on transient `502 Bad Gateway` responses from IRSA (recorded as warnings, not silently
  dropped) across two fetch runs.
- **Licence/terms**: Euclid Q1 public data release (ESA Euclid Consortium), mirrored by IRSA/IPAC;
  standard IRSA archive usage terms (https://irsa.ipac.caltech.edu/holdings/copyright.html).
  Confirmed PUBLIC, no proprietary/embargo restriction.
- Full manifest with SHA-256 checksums, source URLs and retrieval timestamps: `data/manifest.csv`
  (45 rows: 1 catalogue-query product + 1 earlier catalogue-query product + 23 spectrum products,
  reflecting both fetch runs). Raw files are **not committed** (`.gitignore` excludes `data/raw/*`).

## 6. Validation and uncertainty outcomes

- **Synthetic injection-recovery gate**: PASSED. `tests/test_lines.py` confirms a strong injected
  Gaussian line (flux ratio 0.7-1.3 of injected) is recovered and detected above the S/N threshold;
  `tests/test_core_pipeline.py::test_run_pipeline_recovers_known_synthetic_effect` confirms the full
  pipeline recovers a deliberately injected clean/flagged group separation in continuum S/N, line
  recovery and redshift reliability, with the permuted-label negative control correctly rejecting
  the null (p<0.05) for the real injected effect.
- **Null control**: `tests/test_uncertainty.py::test_permuted_label_null_distribution_no_real_effect`
  confirms no spurious effect is found when group labels carry no real signal.
- **Failure-mode tests**: missing/empty catalogue CSV, missing VOTable file, missing VOTable
  columns, non-finite values, too-few-pixels, degenerate contamination split, ill-conditioned fit
  covariance — all raise the documented exceptions or emit a warning rather than crashing.
- **Real-data result**: on tile 102160061 (RGS grism), n=6000 objects with a rank-0 galaxy
  redshift candidate, every released SPE quality flag examined is constant, so the clean/flagged
  split is degenerate (n_clean=6000, n_flagged=0 under `spe_error_flag`). `results/summary.json`
  reports `n_total=6000`; the group-comparison, threshold-sensitivity and negative-control metrics
  are `NaN`/omitted with the reason recorded in `results/warnings.json`
  (`"contamination split is degenerate (n_clean=6000, n_flagged=0)"`). This is a genuine finding
  reported transparently, not a code failure.
- **Benchmarks**: `results/benchmarks.json` records wall time, tracemalloc peak memory, dataset
  size, Python version, platform and package version for both the demo and real runs, separating
  data-download time (in `fetch_data.py`'s own log) from analysis time.

## 7. Figures regenerated

6/6 demo figures (SVG+PNG+JSON) verified non-degenerate (real variation across x, not a single
point/value). 6/6 real-data figures regenerated against the real n=6000 catalogue and the 23 real
spectra; the field-map, continuum-S/N-distribution and line-recovery-curve figures show real
single-group (not two-group) distributions, consistent with the degenerate real split reported
above; the low/high-contamination example-spectrum figure uses one real spectrum (no real
"flagged" example spectrum existed in the downloaded 23 to contrast against, so the figure labels
this honestly rather than fabricating a flagged example).

## 8. Claims safe for public release

- "A synthetic-continuum line-injection/recovery test with a known truth passes for this
  pipeline's line-measurement code, and a synthetic catalogue with a known injected contamination
  effect is correctly recovered by the full pipeline (bootstrap, threshold-sensitivity, permuted-
  label negative control)."
- "On Euclid Q1 tile 102160061 (RGS grism), the released `spe_context_warning_flag`,
  `spe_context_error_flag` and `spe_error_flag` columns are constant across all 6000 objects with a
  rank-0 galaxy redshift candidate — this specific real sample does not support a clean/flagged
  group comparison using these released flags."

## Claims that must NOT be made
- Any claim that this project measured a real clean-vs-flagged contamination effect in Euclid Q1
  NISP data — the real-data result is a documented null/degenerate finding, not a positive result.
- Any claim generalising beyond tile 102160061 / RGS grism / the 6000 queried galaxy-candidate
  objects.
- Any claim about the exact bit semantics of the released SPE quality flags beyond what is
  reported in docs/ASSUMPTIONS_AND_LIMITATIONS.md and this report.

## 9. Unresolved limitations / manual review points

- The real-data null finding (Sec 6) means the original scientific question is answered
  ("no measurable real-data relationship was detectable with these released flags on this bounded
  sample") but not in the originally anticipated positive-correlation form. A follow-on release
  could query a broader tile/grism sample, or a non-galaxy-candidate-gated source of continuum
  S/N, to test whether real flag variance exists elsewhere in Q1.
- Only 23/40 spectrum downloads succeeded (transient IRSA 502s); sufficient for the example-
  spectrum figure but a smaller real pool than intended.
- The generic "Euclid Q1 spectroscopy science papers" literature seed could not be verified as a
  specific citable work; marked `TODO_VERIFY` in `references.bib` with an explanation rather than
  invented or silently dropped.
- `reports/report.tex` could not be compiled to PDF in this environment (no local LaTeX toolchain);
  only structural completeness was verified.

## 10. Summary

The pipeline, tests, figures, web dashboard and report are complete and internally consistent. The
real-data run against IRSA's public Euclid Q1 mirror succeeded operationally (queries, downloads,
checksums, provenance all recorded correctly) but surfaced a genuine, reportable null finding about
the released contamination-flag columns on the sampled tile, which is documented transparently
throughout the report, dashboard and this completion report rather than concealed or worked around
with fabricated data.

# Conditional Euclid redshift-quality audit — 2026-11-21

Issue: #1. Three diagnostics: normalized absolute redshift error, flagged-minus-clean catastrophic-outlier fraction stratified by externally supplied S/N bins with seeded within-cell bootstrap, and conditional empirical interval coverage.

Run `python -m unittest discover -s research/2026-11-21 -p 'test_*.py'`. Requirements: NumPy. To produce a labelled synthetic plot run `python research/2026-11-21/plot_conditional_audit.py` (Matplotlib required).

No external reference redshifts or truth-labelled flags were accessed. The existing Q1 sample with no flagged objects cannot estimate a flagged-versus-clean contrast. This diagnostic intentionally raises on strata without both groups. The bootstrap assumes independent sources within cells and does not handle survey clustering or selection bias. All numbers in the example plot are synthetic.

Scientific author: Biswajit Jana. AI-assisted implementation disclosed; no fabricated co-authors.
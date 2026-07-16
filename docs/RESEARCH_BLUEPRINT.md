# Research Blueprint

## Technical title

Euclid NISP Spectral-Contamination Reliability Audit

## Category

Slitless spectroscopy instrumentation / data science

## Bounded scientific question

How do released contamination indicators relate to continuum S/N, emission-line recovery and redshift reliability in selected NISP spectra?

## Gap statement

An external audit; not a decontamination, extraction or redshift pipeline.

## First-release scope

The first release must be completable as a focused 4–6 hour implementation pass after data access is working. It must deliver one reproducible analysis pipeline, one deterministic example/smoke dataset, tests, 4–6 figures, a concise TeX report and a deployable research webpage.

## Validation and uncertainty

- line injection into real continua
- bootstrap contamination groups
- redshift outlier fraction
- threshold sensitivity
- permuted-label negative control

## Required figures

1. low/high contamination spectra
2. S/N vs contamination
3. redshift error
4. outlier fraction
5. line recovery
6. field map

## Reusable scientific modules

- `io.py`
- `masks.py`
- `continuum.py`
- `lines.py`
- `redshift_metrics.py`
- `uncertainty.py`
- `provenance.py`

## Explicit exclusions

- No novelty claim beyond the bounded dataset/question/method combination.
- No causal claim from descriptive catalogue correlations.
- No hidden manual data editing.
- No unsupported precision beyond the input uncertainties.
- No production-pipeline replacement claim.

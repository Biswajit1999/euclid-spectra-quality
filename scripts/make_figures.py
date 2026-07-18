"""Generate the 6 required figures (docs/FIGURE_AND_UI_SPEC.md) as SVG +
300 dpi PNG, each with a sidecar JSON recording git commit, config hash,
sample size and units.

--demo builds figures from the synthetic, clearly-labelled data model in
`euclid_nisp_contamination_audit.synthetic`. The real-data path reads
data/source_catalog.csv (written by scripts/fetch_data.py) and the real
downloaded 1-D spectra under data/raw/, and must only be run after
scripts/run_analysis.py (real mode) has produced validated results.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401 - importing registers the bundled matplotlib styles

from euclid_nisp_contamination_audit import __version__
from euclid_nisp_contamination_audit.config import load_config
from euclid_nisp_contamination_audit.core import read_catalogue_csv
from euclid_nisp_contamination_audit.io import read_spectrum_votable
from euclid_nisp_contamination_audit.masks import CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED
from euclid_nisp_contamination_audit.plotting import CONTAMINATION_GROUP_COLORS, save_figure_with_sidecar
from euclid_nisp_contamination_audit.provenance import get_git_commit, sha256_config
from euclid_nisp_contamination_audit.synthetic import build_synthetic_catalogue, build_synthetic_spectrum

plt.style.use(["science", "no-latex"])

REPO_ROOT = Path(__file__).resolve().parents[1]


def _save(fig, out_dir: Path, name: str, *, data_kind: str, sample_size: int, units: str, config_path: Path) -> None:
    save_figure_with_sidecar(
        fig, out_dir / name, figure_name=name, data_kind=data_kind, sample_size=sample_size, units=units,
        git_commit=get_git_commit(REPO_ROOT),
        config_sha256=sha256_config(config_path) if config_path.is_file() else "unavailable",
        package_version=__version__,
    )
    plt.close(fig)


def _plot_spectra(ax, wavelength, signal, mask, quality, label, color) -> None:
    good = (mask == 0) & (quality == 0)
    ax.plot(wavelength, signal, color=color, lw=0.8, alpha=0.6, label=f"{label} (all pixels)")
    ax.plot(wavelength[good], signal[good], color=color, lw=1.2, label=f"{label} (good pixels)")


def _fig_snr_vs_contamination(out_dir, config_path, data_kind, cont_snr, group, sample_size):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        values = cont_snr[(group == grp) & np.isfinite(cont_snr)]
        if values.size == 0:
            continue
        ax.hist(values, bins=20, alpha=0.6, color=CONTAMINATION_GROUP_COLORS[grp], label=f"{grp} (n={values.size})")
    ax.set_xlabel("Continuum S/N (spe_cont_snr, dimensionless)")
    ax.set_ylabel("N objects")
    ax.set_title(f"Continuum S/N by contamination group (n={sample_size})")
    ax.legend()
    fig.tight_layout()
    _save(fig, out_dir, "snr_vs_contamination", data_kind=data_kind, sample_size=sample_size, units="dimensionless (S/N)", config_path=config_path)


def _fig_redshift_error(out_dir, config_path, data_kind, z_err, group, sample_size):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        values = z_err[(group == grp) & np.isfinite(z_err)]
        if values.size == 0:
            continue
        ax.hist(values, bins=20, alpha=0.6, color=CONTAMINATION_GROUP_COLORS[grp], label=f"{grp} (n={values.size})")
    ax.set_xlabel("Formal redshift uncertainty spe_z_err (dimensionless)")
    ax.set_ylabel("N objects")
    ax.set_title(f"Redshift error by contamination group (n={sample_size})")
    ax.legend()
    fig.tight_layout()
    _save(fig, out_dir, "redshift_error", data_kind=data_kind, sample_size=sample_size, units="dimensionless (delta_z)", config_path=config_path)


def _fig_outlier_fraction(out_dir, config_path, data_kind, z_rel, group, sample_size, seed, n_resamples):
    from euclid_nisp_contamination_audit.uncertainty import bootstrap_statistic

    fig, ax = plt.subplots(figsize=(6, 4.5))
    labels, estimates, lows, highs = [], [], [], []
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        indicator = (z_rel[(group == grp) & np.isfinite(z_rel)] < 0.5).astype(float)
        if indicator.size < 2:
            continue
        boot = bootstrap_statistic(indicator, np.mean, n_resamples, seed)
        labels.append(grp)
        estimates.append(boot.estimate)
        lows.append(boot.estimate - boot.ci_low)
        highs.append(boot.ci_high - boot.estimate)
    colors = [CONTAMINATION_GROUP_COLORS[g] for g in labels]
    ax.bar(labels, estimates, yerr=[lows, highs], color=colors, capsize=6)
    ax.set_ylabel("Fraction with spe_z_rel < 0.5 (unreliable)")
    ax.set_title(f"Redshift-reliability outlier fraction by group, 95% bootstrap CI (n={sample_size})")
    fig.tight_layout()
    _save(fig, out_dir, "outlier_fraction", data_kind=data_kind, sample_size=sample_size, units="fraction", config_path=config_path)


def _fig_line_recovery(out_dir, config_path, data_kind, best_line_snr, group, sample_size):
    from euclid_nisp_contamination_audit.lines import DEFAULT_LINE_RECOVERY_SNR_THRESHOLD

    thresholds = np.linspace(1, 15, 20)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        values = best_line_snr[(group == grp) & np.isfinite(best_line_snr)]
        if values.size == 0:
            continue
        rates = [float(np.mean(values >= t)) for t in thresholds]
        ax.plot(thresholds, rates, color=CONTAMINATION_GROUP_COLORS[grp], label=f"{grp} (n={values.size})", marker="o", ms=3)
    ax.axvline(DEFAULT_LINE_RECOVERY_SNR_THRESHOLD, color="gray", ls="--", lw=1, label="adopted S/N threshold")
    ax.set_xlabel("Best-line S/N threshold (spe_line_snr_gf)")
    ax.set_ylabel("Line recovery rate (fraction of objects)")
    ax.set_title(f"Line recovery vs S/N threshold, by contamination group (n={sample_size})")
    ax.legend()
    fig.tight_layout()
    _save(fig, out_dir, "line_recovery", data_kind=data_kind, sample_size=sample_size, units="fraction", config_path=config_path)


def _fig_field_map(out_dir, config_path, data_kind, ra, dec, group, sample_size):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for grp in (CONTAMINATION_GROUP_CLEAN, CONTAMINATION_GROUP_FLAGGED):
        mask = group == grp
        ax.scatter(ra[mask], dec[mask], s=14, alpha=0.7, color=CONTAMINATION_GROUP_COLORS[grp], label=f"{grp} (n={int(np.sum(mask))})")
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    ax.set_title(f"Field map by contamination group (n={sample_size})")
    ax.legend()
    ax.invert_xaxis()
    fig.tight_layout()
    _save(fig, out_dir, "field_map", data_kind=data_kind, sample_size=sample_size, units="degrees", config_path=config_path)


def make_demo_figures(out_dir: Path, config_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    data_kind = "synthetic_demo"
    config = load_config(config_path)

    catalogue = build_synthetic_catalogue(n_sources=300, seed=config.execution.seed)
    group = catalogue.contamination_group
    n = catalogue.object_id.size

    clean_spec = build_synthetic_spectrum(object_id=1, seed=101, contaminated=False, true_z=1.0)
    flagged_spec = build_synthetic_spectrum(object_id=2, seed=102, contaminated=True, true_z=1.0)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    _plot_spectra(ax, clean_spec.wavelength_angstrom, clean_spec.signal, clean_spec.mask, clean_spec.quality, "clean", CONTAMINATION_GROUP_COLORS["clean"])
    _plot_spectra(ax, flagged_spec.wavelength_angstrom, flagged_spec.signal, flagged_spec.mask, flagged_spec.quality, "flagged", CONTAMINATION_GROUP_COLORS["flagged"])
    ax.set_xlabel("Observed wavelength (Angstrom)")
    ax.set_ylabel("Signal (1e-16 erg/s/cm^2/Angstrom)")
    ax.set_title("Example low/high-contamination spectra - SYNTHETIC DEMO")
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, out_dir, "low_high_contamination_spectra", data_kind=data_kind, sample_size=2, units="flux density vs Angstrom", config_path=config_path)

    _fig_snr_vs_contamination(out_dir, config_path, data_kind, catalogue.spe_cont_snr, group, n)
    _fig_redshift_error(out_dir, config_path, data_kind, catalogue.spe_z_err, group, n)
    _fig_outlier_fraction(out_dir, config_path, data_kind, catalogue.spe_z_rel, group, n, config.execution.seed, config.validation.bootstrap_resamples)
    _fig_line_recovery(out_dir, config_path, data_kind, catalogue.best_line_snr, group, n)
    _fig_field_map(out_dir, config_path, data_kind, catalogue.ra, catalogue.dec, group, n)

    print(f"Wrote 6 demo figures (SVG+PNG+JSON) to {out_dir}")


def make_real_figures(out_dir: Path, config_path: Path, catalogue_csv: Path, raw_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path)
    data_kind = config.input.data_mode

    rows = read_catalogue_csv(catalogue_csv)
    n = len(rows)
    from euclid_nisp_contamination_audit.masks import contamination_group_labels

    context_warning = np.array([r.context_warning_flag for r in rows])
    context_error = np.array([r.context_error_flag for r in rows])
    group = contamination_group_labels(context_warning, context_error)
    cont_snr = np.array([r.spe_cont_snr for r in rows], dtype=float)
    z_err = np.array([r.spe_z_err for r in rows], dtype=float)
    z_rel = np.array([r.spe_z_rel for r in rows], dtype=float)
    best_line_snr = np.array([r.best_line_snr for r in rows], dtype=float)
    ra = np.array([r.ra for r in rows], dtype=float)
    dec = np.array([r.dec for r in rows], dtype=float)
    object_ids = [r.object_id for r in rows]
    group_by_id = dict(zip(object_ids, group))

    example_clean, example_flagged = None, None
    for vot_path in sorted(raw_dir.glob("euclid_q1_nisp_spectrum_*.vot")):
        object_id = int(vot_path.stem.rsplit("_", 1)[-1])
        if object_id not in group_by_id:
            continue
        spectrum = read_spectrum_votable(vot_path, object_id)
        if group_by_id[object_id] == CONTAMINATION_GROUP_CLEAN and example_clean is None:
            example_clean = spectrum
        elif group_by_id[object_id] == CONTAMINATION_GROUP_FLAGGED and example_flagged is None:
            example_flagged = spectrum
        if example_clean is not None and example_flagged is not None:
            break

    fig, ax = plt.subplots(figsize=(8, 4.5))
    if example_clean is not None:
        _plot_spectra(ax, example_clean.wavelength_angstrom, example_clean.signal, example_clean.mask, example_clean.quality, f"clean (id={example_clean.object_id})", CONTAMINATION_GROUP_COLORS["clean"])
    if example_flagged is not None:
        _plot_spectra(ax, example_flagged.wavelength_angstrom, example_flagged.signal, example_flagged.mask, example_flagged.quality, f"flagged (id={example_flagged.object_id})", CONTAMINATION_GROUP_COLORS["flagged"])
    if example_clean is None and example_flagged is None:
        ax.text(0.5, 0.5, "no downloaded example spectra available in data/raw/", transform=ax.transAxes, ha="center", va="center")
    ax.set_xlabel("Observed wavelength (Angstrom)")
    ax.set_ylabel("Signal (1e-16 erg/s/cm^2/Angstrom)")
    ax.set_title("Example real low/high-contamination NISP spectra")
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, out_dir, "low_high_contamination_spectra", data_kind=data_kind, sample_size=int(example_clean is not None) + int(example_flagged is not None), units="flux density vs Angstrom", config_path=config_path)

    _fig_snr_vs_contamination(out_dir, config_path, data_kind, cont_snr, group, n)
    _fig_redshift_error(out_dir, config_path, data_kind, z_err, group, n)
    _fig_outlier_fraction(out_dir, config_path, data_kind, z_rel, group, n, config.execution.seed, config.validation.bootstrap_resamples)
    _fig_line_recovery(out_dir, config_path, data_kind, best_line_snr, group, n)
    _fig_field_map(out_dir, config_path, data_kind, ra, dec, group, n)

    print(f"Wrote 6 real-data figures (SVG+PNG+JSON) to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    parser.add_argument("--config", type=Path, default=Path("config/analysis.yml"))
    parser.add_argument("--catalogue-csv", type=Path, default=Path("data/source_catalog.csv"))
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    args = parser.parse_args()

    if args.demo:
        make_demo_figures(args.out_dir, args.config)
        return

    make_real_figures(args.out_dir, args.config, args.catalogue_csv, args.raw_dir)


if __name__ == "__main__":
    main()

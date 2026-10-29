"""Reproducible synthetic-only demo; does not use real Euclid tile as truth."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from redshift_robustness import (
    matched_snr_outlier_contrast, redshift_coverage, reliability_calibration,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/redshift_robustness"))
    args = parser.parse_args()
    rng = np.random.default_rng(20261029)
    n = 600
    group = np.array(["clean"] * (n // 2) + ["flagged"] * (n // 2))
    truth = rng.uniform(.7, 1.7, n)
    snr = rng.uniform(4, 16, n)
    sigma = .02 * (1 + truth)
    z = truth + rng.normal(0, sigma)
    injected = (group == "flagged") & (rng.random(n) < .3)
    z[injected] += .5 * (1 + truth[injected])
    reliability = np.where(group == "clean", .9, .7)
    result = {"data_type": "synthetic", "seed": 20261029, "n_objects": n,
              "injected_flagged_catastrophic_fraction": .3,
              "caveat": "No empirical Euclid contamination effect or external truth calibration claimed.",
              "coverage": redshift_coverage(z, truth, sigma, group),
              "reliability": reliability_calibration(z, truth, reliability),
              "matched_snr": matched_snr_outlier_contrast(z, truth, group, snr)}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    rows = result["matched_snr"]["strata"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar([str(row["bin"]) for row in rows], [row["difference"] for row in rows])
    ax.axhline(0, color="black", linewidth=.8)
    ax.set(xlabel="Pooled SNR quantile stratum", ylabel="Flagged − clean outlier fraction",
           title="Synthetic contamination injection: SNR-matched contrast")
    fig.tight_layout()
    fig.savefig(args.output_dir / "matched_snr.png", dpi=160)
    plt.close(fig)
    print(json.dumps({"output": str(args.output_dir), "matched_difference":
                      result["matched_snr"]["outlier_rate_difference_b_minus_a"]}))


if __name__ == "__main__":
    main()

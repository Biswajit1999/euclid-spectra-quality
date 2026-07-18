from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401 - importing registers the bundled matplotlib styles

plt.style.use(["science", "no-latex"])

CONTAMINATION_GROUP_COLORS = {"clean": "#2a6f97", "flagged": "#c1440e"}


def plot_demo(values: np.ndarray, output: str | Path) -> Path:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(np.arange(values.size), values)
    ax.set_xlabel("Synthetic index")
    ax.set_ylabel("Synthetic value")
    ax.set_title("Smoke-test output - not a scientific result")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def save_figure_with_sidecar(
    fig: plt.Figure,
    output_stem: str | Path,
    *,
    figure_name: str,
    data_kind: str,
    sample_size: int,
    units: str,
    git_commit: str,
    config_sha256: str,
    package_version: str,
) -> dict[str, Path]:
    """Save a figure as SVG + 300dpi PNG plus a sidecar JSON of provenance
    metadata, per docs/FIGURE_AND_UI_SPEC.md.
    """
    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    svg_path = stem.with_suffix(".svg")
    png_path = stem.with_suffix(".png")
    json_path = stem.with_suffix(".json")

    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)

    sidecar: dict[str, Any] = {
        "figure_name": figure_name,
        "data_kind": data_kind,
        "sample_size": sample_size,
        "units": units,
        "git_commit": git_commit,
        "config_sha256": config_sha256,
        "package_version": package_version,
    }
    json_path.write_text(json.dumps(sidecar, indent=2), encoding="utf-8")
    return {"svg": svg_path, "png": png_path, "json": json_path}


__all__ = ["CONTAMINATION_GROUP_COLORS", "plot_demo", "save_figure_with_sidecar"]

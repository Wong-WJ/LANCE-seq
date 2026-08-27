"""Small reusable plotting and output helpers for LANCE-seq Figure 4."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import seaborn as sns


EXPECTED_FILES = [
    "Fig4_5A_Pseudobulk_PCA_PC1_PC2.png",
    "Fig4_5B_6module_heatmap.png",
    "Fig4_6A3_MN2_MN4_expression_agreement_labeled.png",
    "Fig4_6B3_APAP_injury_module_comparison.png",
]


def configure_style() -> None:
    """Apply the common Figure 4 plotting style."""
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    sns.set_theme(style="white", font="Arial")


def prepare_output_dir(output_dir: str | Path) -> Path:
    """Create the output folder and clear stale PNG files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in output_dir.glob("*.png"):
        path.unlink()
    return output_dir


def save_png(fig: mpl.figure.Figure, output_dir: str | Path, filename: str) -> Path:
    """Save one publication PNG with the shared export settings."""
    path = Path(output_dir) / filename
    fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def add_descriptive_ellipse(
    ax: mpl.axes.Axes,
    points: np.ndarray,
    color: str,
    global_x_span: float,
    global_y_span: float,
) -> None:
    """Add the descriptive group ellipse used in the PCA panel."""
    center = points.mean(axis=0)
    width = max(np.ptp(points[:, 0]) * 1.30, global_x_span * 0.10)
    height = max(np.ptp(points[:, 1]) * 1.30, global_y_span * 0.10)
    radius = np.sqrt(
        ((points[:, 0] - center[0]) / (width / 2)) ** 2
        + ((points[:, 1] - center[1]) / (height / 2)) ** 2
    )
    scale = max(1.0, float(np.max(radius)) / 0.88)
    ax.add_patch(
        Ellipse(
            center,
            width=width * scale,
            height=height * scale,
            facecolor=mpl.colors.to_rgba(color, 0.06),
            edgecolor=color,
            linewidth=1.2,
            zorder=1,
        )
    )


def density_values(
    x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, PowerNorm, LinearSegmentedColormap]:
    """Return binned point density and its display mapping."""
    x_pad = max(np.ptp(x) * 0.002, 1e-9)
    y_pad = max(np.ptp(y) * 0.002, 1e-9)
    hist, x_edges, y_edges = np.histogram2d(
        x,
        y,
        bins=180,
        range=[
            [x.min() - x_pad, x.max() + x_pad],
            [y.min() - y_pad, y.max() + y_pad],
        ],
    )
    x_bin = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, 179)
    y_bin = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, 179)
    density = hist[x_bin, y_bin]
    positive = density[density > 0]
    dmin, dmax = np.quantile(positive, [0.01, 0.99])
    norm = PowerNorm(
        gamma=0.60, vmin=dmin, vmax=max(dmax, dmin + 1e-12), clip=True
    )
    cmap = LinearSegmentedColormap.from_list(
        "soft_density_final",
        ["#F3F3F3", "#DCE7F2", "#AFC7DE", "#E7CBC7", "#E6B0AA", "#B9625C"],
    )
    return density, norm, cmap


def output_inventory(output_dir: str | Path) -> pd.DataFrame:
    """Summarize the four Figure 4 PNG files."""
    output_dir = Path(output_dir)
    observed = sorted(path.name for path in output_dir.glob("*.png"))
    missing = sorted(set(EXPECTED_FILES) - set(observed))
    extra = sorted(set(observed) - set(EXPECTED_FILES))
    if missing or extra:
        raise ValueError(f"Unexpected PNG inventory; missing={missing}, extra={extra}")
    return pd.DataFrame(
        [
            {"filename": name, "bytes": (output_dir / name).stat().st_size}
            for name in EXPECTED_FILES
        ]
    )

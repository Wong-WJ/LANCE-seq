"""Small reusable plotting and mathematical helpers for Figure 5."""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Mapping, Sequence

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.spatial import cKDTree


EXPECTED_OUTPUTS = [
    "Fig5B_MN2_vs_MN1_Volcano.png",
    "Fig5C1_GO_Up.png",
    "Fig5C2_GO_Down.png",
    "Fig5D_APAP_Injury_Gene_Spatial_Maps_fixed.png",
    "Fig5E_Liver_Zonation_Marker_Spatial_Maps.png",
    "Fig5F_Relative_Zonation_3State_Spatial_Maps_fixed.png",
    "Fig5G_ZoneStratified_Transcriptional_Effect_Sizes.png",
    "Fig5H_APAP_Module_Response_Across_Relative_Zonation_Strata.png",
    "FigS5A1_KEGG_Up.png",
    "FigS5A2_KEGG_Down.png",
    "FigS5B_APAP_Marker_Expression_Distribution_fixed.png",
    "FigS5C_Zonation_Marker_Expression_Distribution.png",
    "FigS5D_Relative_Zonation_Continuous_Spatial_Maps_fixed.png",
    "FigS5E_Healthy_HeldOut_Zonation_Validation.png",
    "FigS5G_HealthyReference_Absolute_Zonation_Program_Remodeling.png",
    "FigS5H_Absolute_CV_Associated_Program_Spatial_Maps.png",
    "FigS5I_ZoneStratified_Supportive_Concordance.png",
    "FigS5J_Spatial_Autocorrelation_MoransI.png",
]


def configure_style(random_seed: int = 0) -> None:
    np.random.seed(random_seed)
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#2F2F2F",
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 9,
        }
    )
    sns.set_style("white")


def prepare_output(output_dir: str | Path) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in output_dir.glob("*.png"):
        path.unlink()
    return output_dir


def save_png(output_dir: str | Path, fig: plt.Figure, filename: str) -> Path:
    if filename not in EXPECTED_OUTPUTS:
        raise ValueError(f"Unexpected Figure 5 filename: {filename}")
    path = Path(output_dir) / filename
    fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
    return path


def gene_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    values = adata[:, gene].X
    return values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel()


def normalize_selected(X: sparse.csr_matrix, indices: np.ndarray) -> np.ndarray:
    """Normalize selected genes to 10,000 counts per spot and apply log1p."""
    library = np.asarray(X.sum(axis=1)).ravel()
    selected = X[:, indices].toarray()
    return np.log1p(
        np.divide(
            selected,
            library[:, None],
            out=np.zeros_like(selected, dtype=float),
            where=library[:, None] > 0,
        )
        * 1e4
    )


def enrichment_metrics(table: pd.DataFrame) -> pd.DataFrame:
    result = table.copy()
    overlap = result["Overlap"].astype(str).str.extract(
        r"(?P<hit>\d+)\s*/\s*(?P<total>\d+)"
    ).astype(float)
    result["Overlap_Ratio"] = overlap["hit"] / overlap["total"].replace(0, np.nan)
    result["minus_log10_padj"] = -np.log10(
        np.maximum(result["Adjusted P-value"].astype(float), 1e-300)
    )
    result["Term_without_ID"] = result["Term"].astype(str).str.replace(
        r"\s*\(GO:\d+\)$", "", regex=True
    )
    return result


def _integer_ratio_refs(values: pd.Series) -> np.ndarray:
    pct = np.asarray(values, float) * 100
    low, high = float(pct.min()), float(pct.max())
    if np.isclose(low, high):
        center = max(1, int(round(low)))
        return np.array([max(1, center - 5), center, center + 5])
    refs = np.rint(np.linspace(low, high, 3)).astype(int)
    refs[0] = max(1, refs[0])
    refs[1] = max(refs[0] + 1, refs[1])
    refs[2] = max(refs[1] + 1, refs[2])
    return refs


def plot_enrichment(
    selected: pd.DataFrame,
    output_dir: str | Path,
    filename: str,
    title: str,
    color: str,
    ylabel: str,
) -> None:
    selected = selected.sort_values("minus_log10_padj").copy()
    selected["Display_Term"] = selected["Term_without_ID"].map(
        lambda value: textwrap.fill(value, width=36)
    )
    rmin, rmax = selected["Overlap_Ratio"].min(), selected["Overlap_Ratio"].max()

    def size(values: Sequence[float]) -> np.ndarray:
        values = np.asarray(values, float)
        if np.isclose(rmin, rmax):
            return np.full(values.shape, 120.0)
        return 50 + (values - rmin) / (rmax - rmin) * 140

    fig, ax = plt.subplots(figsize=(5.8, max(3.6, 0.48 * len(selected) + 1.25)))
    y = np.arange(len(selected))
    x = selected["minus_log10_padj"].to_numpy(float)
    xleft = max(0.0, float(x.min()) - max(0.35, 0.06 * float(np.ptp(x))))
    for yi, xi in zip(y, x):
        ax.hlines(yi, xleft, xi, color="#D0D0D0", linestyle="--", linewidth=0.7, alpha=0.75)
    ax.scatter(x, y, s=size(selected["Overlap_Ratio"]), c=color, alpha=0.9, edgecolors="white", linewidth=0.7)
    refs = _integer_ratio_refs(selected["Overlap_Ratio"])
    handles = [
        ax.scatter([], [], s=float(size([value / 100])[0]), c=color, edgecolors="white", linewidth=0.7)
        for value in refs
    ]
    ax.legend(handles, [f"{value}%" for value in refs], title="Overlap (%)", frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_yticks(y, selected["Display_Term"])
    ax.set_xlim(left=xleft)
    ax.set_xlabel("−log10 adjusted P-value")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10, pad=8)
    ax.tick_params(axis="y", labelsize=8, length=0)
    ax.grid(False)
    fig.tight_layout(pad=0.7)
    save_png(output_dir, fig, filename)
    plt.close(fig)


def plot_spatial_gene_grid(
    frames: Mapping[str, pd.DataFrame],
    coords: Mapping[str, np.ndarray],
    genes: Sequence[str],
    output_dir: str | Path,
    filename: str,
) -> None:
    fig, axes = plt.subplots(len(genes), 2, figsize=(8.8, 3.0 * len(genes)), constrained_layout=True)
    for row, gene in enumerate(genes):
        combined = np.concatenate([frames["MN1"][gene], frames["MN2"][gene]])
        vmin, vmax = np.quantile(combined, [0.01, 0.99])
        if np.isclose(vmin, vmax):
            vmax = max(float(combined.max()), float(vmin) + 1)
        for col, sample in enumerate(["MN1", "MN2"]):
            ax = axes[row, col]
            ax.scatter(coords[sample][:, 0], coords[sample][:, 1], c=frames[sample][gene], s=7.0, cmap="Spectral_r", vmin=vmin, vmax=vmax, linewidths=0)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row == 0:
                ax.set_title(sample, fontsize=10)
            if col == 0:
                ax.text(-0.08, 0.5, gene, transform=ax.transAxes, ha="right", va="center", fontsize=9)
        scalar = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap="Spectral_r")
        fig.colorbar(scalar, ax=axes[row, :], fraction=0.025, pad=0.02).ax.tick_params(labelsize=7)
    save_png(output_dir, fig, filename)
    plt.close(fig)


def plot_box_grid(
    frames: Mapping[str, pd.DataFrame],
    genes: Sequence[str],
    output_dir: str | Path,
    filename: str,
) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 4.2), sharey=False)
    for ax, gene in zip(axes, genes):
        frame = pd.concat(
            [pd.DataFrame({"sample": sample, "expression": frames[sample][gene].to_numpy()}) for sample in ["MN1", "MN2"]],
            ignore_index=True,
        )
        sns.boxplot(data=frame, x="sample", y="expression", hue="sample", order=["MN1", "MN2"], palette={"MN1": "#AFC7E8", "MN2": "#E6B0AA"}, width=0.5, linewidth=2, showfliers=False, legend=False, ax=ax)
        ax.set_title(gene, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("Log-normalized expression")
    fig.tight_layout()
    save_png(output_dir, fig, filename)
    plt.close(fig)


def rank_tertiles(score: pd.Series) -> pd.Series:
    order = np.argsort(score.to_numpy(float), kind="mergesort")
    labels = np.empty(len(score), dtype=object)
    cut1, cut2 = len(score) // 3, 2 * len(score) // 3
    labels[order[:cut1]] = "Relative PV-like"
    labels[order[cut1:cut2]] = "Intermediate"
    labels[order[cut2:]] = "Relative CV-like"
    return pd.Series(labels, index=score.index, dtype="object")


def clean_spatial_axis(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, fontsize=8.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def common_spatial_limits(
    coords: Mapping[str, np.ndarray], samples: Sequence[str]
) -> tuple[dict[str, tuple], dict[str, float]]:
    nn = {
        sample: float(np.median(cKDTree(coords[sample]).query(coords[sample], k=2)[0][:, 1]))
        for sample in samples
    }
    common_span = max(float(np.ptp(coords[sample][:, axis])) for sample in samples for axis in (0, 1)) * 1.1
    limits = {}
    for sample in samples:
        center = coords[sample].mean(axis=0)
        half = common_span / 2
        limits[sample] = ((center[0] - half, center[0] + half), (center[1] - half, center[1] + half))
    return limits, nn


def symmetric_knn(coords: np.ndarray, k: int = 6) -> sparse.csr_matrix:
    _, indices = cKDTree(coords).query(coords, k=min(k + 1, len(coords)))
    rows = np.repeat(np.arange(len(coords)), indices.shape[1] - 1)
    cols = indices[:, 1:].reshape(-1)
    weights = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(coords), len(coords)))
    weights = weights.maximum(weights.T)
    weights.setdiag(0)
    weights.eliminate_zeros()
    return weights


def moran_many(
    matrix: np.ndarray,
    weights: sparse.csr_matrix,
    permutations: int = 999,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(matrix, float)
    z = matrix - matrix.mean(axis=0, keepdims=True)
    denominator = np.sum(z * z, axis=0)
    factor = matrix.shape[0] / float(weights.sum())
    observed = factor * np.sum(z * (weights @ z), axis=0) / denominator
    rng = np.random.default_rng(seed)
    exceed = np.zeros(matrix.shape[1], dtype=int)
    for _ in range(permutations):
        permuted = z[rng.permutation(matrix.shape[0]), :]
        simulated = factor * np.sum(permuted * (weights @ permuted), axis=0) / denominator
        exceed += np.abs(simulated) >= np.abs(observed)
    return observed, (exceed + 1) / (permutations + 1)


def output_inventory(output_dir: str | Path) -> pd.DataFrame:
    output_dir = Path(output_dir)
    observed = sorted(path.name for path in output_dir.glob("*.png"))
    if observed != sorted(EXPECTED_OUTPUTS):
        raise RuntimeError(f"Output mismatch: {observed}")
    return pd.DataFrame(
        [{"filename": name, "bytes": (output_dir / name).stat().st_size} for name in EXPECTED_OUTPUTS]
    )

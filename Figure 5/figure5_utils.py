"""Plot/export helpers for Figure 5. / 图5绘图与导出辅助函数。"""
from pathlib import Path
from typing import Mapping, Sequence
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

EXPECTED_STEMS = ['Fig5B_Longitudinal_MA', 'Fig5C_GO_Up', 'Fig5D_Injury_Spatial', 'Fig5E_Injury_Distributions', 'Fig5F_CV_Reference_Spatial', 'Fig5G_Matched_Marker_Profiles', 'Fig5H_Expression_Moran_Changes', 'Fig5I_Focal_Gene_Spatial', 'FigS5A_GO_Down', 'FigS5B_KEGG_Up', 'FigS5C_KEGG_Down', 'FigS6_Section_Injury_Markers', 'FigS7A_Zonation_Marker_Spatial', 'FigS7B_Zonation_Marker_Distributions', 'FigS8A_Relative_Zonation_Spatial', 'FigS8B_Heldout_Validation', 'FigS9A_CV_Reference_Distributions', 'FigS9B_PV_Reference_Distributions', 'FigS9C_PV_Reference_Spatial']

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


def prepare_output(output_dir):
    """Create output folders without deleting existing files. / 建立目录，保留已有文件。"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    return output_dir


def save_panel(fig, output_dir, stem, *, tight=False, dpi=600):
    """Export PNG/PDF and display a notebook preview. / 导出 PNG/PDF 并显示预览。"""
    output_dir = Path(output_dir)
    options = {"bbox_inches": "tight"} if tight else {}
    fig.canvas.draw()
    fig.savefig(output_dir / (stem + ".png"), dpi=dpi, facecolor="white", **options)
    fig.savefig(output_dir / (stem + ".pdf"), dpi=350, facecolor="white", **options)
    plt.show()
    plt.close(fig)
    return output_dir / (stem + ".png")


def output_inventory(output_dir):
    """Check every expected PNG and PDF. / 核对全部预期图件。"""
    rows = []
    for stem in EXPECTED_STEMS:
        for extension in ("png", "pdf"):
            path = Path(output_dir) / (stem + "." + extension)
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
            rows.append({"file": path.name, "bytes": path.stat().st_size})
    return pd.DataFrame(rows)


def plot_gene_maps(frames, coords, genes, output_dir, stem):
    """Two samples by four genes, with shared per-gene limits. / 两样本四基因共用色标。"""
    fig, axes = plt.subplots(2, len(genes), figsize=(7.5, 3.7))
    fig.subplots_adjust(left=.07, right=.96, bottom=.07, top=.91, hspace=.12, wspace=.17)
    samples = ["MN1", "MN2"]
    span = max(np.ptp(coords[s], axis=0).max() for s in samples) + 4
    for col, gene in enumerate(genes):
        lo, hi = np.quantile(np.r_[frames["MN1"][gene], frames["MN2"][gene]], [.01,.99])
        if np.isclose(lo,hi): hi = max(float(max(frames[s][gene].max() for s in samples)), lo+1)
        for row, sample in enumerate(samples):
            ax=axes[row,col]; xy=coords[sample]; mid=(xy.min(0)+xy.max(0))/2
            artist=ax.scatter(*xy.T,c=frames[sample][gene],s=3.2,linewidths=0,cmap="Spectral_r",vmin=lo,vmax=hi)
            ax.set(xlim=(mid[0]-span/2,mid[0]+span/2),ylim=(mid[1]-span/2,mid[1]+span/2),aspect="equal")
            ax.set_axis_off()
            if row==0: ax.set_title(gene,fontstyle="italic",fontsize=9,pad=4)
            if col==0: ax.text(-.08,.5,sample,transform=ax.transAxes,rotation=90,ha="center",va="center",fontsize=8)
        box=axes[1,col].get_position()
        cax=fig.add_axes([box.x1+.004,box.y0+.05,.006,.10])
        bar=fig.colorbar(artist,cax=cax,ticks=[]);bar.outline.set_linewidth(.4)
    save_panel(fig,output_dir,stem)


def plot_gene_boxes(frames, genes, output_dir, stem):
    """All spots; median/IQR and 1.5-IQR whiskers. / 全部spot，中位数、四分位数和1.5倍IQR须。"""
    fig, axes=plt.subplots(1,len(genes),figsize=(7.5,2.3))
    fig.subplots_adjust(left=.09,right=.98,top=.87,bottom=.21,wspace=.48)
    for col,(ax,gene) in enumerate(zip(axes,genes)):
        box=ax.boxplot([frames[s][gene] for s in ["MN1","MN2"]],patch_artist=True,widths=.5,showfliers=False,whis=1.5,
                       medianprops={"color":"black","linewidth":.8},boxprops={"linewidth":.8},whiskerprops={"linewidth":.8},capprops={"linewidth":.8})
        for patch,color in zip(box["boxes"],["#B5C8E0","#DDB3AE"]):patch.set_facecolor(color)
        ax.set_xticks([1,2],["MN1","MN2"]);ax.spines[["top","right"]].set_visible(False)
        ax.set_title(gene,fontstyle="italic",fontsize=9)
        ax.tick_params(labelsize=8)
        if col==0:ax.set_ylabel("Log-normalized expression",fontsize=8)
    save_panel(fig,output_dir,stem)

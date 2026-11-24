"""Plot/export helpers for Figure 4. / 图4绘图与导出辅助函数。"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

EXPECTED_STEMS = ['Fig4F_Pseudobulk_PCA', 'Fig4G_Five_Module_Heatmap', 'Fig4H_Module_Gene_Changes', 'Fig4I_Gene_Difference_ECDF']

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

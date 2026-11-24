# LANCE-seq Figure 4

Reproduce the final **Figure 4F–I** from the eight liver captures in `../data/LANCE8.h5ad`.

## Files and run

- `LANCE-seq_Figure4.ipynb`: calculation and plotting cells, with Chinese/English explanations.
- `figure4_utils.py`: plotting and export helpers.
- `panel_sources.csv`: final panel-to-source-notebook mapping.
- `outputs/`: four PNG/PDF pairs; `outputs/tables/`: numerical results.

From the repository root run `python run_figures.py --figure 4`, or open the notebook in this folder and run all cells. Dependencies are pinned in `../requirements-figures45.txt`.

## Analysis

1. Sum count-like X across spots in each sample, retaining genes detected in at least three spots. Preserve fractional counts. Compute CPM over these retained genes; keep genes with CPM ≥1 in at least two samples and transform as `log2(CPM+1)`.
2. **F:** centered, unscaled PCA of the eight-sample matrix.
3. **G:** standardize each gene across the eight samples (`ddof=0`) and average within five gene sets: immediate early response, general cellular stress, inflammatory/NF-κB signaling, oxidative stress/NRF2 signaling and cell death/injury. These comprise 50 memberships and 44 distinct genes.
4. **H:** per-gene expression changes for MN3−MN5, MN4−MN5 and MN2−MN1; black lines show medians. A gene belonging to multiple sets appears in each set.
5. **I:** empirical cumulative distributions of absolute expression differences for the same 44 distinct genes in MN2/MN4, MN2/MN3 and MN4/MN3. Each gene is counted once per comparison.

H and I summarize the specified samples. MN2/MN1 is the longitudinal pair; the sham and endpoint contrasts use distinct animals.

## Outputs

| Final panel | Output stem (`.png` / `.pdf`) | Content |
|---|---|---|
| Fig. 4F | `Fig4F_Pseudobulk_PCA` | Sample-level pseudobulk PCA |
| Fig. 4G | `Fig4G_Five_Module_Heatmap` | Five final functional gene sets |
| Fig. 4H | `Fig4H_Module_Gene_Changes` | Gene-level changes in three sample contrasts |
| Fig. 4I | `Fig4I_Gene_Difference_ECDF` | Cumulative absolute differences across 44 genes |

Figure 4A–E are supplied in the manuscript figure artwork and are outside this notebook's scope.

## 中文说明

保留原PCA计算；图4G改为最终图的5个基因集；图4H/I替换为最终版基因表达变化及44个基因的绝对表达差累积分布。全部输入/输出为相对路径；每个cell附中英文说明。PNG/PDF和汇总表达、PCA坐标、基因归属及表达变化表均可重新生成。

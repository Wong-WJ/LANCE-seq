# LANCE-seq Figure 5

Reproduce final **Figure 5B–I and supplementary Figures S5–S9** as 19 PNG/PDF pairs.

## Files and run

- `LANCE-seq_Figure5.ipynb`: full analysis with short Chinese/English explanations for every cell.
- `figure5_utils.py`: reusable normalization, plotting and export helpers.
- `data/`: frozen GO/KEGG mouse annotations and source checksums.
- `panel_sources.csv`: final panel-to-source-notebook mapping.
- `outputs/`: PNG/PDF files; `outputs/tables/`: numerical results, gene lists and enrichment tables.

Place `LANCE8.h5ad` in `../data/`, then run `python run_figures.py --figure 5` from the repository root. Alternatively, run all notebook cells from this folder. All paths are relative, and no online enrichment service is needed. Tested packages are in `../requirements-figures45.txt`.

## Methods implemented

- **5B:** aggregate all genes across spots per sample, compute CPM and `log2(CPM+1)`, then plot mean expression A and MN2−MN1 change M. Detected genes with A>1 and M>1 or M<−1 define the up/down sets. These are sample-level descriptive changes.
- **5C/S5:** intersect the 5B-eligible genes with the frozen mouse annotations to define each database-specific background. Test terms of size 10–500 with a one-sided hypergeometric test; apply BH to all tested terms, including zero overlaps, separately by database and direction. Keep the final displayed term selections and a common bubble-area scale. Complete results, backgrounds, queries and displayed terms are saved as CSVs.
- **5D/E, S6/S7 and 5I:** per-spot normalization to 10,000 followed by natural-log `log1p`; native spatial coordinates. Each gene uses common color limits across compared samples. S6 section maps share the MN1/MN2 color limits of 5D. Boxes show median, IQR and 1.5-IQR whiskers; outlier symbols are hidden.
- **5F/S9:** four CV and four PV core genes. Microneedle scores use the equally weighted MN1/MN5 healthy reference, including within- and between-sample variance. Section scores use the CON3 mean/SD. The mean standardized core-gene expression defines each program score. Spatial colors use pooled 1st–99th percentiles.
- **S8:** within-sample core-gene z scores, mean CV minus mean PV, stable rank tertiles, and five held-out marker medians/Spearman correlations in healthy samples.
- **5G:** retain the revised depth, detection and mean-count matching procedure, including stochastic rounding of fractional values. Sample 30 APAP spots per depth stratum across 20 strata and match healthy spots with sufficient depth, for 300 draws. Retain the original two random streams (20261121 and 22061121). Plot median relative-expression profiles for held-out Glul, Cyp2e1 and Ass1 across within-sample tertiles. These draws resample spots within the specified captures.
- **5H:** include genes detected in ≥10% of spots in both MN1 and MN2 with finite Moran statistics. Compare pooled `log2(CPM+1)` changes with global Moran's I changes calculated on per-spot log-normalized expression and symmetric binary 6-nearest-neighbor graphs. Retain all 11 final labels, including the original protein-coding selection rules for Abi1, Scd1 and Alb. Endpoint and spatial-weight/depth sensitivity checks are saved as tables.

MN1/MN2 is the same-mouse baseline/APAP 6 h pair. MN5/MN4 and CON3/CON2 are endpoint comparisons; microneedle and section measurements come from the corresponding healthy and APAP animals.

## Final panel mapping

| Final panel | Output stem (`.png` / `.pdf`) | Content |
|---|---|---|
| Fig. 5B | `Fig5B_Longitudinal_MA` | Sample-level longitudinal MA plot |
| Fig. 5C | `Fig5C_GO_Up` | Upregulated GO biological processes |
| Fig. 5D | `Fig5D_Injury_Spatial` | Injury marker spatial expression |
| Fig. 5E | `Fig5E_Injury_Distributions` | Injury marker spot distributions |
| Fig. 5F | `Fig5F_CV_Reference_Spatial` | Healthy-reference CV program maps |
| Fig. 5G | `Fig5G_Matched_Marker_Profiles` | Depth/detection/mean-matched marker profiles |
| Fig. 5H | `Fig5H_Expression_Moran_Changes` | Expression changes versus Moran’s I changes |
| Fig. 5I | `Fig5I_Focal_Gene_Spatial` | Ass1, Abi1 and Alb spatial expression |
| Fig. S5A | `FigS5A_GO_Down` | Downregulated GO biological processes |
| Fig. S5B | `FigS5B_KEGG_Up` | Upregulated KEGG pathways |
| Fig. S5C | `FigS5C_KEGG_Down` | Downregulated KEGG pathways |
| Fig. S6A–B | `FigS6_Section_Injury_Markers` | Section injury marker maps and distributions |
| Fig. S7A | `FigS7A_Zonation_Marker_Spatial` | Zonation marker spatial expression |
| Fig. S7B | `FigS7B_Zonation_Marker_Distributions` | Zonation marker spot distributions |
| Fig. S8A | `FigS8A_Relative_Zonation_Spatial` | Within-sample relative zonation maps |
| Fig. S8B | `FigS8B_Heldout_Validation` | Five held-out markers in healthy samples |
| Fig. S9A | `FigS9A_CV_Reference_Distributions` | CV reference score distributions |
| Fig. S9B | `FigS9B_PV_Reference_Distributions` | PV reference score distributions |
| Fig. S9C | `FigS9C_PV_Reference_Spatial` | PV reference score maps |

Figure 5A is the experiment schematic. Final supplementary Figure S10 contains cluster maps and is outside this notebook's scope. Output names use the final figure numbering.

## Annotation provenance

`data/GO_mouse_annotation.json.gz` contains the exact GO Biological Process snapshot used in the revised figures (MSigDB M5, 2026.2.Mm). `data/KEGG_mouse_annotation.json.gz` contains the corresponding KEGG mouse snapshot retrieved on 2026-11-20. `annotation_sources.json` records original source URLs and checksums. The notebook validates the uncompressed snapshot checksum before use. These snapshots replace the old cached enrichment-result CSVs.

## 中文说明

本notebook按最终版生成图5 B–I及补充图S5–S9。原“新S9”中的CV空间图对应图5F，“S5.5”对应补充图S6，“S9-1”对应补充图S9C。图5B/G/H/I采用修改后的计算，保留随机种子、抽样次数、阈值和基因筛选。运行输出19组PNG/PDF，以及完整富集结果、查询基因、背景基因、重采样和空间统计CSV。最终补充图S10是聚类图，不使用修改目录中的旧同名自相关柱状图。

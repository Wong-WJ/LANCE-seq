# LANCE-seq Figure 5 analysis

## Files

- `LANCE-seq_Figure5.ipynb` — complete Figure 5 analysis notebook.
- `figure5_utils.py` — small reusable plotting and mathematical helpers; all core scientific calculations remain visible in the notebook.
- `data/` — GO and KEGG enrichment tables used by the notebook.
## Inputs

- AnnData path: set `INPUT_H5AD` in the configuration cell.
- Local GO/KEGG table directory: set `ENRICHMENT_DIR` in the configuration cell.
- Figure destination: set `OUTPUT_DIR` in the configuration cell.

The h5ad must contain count-like `X`, `obs['batch']`, and `obsm['spatial']`.

## Run

Open `LANCE-seq_Figure5.ipynb`, restart the kernel, and run all cells from top to bottom. The notebook recreates the output directory and produces exactly the 18 PNG files listed below.

## Key methods

- MN2 versus MN1 differential expression: Scanpy Wilcoxon on per-spot 10,000-count normalization followed by log1p; displayed thresholds are `|log2FC| > 1` and adjusted `P < 0.05`.
- Relative zonation: within-sample gene z-scores; mean CV-core score minus mean PV-core score; stable rank tertiles define Relative PV-like, Intermediate, and Relative CV-like states.
- Held-out validation: `Cyp2e1`, `Glul`, `Cyp7a1`, `Ass1`, and `Alb` are excluded from score construction and assessed by stratum medians and Spearman correlation.
- Healthy-reference zonation programs: equally weighted MN1/MN5 healthy reference with both within-sample and between-sample variance.
- Zone effects: sample×stratum count aggregation, CPM, `log2(CPM+1)`, and APAP−baseline differences without replicate-level inference.
- Figure 5h modules: each gene is standardized with fixed MN1-spot mean and SD before module averaging in MN1/MN2.
- Supportive concordance: Spearman comparison of MN1→MN2 effects with MN5→MN4 and CON3→CON2 effects.
- Spatial autocorrelation: symmetric binary 6-nearest-neighbor weights, Global Moran’s I, 999 permutations, and within-sample BH-FDR.

## Outputs and manuscript panels

| Output file | Manuscript panel |
|---|---|
| `Fig5B_MN2_vs_MN1_Volcano.png` | Fig. 5b |
| `Fig5C1_GO_Up.png` | Fig. 5c |
| `Fig5D_APAP_Injury_Gene_Spatial_Maps_fixed.png` | Fig. 5d |
| `FigS5B_APAP_Marker_Expression_Distribution_fixed.png` | Fig. 5e |
| `Fig5F_Relative_Zonation_3State_Spatial_Maps_fixed.png` | Fig. 5f |
| `Fig5G_ZoneStratified_Transcriptional_Effect_Sizes.png` | Fig. 5g |
| `Fig5H_APAP_Module_Response_Across_Relative_Zonation_Strata.png` | Fig. 5h |
| `Fig5C2_GO_Down.png` | Supplementary Fig. 5a |
| `FigS5A1_KEGG_Up.png` | Supplementary Fig. 5b |
| `FigS5A2_KEGG_Down.png` | Supplementary Fig. 5c |
| `Fig5E_Liver_Zonation_Marker_Spatial_Maps.png` | Extended Data Fig. 1a |
| `FigS5C_Zonation_Marker_Expression_Distribution.png` | Extended Data Fig. 1b |
| `FigS5D_Relative_Zonation_Continuous_Spatial_Maps_fixed.png` | Extended Data Fig. 2a |
| `FigS5E_Healthy_HeldOut_Zonation_Validation.png` | Extended Data Fig. 2b |
| `FigS5G_HealthyReference_Absolute_Zonation_Program_Remodeling.png` | Extended Data Fig. 3a,b |
| `FigS5H_Absolute_CV_Associated_Program_Spatial_Maps.png` | Extended Data Fig. 3c |
| `FigS5J_Spatial_Autocorrelation_MoransI.png` | Supplementary Fig. 6 |
| `FigS5I_ZoneStratified_Supportive_Concordance.png` | Supplementary Fig. 7 |

## Dependencies

Python 3.10, scanpy 1.11.5, anndata 0.11.4, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3, matplotlib 3.10.8, seaborn 0.13.2, statsmodels, and Jupyter/nbformat.

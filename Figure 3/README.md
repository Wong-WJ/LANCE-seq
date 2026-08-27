# LANCE-seq Figure 3 cell-type deconvolution

This folder contains the cell2location analysis used for the mouse brain cell-type deconvolution shown in **Figure 3h** of the LANCE-seq manuscript.

The notebook only covers the deconvolution analysis. Other Figure 3 analyses, including clustering, anatomical annotation, marker-gene analysis, and MN/Section comparisons, are not included here.

## File

- `Figure 3 c2loc_mouse_brain.ipynb` — cell2location reference-model training, spatial deconvolution, and visualization

## Reference dataset

Cell-type expression signatures are derived from the mouse whole-brain single-nucleus RNA-seq dataset **E-MTAB-11115**.

The notebook uses:

- `cell_annotation.csv` for cell annotations
- sample-level 10x Genomics `filtered_feature_bc_matrix.h5` files
- `annotation_1` as the cell-type label
- `sample` as the batch variable

Cells annotated as `Unk_1`, `Unk_2`, `LowQ_1`, or `LowQ_2` are excluded before reference-model fitting.

## Spatial input

The spatial transcriptomic input is the filtered LANCE-seq mouse brain AnnData object:

`filtered_adata.h5ad`

The AnnData object must contain the spatial coordinates in:

`adata.obsm["spatial"]`

Input paths in the notebook should be updated to match the local data location before running.

## cell2location workflow

Reference signatures are estimated using `cell2location.models.RegressionModel`.

Main settings used in the notebook:

- reference model training: `max_epochs=250`
- posterior sampling: `num_samples=1000`
- posterior batch size: `2500`

The resulting reference signatures are exported as `inf_aver.csv` and used for spatial mapping.

Spatial deconvolution is performed with `cell2location.models.Cell2location` using:

- `N_cells_per_location=25`
- `detection_alpha=20`
- `detection_mean_per_sample=False`
- spatial model training: `max_epochs=8000`
- posterior sampling: `num_samples=1000`

The estimated cell-type abundances are stored in:

`adata.obsm["q05_cell_abundance_w_sf"]`

and the deconvolved spatial AnnData object is saved as `sp.h5ad`.

## Figure 3h

Twenty-four representative neuronal and glial subpopulations are selected for the final visualization:

- `Ext_Hpc_CA1`
- `Ext_Hpc_CA2`
- `Ext_Hpc_CA3`
- `Ext_Hpc_DG1`
- `Ext_L5_1`
- `Ext_L5_2`
- `Ext_L6`
- `Ext_L6B`
- `Ext_L23`
- `Ext_L25`
- `Ext_L56`
- `Ext_Amy_2`
- `Ext_Pir`
- `Ext_Thal_1`
- `Ext_Thal_2`
- `Inh_1`
- `Inh_4`
- `Inh_5`
- `Inh_6`
- `Inh_Meis2_3`
- `Inh_Meis2_4`
- `Oligo_1`
- `Oligo_2`
- `Micro`

Their inferred abundances are mapped back to the original LANCE-seq spatial coordinates for visualization.

## Run

Open `Figure 3 c2loc_mouse_brain.ipynb`, update the reference-data and spatial-data paths, and run the notebook from top to bottom.

Main dependencies:

`cell2location`, `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`, `scipy`

A CUDA-enabled GPU is recommended for model training.

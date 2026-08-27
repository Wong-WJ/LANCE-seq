# LANCE-seq

Analysis code associated with the LANCE-seq manuscript.

LANCE-seq (Live-tissue Analysis via Needle-encoded Capture Extraction sequencing) is a spatial transcriptomic platform based on spatially encoded microneedle arrays for transcriptome profiling of living tissues.

## Repository structure

- `FASTQ_preprocessing/` — custom Read 1 barcode/UMI extraction and Read 2 trimming scripts used before STARsolo processing.
- `Figure 3/` — mouse brain cell-type deconvolution using cell2location (Fig. 3h)
- `Figure 4/` — pseudobulk transcriptomic analyses for the liver longitudinal-sampling experiment (Fig. 4f–i)
- `Figure 5/` — longitudinal APAP-response, liver zonation, functional-module, and spatial-autocorrelation analyses (Fig. 5 and Supplementary Figs. S5–S10)

Each folder contains its own README with analysis-specific inputs, methods, and output-panel mapping.

## Data availability

Raw sequencing data generated in this study have been deposited in the Genome Sequence Archive (GSA), National Genomics Data Center. The accession number will be added here after assignment.

Processed spatial transcriptomic data, including expression matrices, spot metadata, and spatial coordinates in AnnData format, are available at Zenodo under DOI https://doi.org/10.5281/zenodo.22121908.

The mouse whole-brain single-nucleus RNA-seq reference dataset used for cell-type deconvolution is available from ArrayExpress under accession `E-MTAB-11115`.

## Requirements

The analysis notebooks primarily use Python packages including:

- Scanpy
- AnnData
- NumPy
- pandas
- SciPy
- Matplotlib
- seaborn
- scikit-learn
- statsmodels
- cell2location (Figure 3)

See the README within each figure folder for analysis-specific requirements.

## Usage

Update the input paths defined near the beginning of each notebook, then run the notebook from top to bottom.

Large sequencing and processed-data files are not stored in this repository.

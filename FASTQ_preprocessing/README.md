# LANCE-seq FASTQ preprocessing

This folder contains the custom FASTQ preprocessing scripts used before STARsolo alignment and spatial barcode assignment in the LANCE-seq workflow.

## Files

- `fastq_process_docmap.py` — extracts the spatial barcode and UMI sequence from Read 1.
- `fastq_process_TSO.py` — removes the first 30 bases from Read 2 and retains the transcript-derived sequence used for alignment.

## Requirements

- Python 3
- Biopython

Install Biopython if needed:

```bash
pip install biopython
```

## Read 1 preprocessing

`fastq_process_docmap.py` reads a gzipped FASTQ file and reconstructs a 28-bp barcode/UMI read as:

```text
BC2 (8 bp) + BC1 (8 bp) + UMI (12 bp)
```

The sequence segments extracted by the script are:

```text
Python slice   1-based read positions   Length
seq[15:23]     16–23                    8 bp
seq[53:61]     54–61                    8 bp
seq[61:73]     62–73                    12 bp
```

The corresponding quality-score segments are concatenated in the same order.

Example:

```bash
python fastq_process_docmap.py \
    -i sample_R1.fastq.gz \
    -o sample_R1_processed.fastq
```

## Read 2 preprocessing

`fastq_process_TSO.py` removes the first 30 bases of Read 2 and retains bases 31–150:

```python
seq[30:150]
```

The corresponding quality scores are trimmed identically.

Example:

```bash
python fastq_process_TSO.py \
    -i sample_R2.fastq.gz \
    -o sample_R2_processed.fastq
```

## Downstream processing

The processed paired-end FASTQ files are used for downstream alignment and gene-level quantification with STARsolo. Spatial barcode assignment is performed using the predefined LANCE-seq barcode whitelist.

The nucleotide sequences of the spatial barcodes, linkers and related oligonucleotides are provided in the Supplementary Information of the associated manuscript.

## Notes

- The scripts expect gzipped FASTQ input files.
- Output FASTQ files are written as uncompressed text FASTQ files.
- The preprocessing scheme was designed for the 150-bp paired-end sequencing configuration used in this study.

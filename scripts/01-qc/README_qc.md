# QC scripts

## 01_extract_counts.py

What it does:
- Extracts the expression table from the GEO series matrix file.
- Uses ID_REF as the gene identifier and sums counts across duplicate IDs.
- Writes a tab-separated counts matrix with GSM sample IDs as columns.

How to run:
python3 scripts/01-qc/01_extract_counts.py

Expected output:
- data/processed/counts.tsv

## 03_qc_plots.py

What it does:
- Creates QC summaries and plots from counts and metadata.
- Filters samples to main tumor vs normal comparisons.
- Produces library size and PCA plots.

How to run:
python3 scripts/01-qc/03_qc_plots.py

Expected outputs:
- figures/qc/library_size.png
- figures/qc/pca.png
- results/qc/qc_summary.tsv

## 02_merge_htseq_counts.py

What it does:
- Merges HTSeq count files from data/raw/GSE306117_RAW/.
- Removes rows where gene_id starts with "__".
- Outer-joins samples on gene_id and fills missing values with 0.

How to run:
python3 scripts/01-qc/02_merge_htseq_counts.py

Expected output:
- data/processed/counts.tsv

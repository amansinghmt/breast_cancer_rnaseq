#!/usr/bin/env bash
set -euo pipefail

echo "[v2] 1/7 metadata"
python3 scripts/00-metadata/01_make_metadata.py

echo "[v2] 2/7 merge counts"
python3 scripts/01-qc/02_merge_htseq_counts.py

echo "[v2] 3/7 qc summary"
python3 scripts/01-qc/03_qc_plots.py

echo "[v2] 4/7 sample manifest"
python3 scripts/00-metadata/02_make_sample_manifest_v2.py

echo "[v2] 5/7 paired DE"
Rscript scripts/02-de/01_deseq2_paired_v2.R

echo "[v2] 6/7 enrichment"
Rscript scripts/03-pathways/01_enrichment_paired_v2.R

echo "[v2] 7/7 publication figures"
Rscript scripts/04-figures/01_publication_figures_paired_v2.R

echo "[v2] completed"

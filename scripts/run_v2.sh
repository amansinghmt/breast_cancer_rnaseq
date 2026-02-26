#!/usr/bin/env bash
set -euo pipefail

echo "[v2] 1/6 metadata"
python3 scripts/00-metadata/01_make_metadata.py

echo "[v2] 2/6 merge counts"
python3 scripts/01-qc/02_merge_htseq_counts.py

echo "[v2] 3/6 qc summary"
python3 scripts/01-qc/03_qc_plots.py

echo "[v2] 4/6 sample manifest"
python3 scripts/00-metadata/02_make_sample_manifest_v2.py

echo "[v2] 5/6 paired DE"
Rscript scripts/02-de/01_deseq2_paired_v2.R

echo "[v2] 6/6 enrichment"
Rscript scripts/03-pathways/01_enrichment_paired_v2.R

echo "[v2] completed"

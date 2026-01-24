# Differential expression (DE)

Comparison:
- Tumor vs Normal (main tumor vs normal samples only)

Design formula:
- ~ lfs_status + condition

How to run:
Rscript scripts/02-de/01_deseq2_tumor_vs_normal.R

Outputs:
- results/differential_expression/deseq2_results.tsv
- figures/de/ma_plot.png

# Breast cancer RNA-seq analysis (GSE306117)

Differential expression and pathway analysis of breast cancer RNA-seq data from GEO GSE306117.
Primary comparison: Tumor vs Normal, adjusting for LFS status.
Includes reproducible metadata generation, QC, DE, and enrichment workflows.

## Repository structure
- data/raw (gitignored)
- data/processed (gitignored)
- data/metadata (tracked)
- scripts/00-metadata
- scripts/01-qc
- scripts/02-de
- scripts/03-pathways
- figures/qc, figures/de, figures/enrichment
- results/differential_expression, results/enrichment

## Environment setup
A) Python venv
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.lock
```

B) R package restore (renv)
```r
install.packages("renv")
renv::restore()
```

## Reproduce the analysis
1) metadata: `python scripts/00-metadata/01_make_metadata.py`
2) merge counts (htseq): `python scripts/01-qc/02_merge_htseq_counts.py`
3) qc plots: `python scripts/01-qc/03_qc_plots.py`
4) deseq2: `Rscript scripts/02-de/01_deseq2_tumor_vs_normal.R`
5) volcano + heatmap: `Rscript scripts/02-de/02_volcano_heatmap.R`
6) enrichment: `Rscript scripts/03-pathways/01_enrichment_hallmark_go.R`

## V2 Paired Pipeline (recommended)
Run from scratch:
```bash
python3 scripts/00-metadata/01_make_metadata.py
python3 scripts/01-qc/02_merge_htseq_counts.py
python3 scripts/01-qc/03_qc_plots.py
python3 scripts/00-metadata/02_make_sample_manifest_v2.py
Rscript scripts/02-de/01_deseq2_paired_v2.R
Rscript scripts/03-pathways/01_enrichment_paired_v2.R
Rscript scripts/04-figures/01_publication_figures_paired_v2.R
```

One command:
```bash
./scripts/run_v2.sh
```

V2 outputs are written under:
- `results_v2/differential_expression`
- `results_v2/enrichment`
- `figures_v2/de`
- `figures_v2/final`
- `results_v2/figures`

## Key outputs
- figures/qc/library_size.png
- figures/qc/pca.png
- figures/de/ma_plot.png
- figures/de/volcano.png
- figures/de/heatmap_top50.png
- figures/enrichment/hallmark_gsea_top10.png
- figures/enrichment/go_barplot_top15.png
- results/differential_expression/deseq2_results.tsv
- results/differential_expression/top_genes.tsv
- results/enrichment/hallmark_gsea.tsv
- results/enrichment/go_enrich.tsv

V2 publication-ready figure set:
- figures_v2/final/Fig1A_qc_library_size_pairs.png
- figures_v2/final/Fig1B_qc_pca_paired.png
- figures_v2/final/Fig2A_de_ma_plot_refined.png
- figures_v2/final/Fig2B_de_volcano_refined.png
- figures_v2/final/Fig2C_de_heatmap_top40.png
- figures_v2/final/Fig3A_bio_hallmark_gsea.png
- figures_v2/final/Fig3B_bio_go_bp_dotplot.png
- results_v2/figures/final_figure_manifest_paired_v2.tsv

## Notes / caveats
- Raw GEO series matrix had no expression rows; counts were obtained from GEO supplementary HTSeq count files and merged.
- data/raw and data/processed are gitignored due to size; regenerate with scripts.
- GEO landing page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306117

## Citation
Data source: NCBI GEO accession GSE306117.

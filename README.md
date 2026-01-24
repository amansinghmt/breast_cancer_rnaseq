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
A) Python venv (already used)
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

B) R packages installation
```r
install.packages(readr,dplyr,ggplot2,pheatmap,ggrepel,msigdbr,fgsea,BiocManager)
BiocManager::install(DESeq2, clusterProfiler, org.Hs.eg.db, AnnotationDbi)
```

## Reproduce the analysis
1) metadata: `python scripts/00-metadata/01_make_metadata.py`
2) merge counts (htseq): `python scripts/01-qc/02_merge_htseq_counts.py`
3) qc plots: `python scripts/01-qc/03_qc_plots.py`
4) deseq2: `Rscript scripts/02-de/01_deseq2_tumor_vs_normal.R`
5) volcano + heatmap: `Rscript scripts/02-de/02_volcano_heatmap.R`
6) enrichment: `Rscript scripts/03-pathways/01_enrichment_hallmark_go.R`

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

## Notes / caveats
- Raw GEO series matrix had no expression rows; counts were obtained from GEO supplementary HTSeq count files and merged.
- data/raw and data/processed are gitignored due to size; regenerate with scripts.
- GEO landing page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306117

## Citation
Data source: NCBI GEO accession GSE306117.

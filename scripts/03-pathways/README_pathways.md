# Pathway enrichment

Inputs:
- results/differential_expression/deseq2_results.tsv

Methods:
- Hallmark GSEA with fgsea (MSigDB Hallmark gene sets).
- GO Biological Process over-representation with clusterProfiler::enrichGO.

How to run:
Rscript scripts/03-pathways/01_enrichment_hallmark_go.R

Outputs:
- results/enrichment/hallmark_gsea.tsv
- figures/enrichment/hallmark_gsea_top10.png
- results/enrichment/go_enrich.tsv
- figures/enrichment/go_barplot_top15.png

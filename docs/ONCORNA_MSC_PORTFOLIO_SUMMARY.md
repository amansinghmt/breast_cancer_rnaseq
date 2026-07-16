# OncoRNA: MSc Portfolio Summary

## Scientific question

Which gene-expression differences and predefined biological programs are associated with breast Tumor versus matched Normal tissue in GEO GSE306117?

## Dataset and methods

The project processes 74 public HTSeq samples and selects 21 patients with one QC-passing Tumor and one matched Normal sample each. Raw counts are modeled in DESeq2 using `~ patient_id + condition_main`, with Normal as reference. It uses Wald tests, Benjamini-Hochberg FDR, log2-fold-change shrinkage, VST-based exploratory PCA/heatmaps, Hallmark preranked GSEA, and combined plus directional GO Biological Process ORA.

## Outputs and reproducibility

The repository provides validated manifests, complete DE/enrichment tables, seven PNG/PDF figures, mapping diagnostics, robustness summaries, checksums, session information, timestamped logs, a read-only Streamlit dashboard, and regression tests. The maintained entry point is `bash scripts/run_v2.sh`.

## Main result

Within this 21-patient cohort, 6,315 genes have BH-adjusted p-value below 0.05; 1,636 also have absolute shrunken log2 fold change at least 1. Hallmark and directional GO analyses associate Tumor-higher genes with cell-cycle/division-related programs, while Normal-higher genes include circulation, muscle, extracellular-matrix and other tissue-context annotations.

## Limitations

This is one bulk-tissue cohort without independent replication, protein measurement, functional experiments, clinical outcomes, or patient-level prediction. Mapping is incomplete, thresholds affect selected-list analyses, and cell composition may explain part of the observed signal.

## Skills demonstrated

R, Python, Bash, RNA-seq count handling, paired experimental design, DESeq2, multiple testing, enrichment analysis, scientific visualization, reproducible pipelines, output validation, documentation, and claim-boundary review.

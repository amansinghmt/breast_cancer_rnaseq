# Project Summary v2 — Paired Tumor vs Normal Breast Cancer RNA-seq Analysis

## 1. Research question

This project asks a focused transcriptomics question: how do gene-expression profiles differ between breast tumor tissue and matched normal tissue from the same patients? The maintained workflow is designed to detect differential expression at the gene level and then summarize those changes at the pathway level.

The practical goal is to produce a reproducible analysis that can support careful interpretation of tumor-versus-normal transcriptional differences. In this repository, that means running one maintained pipeline (`bash scripts/run_v2.sh`) and generating validated outputs under `results_v2/` and `figures_v2/`.

## 2. Dataset and study design

The analysis is based on GEO accession **GSE306117**. The repository workflow parses the GEO series information and sample-level metadata, then constructs a paired analysis cohort through an explicit manifest process.

The maintained v2 workflow is a **paired tumor vs normal design**. In the current repository outputs, `results_v2/metadata/paired_manifest.tsv` currently includes 21 paired patients (42 `include_paired` samples), where each included patient contributes one Tumor and one Normal sample. This pairing rule is enforced in pipeline validation logic rather than assumed informally.

The key design implication is that each tumor sample is interpreted relative to its own matched normal counterpart from the same patient. This is more controlled than a simple unpaired group comparison and is well aligned with the biological question of within-patient tumor-associated expression change.

## 3. Why paired analysis was used

In human transcriptomic datasets, between-patient variability can be large because of genetics, microenvironment, sample handling, and other factors unrelated to the main contrast of interest. If those differences are not modeled, they can obscure or confound tumor-versus-normal signals.

A paired framework addresses this by comparing conditions within patient. Conceptually, each patient acts as their own baseline. In this repository, that idea is implemented directly in the maintained differential expression model and in sample-manifest checks that require one Tumor and one Normal sample per included patient.

For a learning-focused project, this is also an important methodological lesson: study design decisions (paired vs unpaired) often have as much impact on interpretability as downstream plotting or enrichment choices.

## 4. Differential expression workflow

The maintained differential expression workflow begins with structured data preparation. `scripts/run_v2.sh` executes metadata and QC scripts, builds a sample manifest, and copies the validated paired cohort definition to:

- `results_v2/metadata/paired_manifest.tsv`

From there, paired differential expression is run by:

- `scripts/02-de/01_deseq2_paired_v2.R`

The maintained model formula is:

- `~ patient_id + condition_main`

This model captures patient-specific baseline effects (`patient_id`) while testing the tumor-versus-normal condition contrast (`condition_main`). The canonical DE outputs are written to:

- `results_v2/deseq2/deseq2_paired_v2_results.tsv`
- `results_v2/deseq2/deseq2_paired_v2_samples_used.tsv`
- `results_v2/deseq2/sessionInfo_paired_v2.txt`

At interpretation time, two common columns are central:

- **log2 fold change (log2FC):** estimated effect size for the modeled Tumor-versus-Normal contrast (direction and magnitude of change)
- **adjusted p-value (FDR / `padj`):** controls false discoveries across many tested genes

The repository validation steps also check that required columns exist and that the `samples_used` table matches the paired manifest exactly. This helps ensure that reported differential expression reflects the intended cohort and design.

## 5. Pathway enrichment workflow

Gene-level differential expression often produces long result tables that are hard to interpret directly. To move from individual genes to biological themes, the maintained pipeline performs two complementary enrichment analyses:

- **Hallmark GSEA** (rank-based enrichment)
- **GO Biological Process ORA** (over-representation analysis)

These analyses are executed through:

- `scripts/03-pathways/01_enrichment_paired_v2.R`

Canonical enrichment outputs are:

- `results_v2/enrichment/hallmark_gsea_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_paired_v2.tsv`
- `results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt`

From a research-communication perspective, this step is useful because pathway-level summaries can reveal coherent biological programs that would be difficult to infer from a ranked gene list alone. At the same time, enrichment remains interpretive and should be read as evidence of pattern-level association, not direct proof of mechanism.

## 6. Figure overview

The maintained workflow produces a fixed, publication-style figure contract in `figures_v2/final/`, with exactly seven panels:

- `F01.png`
- `F02.png`
- `F03.png`
- `F04.png`
- `F05.png`
- `F06.png`
- `F07.png`

At a high level, these panels cover:

- **F01:** paired-cohort library size quality-control summary
- **F02:** paired-cohort PCA quality-control summary
- **F03:** MA plot summarizing differential expression results
- **F04:** volcano plot summarizing effect size and statistical significance
- **F05:** heatmap view of differential expression patterns
- **F06:** Hallmark enrichment results summary
- **F07:** GO Biological Process enrichment results summary

The workflow enforces exact naming (`F01` through `F07`) and archives non-canonical figure files out of the final folder. This creates a stable figure interface for reporting and downstream documentation.

## 7. Biological interpretation

This workflow is best viewed as an interpretation framework for tumor-vs-normal transcriptomic contrast in a paired setting. It supports three linked levels of reasoning:

1. **Gene level:** identify transcripts with statistically supported differences between conditions.
2. **Pattern level:** examine whether sets of related genes shift together.
3. **Pathway level:** summarize changes into biological programs that are easier to discuss and compare.

In many cancer transcriptomics studies, themes related to proliferation and cell-cycle-associated biology are frequently relevant. In this project, pathway analysis is included so that such patterns can be evaluated systematically if they are supported by the dataset-specific results. The intent is cautious interpretation rather than broad claims of novelty.

Because the analysis is paired, any observed tumor-normal differences are interpreted within patient context, which can improve confidence that detected contrasts are related to condition differences rather than only cross-patient heterogeneity.

## 8. Reproducibility and project engineering

A key strength of this repository is that it treats the analysis as an engineered workflow, not just a collection of plots. The maintained entrypoint is explicit:

```bash
bash scripts/run_v2.sh
```

Reproducibility features include:

- `renv.lock` plus `renv::restore()` for R dependency restoration
- deterministic execution controls in the maintained pipeline (single pipeline seed and fixed execution order)
- strict output validation at multiple stages
- canonical output directories (`results_v2/deseq2/`, `results_v2/enrichment/`, `results_v2/metadata/`)
- final figure contract enforcement (`figures_v2/final/F01.png` to `F07.png`)
- run logging in `results_v2/logs/`
- checksum manifests in `results_v2/fig_manifest.tsv` and `results_v2/output_manifest.tsv`

A smoke test is also available:

```bash
bash tests/smoke_test_v2.sh
```

This test runs the maintained workflow and checks required DE outputs, enrichment outputs, manifests, exact figure set, and log presence. That gives a practical verification layer for collaborators and future maintenance.

## 9. Limitations

This project has important analytical limits that should be stated clearly.

First, this is **bulk RNA-seq**, so each sample reflects mixed cell populations rather than single-cell resolution. Differential expression therefore captures composite tissue-level signals.

Second, the maintained workflow is currently centered on **one public dataset** (GSE306117). Findings and interpretation should be considered in that dataset context unless replicated externally.

Third, transcript abundance differences do not by themselves establish mechanism or causality. They provide evidence of association, not direct functional proof.

Fourth, pathway enrichment methods (including GSEA and ORA) are interpretive tools that depend on gene-set definitions and statistical thresholds. They are valuable summaries, but not definitive biological proof in isolation.

## 10. Project value

As a learning and portfolio artifact, this project demonstrates core bioinformatics skills in a coherent way: RNA-seq workflow construction, paired differential-expression analysis, pathway interpretation, reproducible pipeline engineering, and technical communication. It is useful for MSc review because it shows both analytical reasoning and practical reproducibility habits in a real repository context. A technical, step-by-step workflow reference is provided in `docs/workflow_v2.md`.

# OncoRNA Code Learning Map

Study these files in order. The goal is to explain the data contract and scientific decisions, not memorize every line.

## 1. Metadata classification

**File:** `scripts/00-metadata/01_make_metadata.py`  
**Difficulty:** Beginner

- Purpose: convert GEO series metadata text into one structured row per sample.
- Input: `data/raw/GSE306117_series_matrix.txt`.
- Output: `data/metadata/metadata.tsv`.
- Important functions: `parse_geo_values`, `tissue_type_from_tissue`, status-classification helpers.
- Concepts: text parsing, categorical labels, precedence in rule-based classification.
- You should explain: why surrounding/contralateral tissue is not silently treated as Normal.
- Safe exercise: add a dry-run summary that prints category counts without changing output.

## 2. Paired-cohort construction

**File:** `scripts/00-metadata/02_make_sample_manifest_v2.py`  
**Difficulty:** Intermediate

- Purpose: apply QC and matching rules and record every exclusion reason.
- Inputs: metadata and `results/qc/qc_summary.tsv`.
- Output: `data/metadata/sample_manifest.tsv`.
- Important blocks: `QC_MIN`, per-row exclusion reasons, `by_patient`, highest-library candidate selection.
- Concepts: deterministic filtering, grouping, constrained selection, audit trails.
- You should explain: one-million-count rule, exactly one pair per patient, and limitations of highest-depth selection.
- Safe exercise: print how included-patient count would change at a hypothetical threshold without writing the manifest.

## 3. Count-matrix generation

**File:** `scripts/01-qc/02_merge_htseq_counts.py`  
**Difficulty:** Beginner

- Purpose: merge compressed per-sample HTSeq files into a gene-by-sample matrix.
- Input: `data/raw/GSE306117_RAW/*.txt.gz`.
- Output: `data/processed/counts.tsv`.
- Important functions: `list_count_files`, `extract_sample_id`, `load_counts`, `main`.
- Concepts: regular expressions, outer joins, missing gene counts, integer validation.
- You should explain: why `__*` technical rows are removed and missing genes become zero.
- Safe exercise: add a read-only count of technical rows removed per sample.

## 4. QC and transformation

**File:** `scripts/01-qc/03_qc_plots.py`  
**Difficulty:** Intermediate

- Purpose: count-depth summary and initial exploratory PCA.
- Inputs: counts and metadata.
- Outputs: QC summary and legacy QC plots.
- Important functions: `compute_log2_cpm`, `plot_library_sizes`, `plot_pca`.
- Concepts: library sums, CPM, log2 transform, variance selection, PCA.
- You should explain: why log2-CPM is acceptable for exploratory QC but not used as DESeq2's inferential input.
- Safe exercise: compare top 500 versus top 2,000 variable genes in a temporary PCA plot.

## 5. Paired differential expression

**File:** `scripts/02-de/01_deseq2_paired_v2.R`  
**Difficulty:** Advanced

- Purpose: fit the authoritative paired Tumor-versus-Normal model.
- Inputs: raw counts and paired manifest.
- Outputs: DE result, samples used, VST, diagnostics, size factors and session information.
- Important blocks: pair validation, count integrity, factor levels, `DESeqDataSetFromMatrix`, `DESeq`, `results`, `lfcShrink`, VST.
- Concepts: design matrix, negative-binomial GLM, dispersion, Wald test, BH FDR, independent filtering, Cook's distance, shrinkage.
- You should explain: `~ patient_id + condition_main`, Normal reference, positive effect direction, why `alpha=0.1` is retained for filter optimization while reporting uses `padj<0.05`.
- Safe exercise: read the diagnostics TSV and explain every metric; do not alter model parameters.

## 6. Enrichment analysis

**File:** `scripts/03-pathways/01_enrichment_paired_v2.R`  
**Difficulty:** Advanced

- Purpose: run full-rank Hallmark GSEA and selected-list GO BP ORA.
- Input: canonical DE table.
- Outputs: Hallmark table, complete GO table, representative GO table and mapping diagnostics.
- Important blocks: `universe_genes`, `sig_genes`, Wald-stat ranking, symbol duplicate handling, `fgseaMultilevel`, `enrichGO`, `simplify`.
- Concepts: rank statistics, null enrichment, NES, over-representation, hypergeometric reasoning, gene universe, annotation loss, semantic similarity.
- You should explain: why GSEA uses the full ranking, ORA uses 1,636 selected genes, and the GO denominator becomes 1,036/15,233 after annotation.
- Safe exercise: change only the number of representative terms displayed, never the raw GO table or significance calculation.

## 7. Robustness summaries

**File:** `scripts/05-reporting/01_build_robustness_summaries.R`  
**Difficulty:** Advanced

- Purpose: calculate predefined sensitivity and presentation tables without replacing canonical results.
- Inputs: canonical DE/enrichment tables, VST, counts and manifest.
- Outputs: threshold grid, prefilter sensitivity, LFC agreement, PCA outlier summary and annotated top genes.
- Concepts: sensitivity analysis, rank correlation, set overlap, robust distance, distinction between canonical and alternative analysis.
- You should explain: why threshold summaries are not threshold optimization and why the prefiltered model remains secondary.
- Safe exercise: add a displayed summary for `padj<0.1` clearly labelled exploratory, without changing canonical outputs.

## 8. Figure generation

**Files:** `scripts/04-figures/01_publication_figures_paired_v2.R` through `07_go_bp_dotplot_paired_v2.R`  
**Difficulty:** Intermediate

- Purpose: convert canonical tables into seven fixed communication panels and PDFs.
- Inputs: paired manifest, VST, DE and pathway tables.
- Outputs: descriptive figure files that become F01-F07.
- Concepts: log axes, effect/evidence categories, balanced ranking, z-scores, clustering, NES direction, GeneRatio.
- You should explain: each selection rule and what each visual cannot prove.
- Safe exercise: change colors or font sizes only, then verify that selected rows and metrics are unchanged.

## 9. Orchestration and validation

**File:** `scripts/run_v2.sh`  
**Difficulty:** Intermediate

- Purpose: execute the maintained workflow in order and reject incomplete/inconsistent outputs.
- Inputs: raw files, lockfiles and maintained scripts.
- Outputs: results, figures, session information, logs and checksum manifests.
- Important functions/blocks: `run_r_script`, `run_python`, paired-manifest validation, sample-used validation, exact figure contract, MD5 manifest generation.
- Concepts: strict shell mode, deterministic seeds, data contracts, fail-fast validation, reproducibility.
- You should explain: which outputs are authoritative and why legacy directories are excluded.
- Safe exercise: add a non-destructive required-column check to the smoke test.

## 10. Read-only dashboard

**File:** `dashboard/app.py`  
**Difficulty:** Intermediate

- Purpose: present existing canonical outputs without recomputation.
- Inputs: `results_v2` TSVs and `figures_v2/final/F01-F07.png`.
- Output: localhost Streamlit website.
- Important functions: `validate_required_files`, `read_tsv`, `metric_map`, `render_figure`.
- Concepts: presentation/data separation, graceful errors, caching, claim boundaries.
- You should explain: why no analysis result is written or altered by the dashboard.
- Safe exercise: add another explanatory expander using an existing canonical column.

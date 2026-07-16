#!/usr/bin/env Rscript

# Purpose:
#   Run paired differential expression for Tumor vs Normal using DESeq2 with
#   patient matching (21 paired patients, 42 samples in the included cohort).
#
# Inputs:
#   - data/metadata/sample_manifest.tsv
#   - data/processed/counts.tsv
#
# Outputs:
#   - results_v2/deseq2/deseq2_paired_v2_results.tsv
#   - results_v2/deseq2/deseq2_paired_v2_samples_used.tsv
#   - results_v2/deseq2/deseq2_paired_v2_vst.tsv
#   - results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv
#   - results_v2/deseq2/deseq2_paired_v2_size_factors.tsv
#   - results_v2/deseq2/sessionInfo_paired_v2.txt
#   - figures_v2/de/ma_plot_paired_v2.png
#
# Determinism / reproducibility:
#   - Fixed paired design formula: ~ patient_id + condition_main
#   - Fixed contrast direction: Tumor vs Normal
#   - Fixed coefficient selection rule: grep("^condition_main_", resultsNames(dds))
#   - Fixed LFC shrinkage method: lfcShrink(type = "normal")
#   - Exact sample order follows manifest include_paired ordering
#
# How to run:
#   Rscript --vanilla scripts/02-de/01_deseq2_paired_v2.R
#
# Assumptions:
#   - include_paired == TRUE corresponds to one Normal and one Tumor per patient.
#   - counts.tsv contains raw non-negative integer-like counts with unique gene_id.

#### Methods explained ####
# DESeq2 models RNA-seq counts with a negative-binomial GLM:
#   1) size factors normalize library depth,
#   2) dispersions are estimated/shrunk,
#   3) Wald tests evaluate coefficients in the design matrix.
#
# Paired design:
#   design = ~ patient_id + condition_main
# This controls patient-specific baseline expression and estimates the within-patient
# Tumor-vs-Normal effect.
#
# condition_main coefficient:
#   resultsNames(dds) contains model coefficients. We select the unique term matching
#   "^condition_main_" and shrink that effect for stable ranking/visualization.
#
# Shrinkage method:
#   lfcShrink(type = "normal") is used as a fixed, deterministic method with no optional
#   branch logic.

#### 1) Load packages + sanity checks ####
required_packages <- c("DESeq2", "readr", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      "\nRestore with renv::restore() before running this script."
    )
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
})

#### 2) Read inputs ####
paired_manifest_path <- "data/metadata/sample_manifest.tsv"
counts_path <- "data/processed/counts.tsv"

de_results_path <- "results_v2/deseq2/deseq2_paired_v2_results.tsv"
samples_used_path <- "results_v2/deseq2/deseq2_paired_v2_samples_used.tsv"
session_info_path <- "results_v2/deseq2/sessionInfo_paired_v2.txt"
vst_path <- "results_v2/deseq2/deseq2_paired_v2_vst.tsv"
diagnostics_path <- "results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv"
size_factors_path <- "results_v2/deseq2/deseq2_paired_v2_size_factors.tsv"
ma_plot_path <- "figures_v2/de/ma_plot_paired_v2.png"

required_manifest_cols <- c(
  "sample_id",
  "patient_id",
  "condition_main",
  "include_paired",
  "library_size"
)

paired_manifest <- readr::read_tsv(paired_manifest_path, show_col_types = FALSE)
counts <- readr::read_tsv(counts_path, show_col_types = FALSE)

missing_manifest_cols <- setdiff(required_manifest_cols, colnames(paired_manifest))
if (length(missing_manifest_cols) > 0) {
  stop(
    paste0(
      "Manifest missing required columns: ",
      paste(missing_manifest_cols, collapse = ", ")
    )
  )
}

if (any(duplicated(paired_manifest$sample_id))) {
  dup_ids <- unique(paired_manifest$sample_id[duplicated(paired_manifest$sample_id)])
  stop(paste0("Duplicate sample_id values in manifest: ", paste(dup_ids, collapse = ", ")))
}

paired_manifest <- paired_manifest %>%
  mutate(include_paired = toupper(include_paired)) %>%
  filter(include_paired == "TRUE") %>%
  mutate(
    patient_id = as.character(patient_id),
    condition_main = as.character(condition_main)
  )

if (nrow(paired_manifest) == 0) {
  stop("No samples with include_paired == TRUE in manifest.")
}

if (any(duplicated(paired_manifest$sample_id))) {
  dup_ids <- unique(paired_manifest$sample_id[duplicated(paired_manifest$sample_id)])
  stop(
    paste0(
      "Duplicate sample_id values after include_paired filter: ",
      paste(dup_ids, collapse = ", ")
    )
  )
}

invalid_conditions <- setdiff(unique(paired_manifest$condition_main), c("Tumor", "Normal"))
if (length(invalid_conditions) > 0) {
  stop(
    paste0(
      "include_paired samples contain invalid condition_main values: ",
      paste(invalid_conditions, collapse = ", ")
    )
  )
}

pair_check <- paired_manifest %>%
  group_by(patient_id) %>%
  summarize(
    n_samples = n(),
    n_tumor = sum(condition_main == "Tumor"),
    n_normal = sum(condition_main == "Normal"),
    .groups = "drop"
  )

bad_pairs <- pair_check %>%
  filter(!(n_samples == 2 & n_tumor == 1 & n_normal == 1))

if (nrow(bad_pairs) > 0) {
  bad_msg <- apply(as.data.frame(bad_pairs), 1, function(x) {
    paste0(
      x[["patient_id"]], "(n=", x[["n_samples"]],
      ", tumor=", x[["n_tumor"]], ", normal=", x[["n_normal"]], ")"
    )
  })
  stop(
    paste0(
      "Invalid paired structure in include_paired cohort: ",
      paste(bad_msg, collapse = "; ")
    )
  )
}

if (!"gene_id" %in% colnames(counts)) {
  stop("Counts file must include a 'gene_id' column.")
}

if (any(duplicated(counts$gene_id))) {
  dup_gene_ids <- unique(counts$gene_id[duplicated(counts$gene_id)])
  stop(
    paste0(
      "Counts file has duplicate gene_id values. Example IDs: ",
      paste(head(dup_gene_ids, 10), collapse = ", ")
    )
  )
}

count_samples <- setdiff(colnames(counts), "gene_id")
manifest_samples <- paired_manifest$sample_id
missing_in_counts <- setdiff(manifest_samples, count_samples)
if (length(missing_in_counts) > 0) {
  stop(
    paste0(
      "Manifest samples missing in counts matrix: ",
      paste(missing_in_counts, collapse = ", ")
    )
  )
}

# Keep exact manifest order for all downstream objects.
ordered_samples <- manifest_samples

counts_subset <- counts %>%
  select(gene_id, all_of(ordered_samples))

count_matrix_numeric <- as.matrix(counts_subset[, -1])
rownames(count_matrix_numeric) <- counts_subset$gene_id

if (!is.numeric(count_matrix_numeric)) {
  stop("Counts matrix contains non-numeric values.")
}
if (any(is.na(count_matrix_numeric))) {
  stop("Counts matrix contains NA values.")
}
if (any(count_matrix_numeric < 0)) {
  stop("Counts matrix contains negative values.")
}

max_abs_dev <- max(abs(count_matrix_numeric - round(count_matrix_numeric)))
if (!is.finite(max_abs_dev) || max_abs_dev > 1e-6) {
  stop(
    paste0(
      "counts.tsv is not raw counts",
      " (max abs deviation from integer = ",
      signif(max_abs_dev, 6),
      ")"
    )
  )
}

count_matrix <- round(count_matrix_numeric)
storage.mode(count_matrix) <- "integer"

#### 3) Construct DESeqDataSet + design ####
col_data <- paired_manifest %>%
  slice(match(ordered_samples, sample_id)) %>%
  mutate(
    patient_id = factor(patient_id),
    condition_main = factor(condition_main, levels = c("Normal", "Tumor"))
  ) %>%
  as.data.frame()

rownames(col_data) <- col_data$sample_id

stopifnot(identical(ordered_samples, colnames(count_matrix)))
stopifnot(identical(ordered_samples, col_data$sample_id))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ patient_id + condition_main
)

#### 4) Run DESeq() ####
# These are DESeq2 defaults, stated explicitly so the inferential contract is auditable.
dds <- DESeq(
  dds,
  test = "Wald",
  fitType = "parametric",
  sfType = "ratio",
  minReplicatesForReplace = 7
)
res <- results(
  dds,
  contrast = c("condition_main", "Tumor", "Normal"),
  # DESeq2's default alpha=0.1 is retained because it controls independent-filter
  # optimization. Reporting thresholds are applied separately at padj<0.05.
  alpha = 0.1,
  independentFiltering = TRUE,
  cooksCutoff = TRUE,
  pAdjustMethod = "BH"
)

#### 5) Extract results + shrinkage ####
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df %>%
  select(gene_id, everything()) %>%
  arrange(is.na(padj), padj)

# Deterministic coefficient selection for the paired Tumor-vs-Normal effect.
coef_name <- grep("^condition_main_", resultsNames(dds), value = TRUE)
if (length(coef_name) != 1) {
  stop(
    paste0(
      "Could not uniquely identify condition_main coefficient for lfcShrink. ",
      "Found: ", paste(coef_name, collapse = ", ")
    )
  )
}

# Deterministic LFC shrinkage method (fixed to type = "normal").
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "normal")
res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)[res_df$gene_id, "log2FoldChange"]
cat("lfcShrink applied with type='normal' using coef:", coef_name, "\n")

#### 6) Exploratory transformation + diagnostics ####
# VST is for PCA/heatmaps only. The DE model continues to use raw integer counts.
vst_obj <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_df <- as.data.frame(assay(vst_obj))
vst_df$gene_id <- rownames(vst_df)
vst_df <- vst_df %>% select(gene_id, all_of(ordered_samples))

all_zero <- res_df$baseMean == 0
pvalue_na <- is.na(res_df$pvalue)
padj_na <- is.na(res_df$padj)
independent_filtered <- !pvalue_na & padj_na
positive_mean_pvalue_na <- !all_zero & pvalue_na
filter_threshold <- metadata(res)$filterThreshold
if (is.null(filter_threshold) || length(filter_threshold) == 0) {
  filter_threshold <- NA_real_
}

diagnostics_df <- data.frame(
  metric = c(
    "genes_tested",
    "all_zero_genes",
    "genes_with_pvalue_na",
    "positive_mean_genes_with_pvalue_na",
    "genes_removed_by_independent_filtering",
    "genes_with_padj_non_na",
    "independent_filter_baseMean_threshold",
    "wald_test",
    "bh_adjustment",
    "cooks_cutoff_enabled",
    "lfc_shrinkage_method",
    "vst_blind"
  ),
  value = c(
    nrow(res_df),
    sum(all_zero, na.rm = TRUE),
    sum(pvalue_na),
    sum(positive_mean_pvalue_na, na.rm = TRUE),
    sum(independent_filtered),
    sum(!padj_na),
    as.character(signif(as.numeric(filter_threshold)[1], 8)),
    "TRUE",
    "Benjamini-Hochberg",
    "TRUE",
    "normal",
    "FALSE"
  ),
  interpretation = c(
    "Rows in the raw-count DESeq2 result table.",
    "Genes with zero counts across all included samples; no statistical test is possible.",
    "Genes without a nominal p-value.",
    "Potential Cook's-distance/outlier-suppressed genes; zero means none were observed.",
    "Genes with a p-value but NA adjusted p-value after independent filtering.",
    "Genes retained for multiple-testing adjustment.",
    "DESeq2-selected mean-expression threshold used for independent filtering.",
    "DESeq2 coefficient test used for the Tumor-vs-Normal contrast.",
    "Multiple-testing correction applied by results().",
    "DESeq2 Cook's-distance filtering was requested.",
    "Method used to stabilize noisy log2 fold-change estimates.",
    "Transformation uses the fitted design/dispersion trend for exploratory plots."
  ),
  stringsAsFactors = FALSE
)

size_factors_df <- data.frame(
  sample_id = colnames(dds),
  size_factor = as.numeric(sizeFactors(dds)),
  stringsAsFactors = FALSE
)

#### 7) Write outputs ####
dir.create(dirname(de_results_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ma_plot_path), recursive = TRUE, showWarnings = FALSE)

# Main DE table for enrichment + figure scripts.
readr::write_tsv(res_df, de_results_path)
readr::write_tsv(vst_df, vst_path)
readr::write_tsv(diagnostics_df, diagnostics_path)
readr::write_tsv(size_factors_df, size_factors_path)

# Manifest subset actually used by DE; used by run_v2.sh strict cohort validation.
readr::write_tsv(
  paired_manifest %>% arrange(patient_id, condition_main),
  samples_used_path
)

# Legacy diagnostic MA plot (kept for continuity, not part of final F01-F07 set).
png(ma_plot_path, width = 2000, height = 1600, res = 300)
plotMA(res, main = "Paired Tumor vs Normal (v2)", alpha = 0.05)
dev.off()

# Session metadata for reproducibility auditing and output manifests.
capture.output(sessionInfo(), file = session_info_path)

sig_mask <- !is.na(res_df$padj) & res_df$padj < 0.05
top10 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  select(gene_id, log2FoldChange, padj) %>%
  slice_head(n = 10)

cat(
  "counts integrity check passed (max abs deviation from integer):",
  signif(max_abs_dev, 6),
  "\n"
)
cat("number of paired patients used:", n_distinct(paired_manifest$patient_id), "\n")
cat("number of samples used:", nrow(paired_manifest), "\n")
cat("number of genes tested:", nrow(res_df), "\n")
cat("number significant at padj < 0.05:", sum(sig_mask), "\n")
cat("top 10 genes by padj (gene_id, log2FoldChange, padj):\n")
print(top10)

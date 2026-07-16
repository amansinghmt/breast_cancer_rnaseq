#!/usr/bin/env Rscript

# Purpose: derive transparent presentation and sensitivity summaries from the
# canonical paired analysis. The canonical DE result is never replaced here.

required_packages <- c(
  "DESeq2", "readr", "dplyr", "tidyr", "AnnotationDbi", "org.Hs.eg.db"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  # Attach dplyr last so its data-frame verbs are not masked by AnnotationDbi generics.
  library(dplyr)
})

seed_value <- suppressWarnings(as.integer(Sys.getenv("PIPELINE_SEED", "20260227")))
if (is.na(seed_value)) stop("PIPELINE_SEED must be an integer.")
set.seed(seed_value)

de_path <- "results_v2/deseq2/deseq2_paired_v2_results.tsv"
samples_path <- "results_v2/deseq2/deseq2_paired_v2_samples_used.tsv"
vst_path <- "results_v2/deseq2/deseq2_paired_v2_vst.tsv"
counts_path <- "data/processed/counts.tsv"
manifest_path <- "results_v2/metadata/paired_manifest.tsv"
hallmark_path <- "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
go_path <- "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
go_representative_path <- "results_v2/enrichment/go_bp_ora_representative_v2.tsv"
go_tumor_path <- "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv"
go_normal_path <- "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv"
enrichment_diagnostics_path <- "results_v2/enrichment/enrichment_diagnostics_v2.tsv"
out_dir <- "results_v2/robustness"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

de <- readr::read_tsv(de_path, show_col_types = FALSE)
samples <- readr::read_tsv(samples_path, show_col_types = FALSE)
vst <- readr::read_tsv(vst_path, show_col_types = FALSE)
counts <- readr::read_tsv(counts_path, show_col_types = FALSE)
manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
hallmark <- readr::read_tsv(hallmark_path, show_col_types = FALSE)
go <- readr::read_tsv(go_path, show_col_types = FALSE)
go_representative <- readr::read_tsv(go_representative_path, show_col_types = FALSE)
go_tumor <- readr::read_tsv(go_tumor_path, show_col_types = FALSE)
go_normal <- readr::read_tsv(go_normal_path, show_col_types = FALSE)
enrichment_diagnostics <- readr::read_tsv(enrichment_diagnostics_path, show_col_types = FALSE)

required_de <- c(
  "gene_id", "baseMean", "log2FoldChange", "pvalue", "padj",
  "log2FoldChange_shrunk"
)
missing_de <- setdiff(required_de, colnames(de))
if (length(missing_de) > 0) {
  stop("DE table missing columns: ", paste(missing_de, collapse = ", "))
}

threshold_grid <- tidyr::crossing(
  padj_cutoff = c(0.05, 0.01),
  abs_shrunken_log2fc_cutoff = c(0, 1, 1.5)
)

threshold_summary <- threshold_grid %>%
  rowwise() %>%
  mutate(
    tumor_higher = sum(
      !is.na(de$padj) & de$padj < padj_cutoff &
        !is.na(de$log2FoldChange_shrunk) &
        de$log2FoldChange_shrunk >= abs_shrunken_log2fc_cutoff
    ),
    normal_higher = sum(
      !is.na(de$padj) & de$padj < padj_cutoff &
        !is.na(de$log2FoldChange_shrunk) &
        de$log2FoldChange_shrunk <= -abs_shrunken_log2fc_cutoff
    ),
    total = tumor_higher + normal_higher
  ) %>%
  ungroup()

readr::write_tsv(
  threshold_summary,
  file.path(out_dir, "de_threshold_sensitivity_v2.tsv")
)

clean_gene_id <- function(x) sub("\\..*$", "", as.character(x))
de <- de %>% mutate(gene_id_clean = clean_gene_id(gene_id))

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de$gene_id_clean),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
  arrange(ENSEMBL, SYMBOL) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

top_de <- de %>%
  left_join(symbol_map, by = c("gene_id_clean" = "ENSEMBL")) %>%
  filter(!is.na(padj)) %>%
  arrange(padj, desc(abs(log2FoldChange_shrunk))) %>%
  mutate(
    rank_by_padj = row_number(),
    direction = case_when(
      log2FoldChange_shrunk > 0 ~ "Higher in Tumor",
      log2FoldChange_shrunk < 0 ~ "Higher in Normal",
      TRUE ~ "No direction"
    ),
    passes_primary_reporting_rule =
      padj < 0.05 & abs(log2FoldChange_shrunk) >= 1
  ) %>%
  select(
    rank_by_padj, gene_id, SYMBOL, baseMean, log2FoldChange,
    log2FoldChange_shrunk, lfcSE, stat, pvalue, padj, direction,
    passes_primary_reporting_rule
  ) %>%
  slice_head(n = 500)

readr::write_tsv(top_de, file.path(out_dir, "top_de_genes_v2.tsv"))

complete_lfc <- complete.cases(de$log2FoldChange, de$log2FoldChange_shrunk)
lfc_spearman <- cor(
  de$log2FoldChange[complete_lfc],
  de$log2FoldChange_shrunk[complete_lfc],
  method = "spearman"
)

lfc_agreement <- data.frame(
  metric = "spearman_unshrunk_vs_shrunk_lfc",
  comparison_size = sum(complete_lfc),
  value = lfc_spearman,
  stringsAsFactors = FALSE
)

for (n in c(100, 500, 1000)) {
  unshrunk_top <- head(
    de$gene_id[order(abs(de$log2FoldChange), decreasing = TRUE, na.last = NA)],
    n
  )
  shrunk_top <- head(
    de$gene_id[order(abs(de$log2FoldChange_shrunk), decreasing = TRUE, na.last = NA)],
    n
  )
  lfc_agreement <- bind_rows(
    lfc_agreement,
    data.frame(
      metric = paste0("top_", n, "_absolute_lfc_overlap"),
      comparison_size = n,
      value = length(intersect(unshrunk_top, shrunk_top)),
      stringsAsFactors = FALSE
    )
  )
}

readr::write_tsv(lfc_agreement, file.path(out_dir, "lfc_agreement_v2.tsv"))

sample_order <- samples$sample_id
count_matrix <- as.matrix(counts[, sample_order, drop = FALSE])
rownames(count_matrix) <- counts$gene_id
storage.mode(count_matrix) <- "integer"

# Predefined sensitivity rule: retain genes with count >=10 in at least two samples.
# This is a robustness check only; it does not replace the canonical all-gene model.
prefilter_keep <- rowSums(count_matrix >= 10) >= 2
col_data <- samples %>%
  mutate(
    patient_id = factor(as.character(patient_id)),
    condition_main = factor(as.character(condition_main), levels = c("Normal", "Tumor"))
  ) %>%
  as.data.frame()
rownames(col_data) <- col_data$sample_id

prefilter_dds <- DESeqDataSetFromMatrix(
  countData = count_matrix[prefilter_keep, , drop = FALSE],
  colData = col_data,
  design = ~ patient_id + condition_main
)
prefilter_dds <- DESeq(
  prefilter_dds,
  test = "Wald",
  fitType = "parametric",
  sfType = "ratio",
  minReplicatesForReplace = 7,
  quiet = TRUE
)
prefilter_res <- results(
  prefilter_dds,
  contrast = c("condition_main", "Tumor", "Normal"),
  alpha = 0.1,
  independentFiltering = TRUE,
  cooksCutoff = TRUE,
  pAdjustMethod = "BH"
)
prefilter_coef <- grep(
  "^condition_main_",
  resultsNames(prefilter_dds),
  value = TRUE
)
if (length(prefilter_coef) != 1) stop("Could not identify prefilter condition coefficient.")
prefilter_shrunk <- lfcShrink(
  prefilter_dds,
  coef = prefilter_coef,
  type = "normal",
  quiet = TRUE
)

prefilter_df <- data.frame(
  gene_id = rownames(prefilter_res),
  padj_prefilter = prefilter_res$padj,
  log2fc_shrunk_prefilter = prefilter_shrunk$log2FoldChange,
  stringsAsFactors = FALSE
)

prefilter_comparison <- de %>%
  select(gene_id, padj, log2FoldChange_shrunk) %>%
  inner_join(prefilter_df, by = "gene_id")

canonical_sig <- !is.na(prefilter_comparison$padj) & prefilter_comparison$padj < 0.05
prefilter_sig <- !is.na(prefilter_comparison$padj_prefilter) &
  prefilter_comparison$padj_prefilter < 0.05
prefilter_lfc_complete <- complete.cases(
  prefilter_comparison$log2FoldChange_shrunk,
  prefilter_comparison$log2fc_shrunk_prefilter
)

prefilter_summary <- data.frame(
  metric = c(
    "canonical_genes",
    "prefilter_genes_retained",
    "canonical_padj_lt_0.05_among_retained",
    "prefilter_padj_lt_0.05",
    "significant_in_both",
    "canonical_only_significant",
    "prefilter_only_significant",
    "canonical_padj_lt_0.05_abs_shrunk_lfc_ge_1_among_retained",
    "prefilter_padj_lt_0.05_abs_shrunk_lfc_ge_1",
    "shrunk_lfc_spearman"
  ),
  value = c(
    nrow(de),
    sum(prefilter_keep),
    sum(canonical_sig),
    sum(prefilter_sig),
    sum(canonical_sig & prefilter_sig),
    sum(canonical_sig & !prefilter_sig),
    sum(!canonical_sig & prefilter_sig),
    sum(
      canonical_sig &
        abs(prefilter_comparison$log2FoldChange_shrunk) >= 1,
      na.rm = TRUE
    ),
    sum(
      prefilter_sig &
        abs(prefilter_comparison$log2fc_shrunk_prefilter) >= 1,
      na.rm = TRUE
    ),
    cor(
      prefilter_comparison$log2FoldChange_shrunk[prefilter_lfc_complete],
      prefilter_comparison$log2fc_shrunk_prefilter[prefilter_lfc_complete],
      method = "spearman"
    )
  ),
  stringsAsFactors = FALSE
)

readr::write_tsv(
  prefilter_summary,
  file.path(out_dir, "low_count_prefilter_sensitivity_v2.tsv")
)

vst_matrix <- as.matrix(vst[, sample_order, drop = FALSE])
rownames(vst_matrix) <- vst$gene_id
storage.mode(vst_matrix) <- "numeric"
gene_variance <- apply(vst_matrix, 1, var)
top_index <- order(gene_variance, decreasing = TRUE, na.last = NA)[
  seq_len(min(2000, sum(is.finite(gene_variance) & gene_variance > 0)))
]
pca <- prcomp(t(vst_matrix[top_index, , drop = FALSE]), center = TRUE, scale. = FALSE)
scores <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(
    samples %>% select(sample_id, patient_id, condition_main),
    by = "sample_id"
  )

robust_scale <- function(x) {
  value <- mad(x, constant = 1.4826)
  if (!is.finite(value) || value == 0) sd(x) else value
}

pc1_scale <- robust_scale(scores$PC1)
pc2_scale <- robust_scale(scores$PC2)
scores <- scores %>%
  mutate(
    robust_pc_distance = sqrt(
      ((PC1 - median(PC1)) / pc1_scale)^2 +
        ((PC2 - median(PC2)) / pc2_scale)^2
    ),
    exploratory_outlier_flag = robust_pc_distance > 3.5
  ) %>%
  arrange(desc(robust_pc_distance))

readr::write_tsv(scores, file.path(out_dir, "pca_outlier_summary_v2.tsv"))

manifest_summary <- manifest %>%
  mutate(
    include_paired = toupper(as.character(include_paired)),
    cohort_status = ifelse(include_paired == "TRUE", "Included", "Excluded"),
    reason = ifelse(
      cohort_status == "Included",
      "included_paired_cohort",
      ifelse(is.na(exclude_reason) | exclude_reason == "", "unspecified", exclude_reason)
    )
  ) %>%
  count(cohort_status, reason, name = "samples") %>%
  arrange(cohort_status, desc(samples), reason)

readr::write_tsv(
  manifest_summary,
  file.path(out_dir, "cohort_inclusion_summary_v2.tsv")
)

metrics <- data.frame(
  metric = c(
    "manifest_rows", "included_samples", "paired_patients", "tumor_samples",
    "normal_samples", "genes_tested", "genes_padj_lt_0.05",
    "genes_padj_lt_0.05_abs_shrunk_lfc_ge_1", "hallmark_sets_tested",
    "hallmark_sets_padj_lt_0.05", "go_terms_tested",
    "go_terms_padj_lt_0.05", "go_representative_terms",
    "go_tumor_higher_terms_tested", "go_tumor_higher_terms_padj_lt_0.05",
    "go_normal_higher_terms_tested", "go_normal_higher_terms_padj_lt_0.05",
    "pca_explained_pc1_percent", "pca_explained_pc2_percent",
    "pca_exploratory_outliers"
  ),
  value = c(
    nrow(manifest),
    sum(toupper(as.character(manifest$include_paired)) == "TRUE"),
    n_distinct(samples$patient_id),
    sum(samples$condition_main == "Tumor"),
    sum(samples$condition_main == "Normal"),
    nrow(de),
    sum(de$padj < 0.05, na.rm = TRUE),
    sum(de$padj < 0.05 & abs(de$log2FoldChange_shrunk) >= 1, na.rm = TRUE),
    nrow(hallmark),
    sum(hallmark$padj < 0.05, na.rm = TRUE),
    nrow(go),
    sum(go$p.adjust < 0.05, na.rm = TRUE),
    nrow(go_representative),
    nrow(go_tumor),
    sum(go_tumor$p.adjust < 0.05, na.rm = TRUE),
    nrow(go_normal),
    sum(go_normal$p.adjust < 0.05, na.rm = TRUE),
    (pca$sdev[1]^2 / sum(pca$sdev^2)) * 100,
    (pca$sdev[2]^2 / sum(pca$sdev^2)) * 100,
    sum(scores$exploratory_outlier_flag)
  ),
  stringsAsFactors = FALSE
)

readr::write_tsv(metrics, file.path(out_dir, "analysis_metrics_v2.tsv"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo_robustness_v2.txt"))

cat("Robustness summaries written to", out_dir, "\n")
cat("Prefilter retained genes:", sum(prefilter_keep), "\n")
cat("Canonical/prefilter significant overlap:", sum(canonical_sig & prefilter_sig), "\n")
cat("PCA exploratory outlier flags:", sum(scores$exploratory_outlier_flag), "\n")

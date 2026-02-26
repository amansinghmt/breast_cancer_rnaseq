#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
})

manifest_path <- "data/metadata/sample_manifest.tsv"
counts_path <- "data/processed/counts.tsv"

results_path <- "results_v2/differential_expression/deseq2_paired_v2_results.tsv"
samples_used_path <- "results_v2/differential_expression/deseq2_paired_v2_samples_used.tsv"
session_info_path <- "results_v2/differential_expression/sessionInfo_paired_v2.txt"
ma_plot_path <- "figures_v2/de/ma_plot_paired_v2.png"

required_manifest_cols <- c(
  "sample_id",
  "patient_id",
  "condition_main",
  "include_paired",
  "library_size"
)

manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))
if (length(missing_manifest_cols) > 0) {
  stop(
    paste0(
      "Manifest missing required columns: ",
      paste(missing_manifest_cols, collapse = ", ")
    )
  )
}

if (any(duplicated(manifest$sample_id))) {
  dup_ids <- unique(manifest$sample_id[duplicated(manifest$sample_id)])
  stop(paste0("Duplicate sample_id values in manifest: ", paste(dup_ids, collapse = ", ")))
}

manifest_paired <- manifest %>%
  mutate(include_paired = toupper(include_paired)) %>%
  filter(include_paired == "TRUE") %>%
  mutate(
    patient_id = as.character(patient_id),
    condition_main = as.character(condition_main)
  )

if (nrow(manifest_paired) == 0) {
  stop("No samples with include_paired == TRUE in manifest.")
}

if (any(duplicated(manifest_paired$sample_id))) {
  dup_ids <- unique(manifest_paired$sample_id[duplicated(manifest_paired$sample_id)])
  stop(
    paste0(
      "Duplicate sample_id values after include_paired filter: ",
      paste(dup_ids, collapse = ", ")
    )
  )
}

invalid_conditions <- setdiff(unique(manifest_paired$condition_main), c("Tumor", "Normal"))
if (length(invalid_conditions) > 0) {
  stop(
    paste0(
      "include_paired samples contain invalid condition_main values: ",
      paste(invalid_conditions, collapse = ", ")
    )
  )
}

pair_check <- manifest_paired %>%
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

counts <- readr::read_tsv(counts_path, show_col_types = FALSE)
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
manifest_samples <- manifest_paired$sample_id
missing_in_counts <- setdiff(manifest_samples, count_samples)
if (length(missing_in_counts) > 0) {
  stop(
    paste0(
      "Manifest samples missing in counts matrix: ",
      paste(missing_in_counts, collapse = ", ")
    )
  )
}

# Keep the exact manifest order for all downstream objects.
ordered_samples <- manifest_samples

counts_subset <- counts %>%
  select(gene_id, all_of(ordered_samples))

count_matrix <- as.matrix(counts_subset[, -1])
rownames(count_matrix) <- counts_subset$gene_id
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

col_data <- manifest_paired %>%
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

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition_main", "Tumor", "Normal"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df %>%
  select(gene_id, everything()) %>%
  arrange(is.na(padj), padj)

# Optional shrinkage with apeglm.
if (requireNamespace("apeglm", quietly = TRUE)) {
  coef_name <- grep("^condition_main_", resultsNames(dds), value = TRUE)
  if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)[res_df$gene_id, "log2FoldChange"]
    cat("lfcShrink applied with apeglm using coef:", coef_name, "\n")
  } else {
    res_df$log2FoldChange_shrunk <- NA_real_
    cat("Skipping lfcShrink: could not uniquely identify condition_main coefficient.\n")
  }
} else {
  res_df$log2FoldChange_shrunk <- NA_real_
  cat("Skipping lfcShrink: package 'apeglm' is not available.\n")
}

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ma_plot_path), recursive = TRUE, showWarnings = FALSE)

readr::write_tsv(res_df, results_path)
readr::write_tsv(
  manifest_paired %>% arrange(patient_id, condition_main),
  samples_used_path
)

png(ma_plot_path, width = 2000, height = 1600, res = 300)
plotMA(res, main = "Paired Tumor vs Normal (v2)", alpha = 0.05)
dev.off()

capture.output(sessionInfo(), file = session_info_path)

sig_mask <- !is.na(res_df$padj) & res_df$padj < 0.05
top10 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  select(gene_id, log2FoldChange, padj) %>%
  slice_head(n = 10)

cat("number of paired patients used:", n_distinct(manifest_paired$patient_id), "\n")
cat("number of samples used:", nrow(manifest_paired), "\n")
cat("number of genes tested:", nrow(res_df), "\n")
cat("number significant at padj < 0.05:", sum(sig_mask), "\n")
cat("top 10 genes by padj (gene_id, log2FoldChange, padj):\n")
print(top10)

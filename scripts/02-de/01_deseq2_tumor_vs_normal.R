#!/usr/bin/env Rscript
# LEGACY SCRIPT: retained for reference; not used by scripts/run_v2.sh.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

counts_path <- "data/processed/counts.tsv"
metadata_path <- "data/metadata/metadata.tsv"
results_path <- "results/differential_expression/deseq2_results.tsv"
ma_plot_path <- "figures/de/ma_plot.png"

required_columns <- c(
  "sample_id",
  "condition",
  "lfs_status",
  "use_main_tumor_vs_normal"
)

metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE)
missing_cols <- setdiff(required_columns, colnames(metadata))
if (length(missing_cols) > 0) {
  stop(paste0("Metadata missing required columns: ", paste(missing_cols, collapse = ", ")))
}

metadata_filtered <- metadata %>%
  filter(use_main_tumor_vs_normal == "yes") %>%
  filter(condition %in% c("Tumor", "Normal"))

if (any(duplicated(metadata_filtered$sample_id))) {
  dup_ids <- unique(metadata_filtered$sample_id[duplicated(metadata_filtered$sample_id)])
  stop(paste0("Duplicate sample_id values in metadata: ", paste(dup_ids, collapse = ", ")))
}

tumor_count <- sum(metadata_filtered$condition == "Tumor")
normal_count <- sum(metadata_filtered$condition == "Normal")
if (tumor_count < 10 || normal_count < 5) {
  stop(paste0(
    "Insufficient samples after filtering (Tumor: ", tumor_count,
    ", Normal: ", normal_count, ")"
  ))
}

metadata_filtered <- metadata_filtered %>%
  mutate(
    condition = factor(condition, levels = c("Normal", "Tumor")),
    lfs_status = factor(lfs_status)
  )

counts <- readr::read_tsv(counts_path, show_col_types = FALSE)
if (!"gene_id" %in% colnames(counts)) {
  stop("Counts file must include a 'gene_id' column.")
}

count_samples <- setdiff(colnames(counts), "gene_id")
metadata_samples <- metadata_filtered$sample_id

missing_in_counts <- setdiff(metadata_samples, count_samples)
if (length(missing_in_counts) > 0) {
  stop(paste0(
    "Missing samples in counts: ",
    paste(missing_in_counts, collapse = ", ")
  ))
}

common_samples <- intersect(metadata_samples, count_samples)
counts_subset <- counts %>%
  select(gene_id, all_of(common_samples))

gene_ids <- counts_subset$gene_id
count_matrix <- as.matrix(counts_subset[, -1])
rownames(count_matrix) <- gene_ids
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

metadata_filtered <- metadata_filtered %>%
  filter(sample_id %in% common_samples) %>%
  slice(match(common_samples, sample_id))
col_data <- as.data.frame(metadata_filtered)
rownames(col_data) <- col_data$sample_id

stopifnot(identical(common_samples, colnames(count_matrix)))
stopifnot(identical(common_samples, col_data$sample_id))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ lfs_status + condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df %>%
  select(gene_id, everything()) %>%
  arrange(is.na(padj), padj)

dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(res_df, results_path)

dir.create(dirname(ma_plot_path), recursive = TRUE, showWarnings = FALSE)
png(ma_plot_path, width = 2000, height = 1600, res = 300)
plotMA(res, main = "Tumor vs Normal", alpha = 0.05)
dev.off()

sig_mask <- !is.na(res_df$padj) & res_df$padj < 0.05
up_mask <- sig_mask & res_df$log2FoldChange > 0
down_mask <- sig_mask & res_df$log2FoldChange < 0

cat("number of genes tested:", nrow(res_df), "\n")
cat("number of significant genes (padj < 0.05):", sum(sig_mask), "\n")
cat("number upregulated in Tumor:", sum(up_mask), "\n")
cat("number downregulated in Tumor:", sum(down_mask), "\n")

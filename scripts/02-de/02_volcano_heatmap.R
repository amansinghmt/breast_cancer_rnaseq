#!/usr/bin/env Rscript
# LEGACY SCRIPT: retained for reference; not used by scripts/run_v2.sh.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

results_path <- "results/differential_expression/deseq2_results.tsv"
counts_path <- "data/processed/counts.tsv"
metadata_path <- "data/metadata/metadata.tsv"

volcano_path <- "figures/de/volcano.png"
heatmap_path <- "figures/de/heatmap_top50.png"
top_genes_path <- "results/differential_expression/top_genes.tsv"

required_metadata_columns <- c(
  "sample_id",
  "condition",
  "tp53_status",
  "use_main_tumor_vs_normal"
)

results <- readr::read_tsv(results_path, show_col_types = FALSE)
if (!"padj" %in% colnames(results)) {
  stop("DESeq2 results must include a 'padj' column.")
}
if (!"gene_id" %in% colnames(results)) {
  stop("DESeq2 results must include a 'gene_id' column.")
}
if (!"log2FoldChange" %in% colnames(results)) {
  stop("DESeq2 results must include a 'log2FoldChange' column.")
}

results <- results %>%
  filter(!is.na(padj))

sig_mask <- results$padj < 0.05
n_sig <- sum(sig_mask)

results <- results %>%
  mutate(
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    is_significant = ifelse(padj < 0.05, "significant", "not_significant")
  )

volcano <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = is_significant), alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666") +
  scale_color_manual(values = c("significant" = "#D55E00", "not_significant" = "#0072B2")) +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = "Significance") +
  theme_minimal()

if (requireNamespace("ggrepel", quietly = TRUE)) {
  top10 <- results %>%
    arrange(padj) %>%
    slice_head(n = 10)
  volcano <- volcano +
    ggrepel::geom_text_repel(
      data = top10,
      aes(label = gene_id),
      size = 3,
      max.overlaps = 50
    )
}

dir.create(dirname(volcano_path), recursive = TRUE, showWarnings = FALSE)
ggsave(volcano_path, plot = volcano, width = 7, height = 5, dpi = 300)

sig_results <- results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

top_genes <- sig_results %>%
  select(gene_id, log2FoldChange, padj) %>%
  slice_head(n = 50)

dir.create(dirname(top_genes_path), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(top_genes, top_genes_path)

if (nrow(top_genes) == 0) {
  stop("No significant genes found (padj < 0.05); heatmap cannot be generated.")
}

metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE)
missing_cols <- setdiff(required_metadata_columns, colnames(metadata))
if (length(missing_cols) > 0) {
  stop(paste0("Metadata missing required columns: ", paste(missing_cols, collapse = ", ")))
}

metadata_filtered <- metadata %>%
  filter(use_main_tumor_vs_normal == "yes") %>%
  filter(condition %in% c("Tumor", "Normal"))

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

count_matrix <- as.matrix(counts_subset[, -1])
rownames(count_matrix) <- counts_subset$gene_id

metadata_filtered <- metadata_filtered %>%
  filter(sample_id %in% common_samples) %>%
  slice(match(common_samples, sample_id))

stopifnot(identical(common_samples, colnames(count_matrix)))
stopifnot(identical(common_samples, metadata_filtered$sample_id))

library_sizes <- colSums(count_matrix)
if (any(library_sizes == 0)) {
  stop("One or more samples have zero library size.")
}

cpm <- sweep(count_matrix, 2, library_sizes, "/") * 1e6
log2_cpm <- log2(cpm + 1)

heatmap_genes <- top_genes$gene_id
missing_genes <- setdiff(heatmap_genes, rownames(log2_cpm))
if (length(missing_genes) > 0) {
  stop(paste0(
    "Missing genes in counts for heatmap: ",
    paste(missing_genes, collapse = ", ")
  ))
}

heatmap_matrix <- log2_cpm[heatmap_genes, , drop = FALSE]
heatmap_scaled <- t(scale(t(heatmap_matrix)))
heatmap_scaled[is.na(heatmap_scaled)] <- 0

annotation_col <- metadata_filtered %>%
  select(sample_id, condition, tp53_status) %>%
  as.data.frame()
rownames(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL


dir.create(dirname(heatmap_path), recursive = TRUE, showWarnings = FALSE)
pheatmap::pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  fontsize_row = 6,
  fontsize_col = 6,
  filename = heatmap_path
)

cat("number of significant genes:", n_sig, "\n")
cat("created:", volcano_path, "\n")
cat("created:", heatmap_path, "\n")
cat("created:", top_genes_path, "\n")

#!/usr/bin/env Rscript

required_cran <- c("readr", "dplyr", "ggplot2", "msigdbr")
required_bioc <- c("fgsea", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0 || length(missing_bioc) > 0) {
  install_lines <- character(0)
  if (length(missing_cran) > 0) {
    install_lines <- c(
      install_lines,
      paste0(
        "install.packages(c(",
        paste(sprintf("'%s'", missing_cran), collapse = ", "),
        "))"
      )
    )
  }
  if (length(missing_bioc) > 0) {
    install_lines <- c(
      install_lines,
      "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')",
      paste0(
        "BiocManager::install(c(",
        paste(sprintf("'%s'", missing_bioc), collapse = ", "),
        "))"
      )
    )
  }
  stop(
    paste(
      "Missing packages:",
      paste(c(missing_cran, missing_bioc), collapse = ", "),
      "\nInstall with:\n",
      paste(install_lines, collapse = "\n"),
      sep = " "
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

results_path <- "results/differential_expression/deseq2_results.tsv"

hallmark_out <- "results/enrichment/hallmark_gsea.tsv"
hallmark_plot <- "figures/enrichment/hallmark_gsea_top10.png"

go_out <- "results/enrichment/go_enrich.tsv"
go_plot <- "figures/enrichment/go_barplot_top15.png"

results <- readr::read_tsv(results_path, show_col_types = FALSE)
required_cols <- c("gene_id", "log2FoldChange", "pvalue", "padj")
missing_cols <- setdiff(required_cols, colnames(results))
if (length(missing_cols) > 0) {
  stop(paste0("DESeq2 results missing columns: ", paste(missing_cols, collapse = ", ")))
}

results_clean <- results %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(pvalue)) %>%
  arrange(pvalue) %>%
  distinct(gene_id, .keep_all = TRUE)

n_genes <- nrow(results_clean)

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(results_clean$gene_id),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, SYMBOL)

results_symbol <- results_clean %>%
  left_join(symbol_map, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL)) %>%
  arrange(pvalue) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  mutate(rank_metric = sign(log2FoldChange) * -log10(pmax(pvalue, .Machine$double.xmin)))

n_symbol <- nrow(results_symbol)

hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

hallmark_pathways <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

hallmark_tbl <- data.frame()
if (n_symbol > 0) {
  stats <- results_symbol$rank_metric
  names(stats) <- results_symbol$SYMBOL
  stats <- sort(stats, decreasing = TRUE)

  hallmark_res <- fgsea::fgsea(
    pathways = hallmark_pathways,
    stats = stats,
    minSize = 15,
    maxSize = 500,
    nperm = 10000
  )

  hallmark_tbl <- as.data.frame(hallmark_res) %>%
    arrange(padj)
}

dir.create(dirname(hallmark_out), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(hallmark_tbl, hallmark_out)

hallmark_sig_count <- 0
if (nrow(hallmark_tbl) > 0) {
  hallmark_sig_count <- sum(hallmark_tbl$padj < 0.05, na.rm = TRUE)
}
if (hallmark_sig_count == 0) {
  warning("No significant Hallmark pathways (padj < 0.05).")
}

if (nrow(hallmark_tbl) > 0) {
  top10 <- hallmark_tbl %>%
    arrange(padj) %>%
    slice_head(n = 10)
  padj_vals <- top10$padj
  if (all(is.na(padj_vals))) {
    subtitle_text <- "padj: NA"
  } else {
    subtitle_text <- sprintf(
      "padj range: %.3g - %.3g",
      min(padj_vals, na.rm = TRUE),
      max(padj_vals, na.rm = TRUE)
    )
  }
  hallmark_plot_obj <- ggplot(top10, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(fill = "#0072B2") +
    coord_flip() +
    labs(
      title = "Hallmark GSEA (Top 10)",
      subtitle = subtitle_text,
      x = NULL,
      y = "Normalized enrichment score"
    ) +
    theme_minimal()
} else {
  hallmark_plot_obj <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No Hallmark results", size = 4) +
    theme_void()
}

dir.create(dirname(hallmark_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(hallmark_plot, plot = hallmark_plot_obj, width = 7, height = 5, dpi = 300)

sig_genes <- results_clean %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1)

entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sig_genes$gene_id),
  keytype = "ENSEMBL",
  columns = "ENTREZID"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENSEMBL, ENTREZID)

entrez_ids <- unique(entrez_map$ENTREZID)

go_df <- data.frame()
if (length(entrez_ids) > 0) {
  go_res <- clusterProfiler::enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  go_df <- as.data.frame(go_res)
}

dir.create(dirname(go_out), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(go_df, go_out)

go_sig_count <- 0
if (nrow(go_df) > 0) {
  go_sig_count <- sum(go_df$p.adjust < 0.05, na.rm = TRUE)
}
if (go_sig_count == 0) {
  warning("No significant GO BP terms (p.adjust < 0.05).")
}

if (nrow(go_df) > 0) {
  top_go <- go_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    mutate(neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin)))
  top_go$Description <- factor(top_go$Description, levels = rev(top_go$Description))

  go_plot_obj <- ggplot(top_go, aes(x = Description, y = neg_log10_padj)) +
    geom_col(fill = "#D55E00") +
    coord_flip() +
    labs(
      title = "GO BP enrichment (Top 15)",
      x = NULL,
      y = "-log10 adjusted p-value"
    ) +
    theme_minimal()
} else {
  go_plot_obj <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No GO results", size = 4) +
    theme_void()
}

dir.create(dirname(go_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(go_plot, plot = go_plot_obj, width = 7, height = 5, dpi = 300)

cat("number of genes in DE table:", n_genes, "\n")
cat("number of genes with SYMBOL for GSEA:", n_symbol, "\n")
cat("number of significant genes for GO:", length(entrez_ids), "\n")
cat("number of significant Hallmark pathways:", hallmark_sig_count, "\n")
cat("number of significant GO BP terms:", go_sig_count, "\n")

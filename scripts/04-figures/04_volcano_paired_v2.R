#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall with:\ninstall.packages(c(",
      paste(sprintf("'%s'", missing_packages), collapse = ", "),
      "))"
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

de_path <- "results_v2/deseq2/deseq2_paired_v2_results.tsv"
figure_dir <- "figures_v2/final"
vector_dir <- "figures_v2/vector"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F04_de_volcano_paired_v2.png")
output_pdf_path <- file.path(vector_dir, "F04_de_volcano_paired_v2.pdf")
fig_manifest_path <- file.path(results_dir, "fig_manifest.tsv")

assert_columns <- function(df, cols, label) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        label,
        " missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

clean_gene_id <- function(x) {
  sub("\\..*$", "", as.character(x))
}

de <- readr::read_tsv(de_path, show_col_types = FALSE)
assert_columns(de, c("gene_id", "padj", "log2FoldChange"), "deseq2_paired_v2_results.tsv")

de <- de %>%
  mutate(gene_id_clean = clean_gene_id(gene_id))

lfc_col <- "log2FoldChange"
if ("log2FoldChange_shrunk" %in% colnames(de) && !all(is.na(de$log2FoldChange_shrunk))) {
  lfc_col <- "log2FoldChange_shrunk"
}

label_col_priority <- c("gene_symbol", "gene_name", "symbol", "external_gene_name", "hgnc_symbol")
label_col_hits <- label_col_priority[label_col_priority %in% colnames(de)]
label_col <- if (length(label_col_hits) > 0) label_col_hits[[1]] else NA_character_

if (!is.na(label_col)) {
  de <- de %>%
    mutate(
      SYMBOL = trimws(as.character(.data[[label_col]])),
      SYMBOL = ifelse(SYMBOL == "", NA_character_, SYMBOL)
    )
} else if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
           requireNamespace("AnnotationDbi", quietly = TRUE)) {
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(de$gene_id_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
  map_df <- map_df %>%
    dplyr::select(ENSEMBL, SYMBOL) %>%
    filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  de <- de %>%
    left_join(map_df, by = c("gene_id_clean" = "ENSEMBL"))
} else {
  de <- de %>%
    mutate(SYMBOL = NA_character_)
}

volcano_df <- de %>%
  mutate(
    log2FC = .data[[lfc_col]],
    SYMBOL = trimws(as.character(SYMBOL))
  ) %>%
  filter(!is.na(log2FC), !is.na(padj)) %>%
  mutate(
    SYMBOL = ifelse(SYMBOL == "", NA_character_, SYMBOL),
    label = ifelse(!is.na(SYMBOL), SYMBOL, gene_id_clean),
    padj_floor = pmax(padj, .Machine$double.xmin),
    neg_log10_padj = -log10(padj_floor),
    lfc_plot = pmax(pmin(log2FC, 8), -8),
    de_status = case_when(
      padj < 0.05 & log2FC >= 1 ~ "Tumor-higher",
      padj < 0.05 & log2FC <= -1 ~ "Normal-higher",
      padj < 0.05 ~ "FDR only",
      TRUE ~ "Not significant"
    )
  )

if (nrow(volcano_df) == 0) {
  stop("No rows available for volcano plot after filtering non-missing padj/log2FC.")
}

sig_df <- volcano_df %>%
  filter(padj < 0.05, abs(log2FC) >= 1) %>%
  arrange(padj, desc(abs(log2FC)))

sig_count <- nrow(sig_df)
up_count <- sum(sig_df$de_status == "Tumor-higher")
down_count <- sum(sig_df$de_status == "Normal-higher")

label_candidates <- sig_df
if ("baseMean" %in% colnames(label_candidates)) {
  label_candidates_high_expr <- label_candidates %>%
    filter(!is.na(baseMean), baseMean >= 10)
  if (nrow(label_candidates_high_expr) > 0) {
    label_candidates <- label_candidates_high_expr
  }
}

label_df <- bind_rows(
  label_candidates %>%
    filter(de_status == "Tumor-higher") %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice_head(n = 4),
  label_candidates %>%
    filter(de_status == "Normal-higher") %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice_head(n = 4)
) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(label_display = label)

labeled_genes <- label_df$label_display
if (length(labeled_genes) == 0) {
  labeled_genes <- character(0)
}

mapped_symbol_count <- sum(!is.na(volcano_df$SYMBOL) & volcano_df$SYMBOL != "")

subtitle_text <- paste0(
  "Primary reporting rule: padj<0.05 and |shrunken log2FC|>=1; ",
  "Tumor higher: ", up_count, "; Normal higher: ", down_count
)

help_items <- c(
  "Dot = one gene",
  "Right = higher in Tumor; Left = higher in Normal (effect size)",
  "Up = more significant (-log10 FDR)",
  "padj = FDR (false discovery rate)",
  "|log2FC|>=1 = >=2x change",
  "Labels use gene symbols when available",
  "Labels = top significant genes",
  "Points beyond +/-8 are clipped for display."
)

help_text <- paste(
  vapply(
    help_items,
    function(item) {
      wrapped <- strwrap(item, width = 34, initial = "- ", exdent = 2)
      paste(wrapped, collapse = "\n")
    },
    character(1)
  ),
  collapse = "\n"
)

palette_status <- c(
  "Not significant" = "#BDBDBD",
  "FDR only" = "#7F8C8D",
  "Normal-higher" = "#1B9E77",
  "Tumor-higher" = "#D95F02"
)

main_plot <- ggplot(volcano_df, aes(x = lfc_plot, y = neg_log10_padj)) +
  geom_vline(xintercept = c(-1, 1), color = "#6D6D6D", linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = -log10(0.05), color = "#6D6D6D", linetype = "dashed", linewidth = 0.45) +
  geom_point(
    data = volcano_df %>% filter(de_status == "Not significant"),
    color = palette_status[["Not significant"]],
    alpha = 0.25,
    size = 0.6
  ) +
  geom_point(
    data = volcano_df %>% filter(de_status != "Not significant"),
    aes(color = de_status),
    alpha = 0.7,
    size = 0.8
  ) +
  scale_color_manual(
    values = palette_status[c("FDR only", "Normal-higher", "Tumor-higher")],
    name = "DE status"
  ) +
  labs(
    title = "Effect size and adjusted statistical evidence",
    subtitle = subtitle_text,
    x = "log2 fold change (Tumor vs Normal)",
    y = "-log10(adjusted p-value)",
    caption = paste(
      "Labels are selected symmetrically from both directions using statistical evidence,\n",
      "not prior biological preference; labels are not validated biomarkers. Values beyond +/-8 are clipped."
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 8, margin = margin(t = 8)),
    plot.margin = margin(10, 8, 10, 10)
  ) +
  coord_cartesian(clip = "off")

if (requireNamespace("ggrepel", quietly = TRUE) && nrow(label_df) > 0) {
  main_plot <- main_plot +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label_display),
      size = 2.3,
      box.padding = 0.25,
      point.padding = 0.1,
      min.segment.length = 0,
      max.overlaps = 20
    )
}

if (FALSE && requireNamespace("patchwork", quietly = TRUE)) {
  help_panel <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 1,
      label = "How to read this figure",
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 4.0
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.9,
      label = help_text,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      lineheight = 1.1
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(10, 8, 10, 0))

  final_plot <- main_plot + help_panel + patchwork::plot_layout(widths = c(5.0, 1.2))
} else {
  final_plot <- main_plot
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 10, height = 6, dpi = 320)
ggsave(output_pdf_path, plot = final_plot, width = 10, height = 6, device = "pdf")

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F04",
  filename = "F04_de_volcano_paired_v2.png",
  purpose = "DE summary: volcano plot (effect size vs significance) for paired Tumor vs Normal",
  inputs = "results_v2/deseq2/deseq2_paired_v2_results.tsv",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F04")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("LFC column used:", lfc_col, "\n")
cat("number of mapped symbols (non-NA SYMBOL):", mapped_symbol_count, "\n")
cat(
  "labeled gene names list (as displayed):",
  if (length(labeled_genes) > 0) paste(labeled_genes, collapse = ", ") else "none",
  "\n"
)
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

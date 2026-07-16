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
output_path <- file.path(figure_dir, "F03_de_ma_paired_v2.png")
output_pdf_path <- file.path(vector_dir, "F03_de_ma_paired_v2.pdf")
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

de <- readr::read_tsv(de_path, show_col_types = FALSE)
assert_columns(de, c("gene_id", "baseMean", "padj", "log2FoldChange"), "deseq2_paired_v2_results.tsv")

lfc_col <- "log2FoldChange"
if ("log2FoldChange_shrunk" %in% colnames(de) && !all(is.na(de$log2FoldChange_shrunk))) {
  lfc_col <- "log2FoldChange_shrunk"
}

ma_df <- de %>%
  mutate(log2FC = .data[[lfc_col]]) %>%
  filter(!is.na(baseMean), !is.na(log2FC)) %>%
  mutate(
    x = log10(baseMean + 1),
    log2FC_plot = pmax(pmin(log2FC, 8), -8),
    de_status = case_when(
      !is.na(padj) & padj < 0.05 & log2FC >= 1 ~ "Up in tumor",
      !is.na(padj) & padj < 0.05 & log2FC <= -1 ~ "Up in normal",
      TRUE ~ "Not significant"
    )
  )

if (nrow(ma_df) == 0) {
  stop("No rows available for MA plot after filtering non-missing baseMean/log2FC.")
}

sig_count <- ma_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) >= 1) %>%
  nrow()

up_count <- sum(ma_df$de_status == "Up in tumor")
down_count <- sum(ma_df$de_status == "Up in normal")

subtitle_text <- paste0(
  "Primary reporting rule: padj<0.05 and |shrunken log2FC|>=1; ",
  "Tumor higher: ", up_count, "; Normal higher: ", down_count
)

help_items <- c(
  "Dot = one gene",
  "Up/down = higher in Tumor/Normal",
  "baseMean = mean expression (average abundance)",
  "log2FC = effect size (Tumor vs Normal)",
  "padj = FDR (false discovery rate)",
  "|log2FC|>=1 = >=2x change"
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
  "Up in normal" = "#1B9E77",
  "Up in tumor" = "#D95F02"
)

main_plot <- ggplot(ma_df, aes(x = x, y = log2FC_plot)) +
  geom_hline(yintercept = 0, color = "#4A4A4A", linewidth = 0.4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "#6D6D6D", linewidth = 0.4) +
  geom_point(
    data = ma_df %>% filter(de_status == "Not significant"),
    color = palette_status[["Not significant"]],
    alpha = 0.25,
    size = 0.6
  ) +
  geom_point(
    data = ma_df %>% filter(de_status != "Not significant"),
    aes(color = de_status),
    alpha = 0.7,
    size = 0.75
  ) +
  scale_color_manual(
    values = palette_status[c("Up in normal", "Up in tumor")],
    name = "DE status"
  ) +
  labs(
    title = "MA plot: paired Tumor vs Normal",
    subtitle = subtitle_text,
    x = "Mean expression (baseMean; log10 scale)",
    y = "Shrunken log2 fold change (Tumor vs Normal)",
    caption = paste(
      "Display values are clipped at +/-8. Statistical association does not establish\n",
      "biological mechanism, clinical relevance, or biomarker validity."
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 8, margin = margin(t = 8)),
    plot.margin = margin(10, 8, 10, 10)
  ) +
  coord_cartesian(clip = "off")

if (requireNamespace("patchwork", quietly = TRUE)) {
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
  xr <- range(ma_df$x, na.rm = TRUE)
  yr <- range(ma_df$log2FC_plot, na.rm = TRUE)
  x_annot <- xr[2] + 0.28 * (xr[2] - xr[1])
  y_annot <- yr[2]

  final_plot <- main_plot +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.42))) +
    annotate(
      "label",
      x = x_annot,
      y = y_annot,
      label = paste("How to read this figure\n", help_text),
      hjust = 0,
      vjust = 1,
      size = 3.0,
      lineheight = 1.1,
      label.size = 0.2,
      fill = "white",
      color = "#222222"
    ) +
    theme(plot.margin = margin(10, 140, 10, 10))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 10, height = 6, dpi = 320)
ggsave(output_pdf_path, plot = final_plot, width = 10, height = 6, device = "pdf")

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F03",
  filename = "F03_de_ma_paired_v2.png",
  purpose = "DE summary: MA plot (log2FC vs mean expression) for paired Tumor vs Normal",
  inputs = "results_v2/deseq2/deseq2_paired_v2_results.tsv",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F03")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("significant genes count (padj<0.05 & |LFC|>=1):", sig_count, "\n")
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

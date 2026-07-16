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

gsea_path <- "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
figure_dir <- "figures_v2/final"
vector_dir <- "figures_v2/vector"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F06_bio_hallmark_nes_paired_v2.png")
output_pdf_path <- file.path(vector_dir, "F06_bio_hallmark_nes_paired_v2.pdf")
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

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", as.character(x))
  x <- gsub("_", " ", x)
  x
}

gsea <- readr::read_tsv(gsea_path, show_col_types = FALSE)
assert_columns(gsea, c("pathway", "NES", "padj"), "hallmark_gsea_paired_v2.tsv")

gsea_clean <- gsea %>%
  filter(!is.na(NES), !is.na(padj))

if (nrow(gsea_clean) == 0) {
  stop("No Hallmark rows available after removing NA NES/padj.")
}

sig_tbl <- gsea_clean %>% filter(padj < 0.05)
n_sig <- nrow(sig_tbl)

fallback_mode <- FALSE
if (n_sig >= 16) {
  pos <- sig_tbl %>%
    filter(NES > 0) %>%
    arrange(padj, desc(NES)) %>%
    slice_head(n = 8)
  neg <- sig_tbl %>%
    filter(NES < 0) %>%
    arrange(padj, NES) %>%
    slice_head(n = 8)

  plot_tbl <- bind_rows(pos, neg) %>%
    distinct(pathway, .keep_all = TRUE)

  if (nrow(plot_tbl) < 16) {
    fill_tbl <- sig_tbl %>%
      filter(!pathway %in% plot_tbl$pathway) %>%
      arrange(padj, desc(abs(NES))) %>%
      slice_head(n = 16 - nrow(plot_tbl))
    plot_tbl <- bind_rows(plot_tbl, fill_tbl)
  }
} else {
  fallback_mode <- TRUE
  plot_tbl <- gsea_clean %>%
    arrange(desc(abs(NES)), padj) %>%
    slice_head(n = 16)
}

if (nrow(plot_tbl) == 0) {
  stop("No pathways selected for plotting.")
}

plot_tbl <- plot_tbl %>%
  mutate(
    direction = ifelse(NES > 0, "Toward Tumor-higher", "Toward Normal-higher"),
    pathway_clean = clean_pathway(pathway),
    fdr_label = paste0("FDR ", scales::scientific(padj, digits = 2))
  ) %>%
  arrange(NES) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

subtitle_suffix <- if (fallback_mode) {
  "Few significant pathways; showing top pathways by |NES|."
} else {
  "Showing significant pathways only (balanced by direction)."
}

subtitle_text <- paste(
  "NES (enrichment strength); padj (FDR); positive = Tumor-higher side, negative = Normal-higher side;\n",
  subtitle_suffix
)

help_items <- c(
  "Bar = one Hallmark pathway",
  "NES = enrichment strength (pathway shift)",
  "NES>0 = toward Tumor-higher; NES<0 = toward Normal-higher",
  "padj = FDR (false discovery rate)",
  "Bigger |NES| = stronger program-level difference"
)
if (fallback_mode) {
  help_items <- c(help_items, "If few are significant, we show top pathways by |NES|.")
}

help_text_raw <- paste(
  vapply(
    help_items,
    function(item) {
      wrapped <- strwrap(item, width = 30, initial = "- ", exdent = 2)
      paste(wrapped, collapse = "\n")
    },
    character(1)
  ),
  collapse = "\n"
)
help_text <- help_text_raw

dir_colors <- c(
  "Toward Tumor-higher" = "#D95F02",
  "Toward Normal-higher" = "#1B9E77"
)

main_plot <- ggplot(plot_tbl, aes(x = pathway_clean, y = NES, fill = direction)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = fdr_label), hjust = ifelse(plot_tbl$NES > 0, -0.05, 1.05), size = 2.5) +
  geom_hline(yintercept = 0, color = "#4A4A4A", linewidth = 0.5) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = dir_colors, name = "Direction") +
  labs(
    title = "Hallmark pathways enriched along the paired DE ranking",
    subtitle = subtitle_text,
    x = NULL,
    y = "Normalized enrichment score (NES)",
    caption = paste(
      "Positive NES indicates enrichment toward Tumor-higher genes; negative NES indicates\n",
      "enrichment toward Normal-higher genes. Enrichment supports hypotheses, not mechanism."
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 8, margin = margin(t = 8)),
    plot.margin = margin(10, 8, 10, 10)
  )

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
      size = 4.3
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.9,
      label = help_text,
      hjust = 0,
      vjust = 1,
      size = 3.3,
      lineheight = 1.1
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(8, 8, 8, 0))

  final_plot <- main_plot + help_panel + patchwork::plot_layout(widths = c(5.0, 1.5))
  final_plot <- final_plot & theme(plot.margin = margin(10, 20, 10, 10))
} else {
  final_plot <- main_plot + scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 11, height = 6.5, dpi = 320)
ggsave(output_pdf_path, plot = final_plot, width = 11, height = 6.5, device = "pdf")

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F06",
  filename = "F06_bio_hallmark_nes_paired_v2.png",
  purpose = "Biology summary: Hallmark GSEA NES barplot (Tumor vs Normal, paired cohort)",
  inputs = "results_v2/enrichment/hallmark_gsea_paired_v2.tsv",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F06")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("pathways plotted:", nrow(plot_tbl), "\n")
cat("number significant (padj<0.05):", n_sig, "\n")
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

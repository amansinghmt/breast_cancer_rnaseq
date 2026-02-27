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

go_path <- "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
figure_dir <- "figures_v2/final"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F07_bio_go_bp_dotplot_paired_v2.png")
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

parse_gene_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  num <- suppressWarnings(vapply(parts, function(p) as.numeric(p[1]), numeric(1)))
  den <- suppressWarnings(vapply(parts, function(p) as.numeric(p[2]), numeric(1)))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

shorten_term <- function(x, max_chars = 60) {
  x <- as.character(x)
  ifelse(
    nchar(x) > max_chars,
    paste0(substr(x, 1, max_chars - 3), "..."),
    x
  )
}

go_df <- readr::read_tsv(go_path, show_col_types = FALSE)
assert_columns(go_df, c("Description", "GeneRatio", "Count", "p.adjust"), "go_bp_ora_paired_v2.tsv")

go_clean <- go_df %>%
  filter(
    !is.na(Description),
    !is.na(GeneRatio),
    !is.na(Count),
    !is.na(p.adjust)
  ) %>%
  mutate(
    gene_ratio = parse_gene_ratio(GeneRatio),
    neg_log10_fdr = -log10(pmax(p.adjust, .Machine$double.xmin))
  ) %>%
  filter(!is.na(gene_ratio), is.finite(gene_ratio), is.finite(neg_log10_fdr))

if (nrow(go_clean) == 0) {
  stop("No GO terms remain after required-field and GeneRatio parsing filters.")
}

n_sig <- sum(go_clean$p.adjust < 0.05, na.rm = TRUE)
sig_tbl <- go_clean %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)

if (nrow(sig_tbl) >= 15) {
  plot_tbl <- sig_tbl %>%
    slice_head(n = 15)
} else {
  plot_tbl <- go_clean %>%
    arrange(p.adjust) %>%
    slice_head(n = 15)
}

plot_tbl <- plot_tbl %>%
  mutate(
    term_short = shorten_term(Description, max_chars = 60)
  ) %>%
  arrange(p.adjust, desc(neg_log10_fdr)) %>%
  mutate(
    term_short = factor(term_short, levels = rev(unique(term_short)))
  )

subtitle_text <- "GeneRatio (fraction of DE genes); Count (#genes); p.adjust (FDR)"

main_plot <- ggplot(
  plot_tbl,
  aes(x = gene_ratio, y = term_short, size = Count, color = neg_log10_fdr)
) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_gradient(
    low = "#9ECAE1",
    high = "#08519C",
    name = "-log10(FDR)"
  ) +
  labs(
    title = "GO Biological Process enrichment (paired cohort)",
    subtitle = subtitle_text,
    x = "GeneRatio (fraction of DE genes)",
    y = NULL,
    size = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 8, 10, 10)
  ) +
  coord_cartesian(clip = "off")

help_items <- c(
  "Dot = one GO term (biological process)",
  "GeneRatio = fraction of DE genes in the term",
  "Count = number of genes in the term",
  "Color = -log10(FDR) (higher = more signif.)",
  "This summarizes functions enriched in tumor vs normal"
)

help_text <- paste(
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
help_text <- paste(strwrap(help_text, width = 26), collapse = "\n")

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
      size = 4.3
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.9,
      label = help_text,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      lineheight = 1.1
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(8, 8, 8, 0))

  final_plot <- main_plot + help_panel + patchwork::plot_layout(widths = c(5.0, 1.6))
  final_plot <- final_plot & theme(plot.margin = margin(10, 24, 10, 10))
} else {
  x_range <- range(plot_tbl$gene_ratio, na.rm = TRUE)
  x_annot <- x_range[2] + 0.18 * (x_range[2] - x_range[1])

  final_plot <- main_plot +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.45))
    ) +
    annotate(
      "label",
      x = x_annot,
      y = 1,
      label = paste("How to read this figure\n", help_text),
      hjust = 0,
      vjust = 1,
      size = 3.0,
      lineheight = 1.1,
      label.size = 0.2,
      fill = "white",
      color = "#222222"
    ) +
    theme(plot.margin = margin(10, 130, 10, 10))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 11, height = 6.5, dpi = 320)

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F07",
  filename = "F07_bio_go_bp_dotplot_paired_v2.png",
  purpose = "Biology summary: GO BP over-representation dotplot for paired cohort DE genes",
  inputs = "results_v2/enrichment/go_bp_ora_paired_v2.tsv",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F07")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("terms plotted:", nrow(plot_tbl), "\n")
cat("significant terms (p.adjust<0.05):", n_sig, "\n")
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

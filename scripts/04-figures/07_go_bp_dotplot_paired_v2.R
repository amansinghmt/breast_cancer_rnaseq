#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

tumor_path <- "results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv"
normal_path <- "results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv"
figure_dir <- "figures_v2/final"
vector_dir <- "figures_v2/vector"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F07_bio_go_bp_dotplot_paired_v2.png")
output_pdf_path <- file.path(vector_dir, "F07_bio_go_bp_dotplot_paired_v2.pdf")
fig_manifest_path <- file.path(results_dir, "fig_manifest.tsv")

assert_columns <- function(df, cols, label) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(label, " missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  numerator <- suppressWarnings(vapply(parts, function(p) as.numeric(p[[1]]), numeric(1)))
  denominator <- suppressWarnings(vapply(parts, function(p) as.numeric(p[[2]]), numeric(1)))
  ratio <- numerator / denominator
  ratio[!is.finite(ratio)] <- NA_real_
  ratio
}

prepare_direction <- function(path, direction, panel_label) {
  tbl <- readr::read_tsv(path, show_col_types = FALSE)
  assert_columns(tbl, c("Description", "GeneRatio", "Count", "p.adjust"), basename(path))
  tbl %>%
    filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    mutate(
      direction = direction,
      panel = panel_label,
      gene_ratio = parse_ratio(GeneRatio),
      evidence = -log10(pmax(p.adjust, .Machine$double.xmin)),
      presentation_family = case_when(
        grepl("nuclear division|chromosome segregation|chromosome separation", Description) ~ "division and segregation",
        grepl("spindle|microtubule", Description) ~ "spindle and microtubule",
        grepl("nucleosome|protein-DNA|chromosome organization|chromosome condensation", Description) ~ "chromatin organization",
        grepl("cell cycle|G2/M", Description) ~ "cell cycle",
        grepl("circulation|blood pressure", Description) ~ "circulation",
        grepl("vasculature|vascular|endothelial", Description) ~ "vascular biology",
        grepl("muscle|action potential", Description) ~ "muscle and excitability",
        grepl("lipid|fat cell|triglyceride|ketone", Description) ~ "lipid metabolism",
        grepl("development|differentiation|maturation|gliogenesis", Description) ~ "development and differentiation",
        grepl("response to", Description) ~ "response processes",
        TRUE ~ Description
      )
    ) %>%
    filter(is.finite(gene_ratio), is.finite(evidence)) %>%
    arrange(p.adjust, desc(Count), Description) %>%
    group_by(presentation_family) %>%
    slice_head(n = 2) %>%
    ungroup() %>%
    arrange(p.adjust, desc(Count), Description) %>%
    slice_head(n = 10)
}

tumor <- prepare_direction(tumor_path, "Tumor-higher", "A  Tumor-higher gene set")
normal <- prepare_direction(normal_path, "Normal-higher", "B  Normal-higher gene set")
if (nrow(tumor) == 0 && nrow(normal) == 0) {
  stop("Neither directional GO analysis has significant representative terms.")
}

plot_tbl <- bind_rows(tumor, normal) %>%
  mutate(term_key = paste(panel, Description, sep = "___")) %>%
  arrange(panel, desc(p.adjust)) %>%
  mutate(term_key = factor(term_key, levels = unique(term_key)))

direction_colors <- c("Tumor-higher" = "#D95F02", "Normal-higher" = "#1B9E77")

final_plot <- ggplot(
  plot_tbl,
  aes(x = gene_ratio, y = term_key, size = Count, color = direction, alpha = evidence)
) +
  geom_point() +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = direction_colors, guide = "none") +
  scale_alpha_continuous(name = "-log10(FDR)", range = c(0.45, 1)) +
  scale_size_continuous(name = "Gene count", range = c(2.5, 7)) +
  labs(
    title = "Directional GO Biological Process over-representation",
    subtitle = paste(
      "Ten significant representatives per direction; Wang cutoff 0.7; max two terms per",
      "documented keyword family; BH FDR < 0.05"
    ),
    x = "Gene ratio within the directional input list",
    y = NULL,
    caption = paste(
      "Both analyses use the same 30,244-gene tested universe (15,233 GO-annotated).",
      "ORA identifies over-representation, not mechanism."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 8),
    legend.position = "top",
    plot.caption = element_text(hjust = 0, size = 8, margin = margin(t = 8)),
    plot.margin = margin(10, 10, 10, 10)
  )

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vector_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 13, height = 7.5, dpi = 320)
ggsave(output_pdf_path, plot = final_plot, width = 13, height = 7.5, device = "pdf")

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F07",
  filename = "F07_bio_go_bp_dotplot_paired_v2.png",
  purpose = "Directional GO BP over-representation for Tumor-higher and Normal-higher genes",
  inputs = paste(tumor_path, normal_path, sep = "; "),
  stringsAsFactors = FALSE
)
if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F07") %>%
    bind_rows(new_row)
} else {
  fig_manifest <- new_row
}
readr::write_tsv(fig_manifest, fig_manifest_path)

cat("Tumor-higher representative terms plotted:", nrow(tumor), "\n")
cat("Normal-higher representative terms plotted:", nrow(normal), "\n")
cat("output file size:", file.info(output_path)$size, "\n")

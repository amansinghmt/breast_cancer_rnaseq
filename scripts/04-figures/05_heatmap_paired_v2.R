#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "scales")
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
  library(tidyr)
  library(ggplot2)
  library(scales)
})

manifest_path <- "data/metadata/sample_manifest.tsv"
de_path <- "results_v2/differential_expression/deseq2_paired_v2_results.tsv"
vst_path <- "results_v2/differential_expression/deseq2_paired_v2_vst.tsv"
log2cpm_path <- "results_v2/differential_expression/deseq2_paired_v2_log2cpm.tsv"
counts_path <- "data/processed/counts.tsv"

figure_dir <- "figures_v2/final"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F05_de_heatmap_top40_paired_v2.png")
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

zscore_row <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x)) / s
}

manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
assert_columns(manifest, c("sample_id", "patient_id", "condition_main", "include_paired"), "sample_manifest.tsv")

paired <- manifest %>%
  mutate(include_paired_chr = toupper(as.character(include_paired))) %>%
  filter(include_paired_chr == "TRUE", condition_main %in% c("Normal", "Tumor")) %>%
  select(sample_id, patient_id, condition_main)

pair_check <- paired %>%
  group_by(patient_id) %>%
  summarize(
    n_samples = n(),
    n_normal = sum(condition_main == "Normal"),
    n_tumor = sum(condition_main == "Tumor"),
    .groups = "drop"
  )

bad_pairs <- pair_check %>% filter(!(n_samples == 2 & n_normal == 1 & n_tumor == 1))
if (nrow(bad_pairs) > 0) {
  bad_msg <- paste(
    apply(as.data.frame(bad_pairs), 1, function(x) {
      paste0(
        x[["patient_id"]],
        "(n=", x[["n_samples"]],
        ", normal=", x[["n_normal"]],
        ", tumor=", x[["n_tumor"]],
        ")"
      )
    }),
    collapse = "; "
  )
  stop(paste0("Invalid paired structure: ", bad_msg))
}

patient_tbl <- paired %>%
  distinct(patient_id) %>%
  mutate(patient_num = suppressWarnings(as.numeric(patient_id)))

if (all(!is.na(patient_tbl$patient_num))) {
  patient_levels <- patient_tbl %>%
    arrange(patient_num) %>%
    pull(patient_id)
} else {
  patient_levels <- patient_tbl %>%
    arrange(patient_id) %>%
    pull(patient_id)
}

paired <- paired %>%
  mutate(
    patient_id = factor(as.character(patient_id), levels = patient_levels),
    condition_main = factor(as.character(condition_main), levels = c("Normal", "Tumor"))
  ) %>%
  arrange(patient_id, condition_main)

sample_order <- paired$sample_id
n_patients <- n_distinct(paired$patient_id)
n_samples <- nrow(paired)

de <- readr::read_tsv(de_path, show_col_types = FALSE)
assert_columns(de, c("gene_id", "padj"), "deseq2_paired_v2_results.tsv")

top_de <- de %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 40) %>%
  mutate(gene_id_clean = clean_gene_id(gene_id))

if (nrow(top_de) == 0) {
  stop("No genes available from DE results after filtering non-NA padj.")
}

label_priority <- c("gene_symbol", "gene_name", "symbol", "external_gene_name", "hgnc_symbol")
label_hits <- label_priority[label_priority %in% colnames(top_de)]
label_col <- if (length(label_hits) > 0) label_hits[[1]] else NA_character_

if (!is.na(label_col)) {
  top_de <- top_de %>%
    mutate(SYMBOL = trimws(as.character(.data[[label_col]])))
} else if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
           requireNamespace("AnnotationDbi", quietly = TRUE)) {
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(top_de$gene_id_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
  map_df <- map_df %>%
    select(ENSEMBL, SYMBOL) %>%
    filter(!is.na(SYMBOL), nzchar(SYMBOL)) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  top_de <- top_de %>%
    left_join(map_df, by = c("gene_id_clean" = "ENSEMBL"))
} else {
  top_de <- top_de %>%
    mutate(SYMBOL = NA_character_)
}

top_de <- top_de %>%
  mutate(
    SYMBOL = ifelse(trimws(as.character(SYMBOL)) == "", NA_character_, trimws(as.character(SYMBOL))),
    label = ifelse(!is.na(SYMBOL), SYMBOL, gene_id_clean)
  )

expr_source_used <- NA_character_
if (file.exists(vst_path)) {
  expr_source_used <- vst_path
} else if (file.exists(log2cpm_path)) {
  expr_source_used <- log2cpm_path
} else if (file.exists(counts_path)) {
  expr_source_used <- counts_path
} else {
  stop(
    paste0(
      "No expression source found. Checked: ",
      paste(c(vst_path, log2cpm_path, counts_path), collapse = ", ")
    )
  )
}

if (expr_source_used %in% c(vst_path, log2cpm_path)) {
  expr_df <- readr::read_tsv(expr_source_used, show_col_types = FALSE)
  assert_columns(expr_df, c("gene_id"), basename(expr_source_used))

  missing_samples <- setdiff(sample_order, colnames(expr_df))
  if (length(missing_samples) > 0) {
    stop(
      paste0(
        "Expression matrix missing paired sample columns: ",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  expr_gene <- clean_gene_id(expr_df$gene_id)
  keep_gene <- !duplicated(expr_gene)
  expr_mat <- as.matrix(expr_df[keep_gene, sample_order, drop = FALSE])
  rownames(expr_mat) <- expr_gene[keep_gene]
  storage.mode(expr_mat) <- "numeric"

  if (any(!is.finite(expr_mat))) {
    stop("Expression matrix contains NA/Inf values after conversion.")
  }
} else {
  counts_df <- readr::read_tsv(counts_path, show_col_types = FALSE)
  assert_columns(counts_df, c("gene_id"), "counts.tsv")

  missing_samples <- setdiff(sample_order, colnames(counts_df))
  if (length(missing_samples) > 0) {
    stop(
      paste0(
        "counts.tsv missing paired sample columns: ",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  counts_gene <- clean_gene_id(counts_df$gene_id)
  keep_gene <- !duplicated(counts_gene)
  counts_mat <- as.matrix(counts_df[keep_gene, sample_order, drop = FALSE])
  rownames(counts_mat) <- counts_gene[keep_gene]
  storage.mode(counts_mat) <- "numeric"

  if (any(!is.finite(counts_mat))) {
    stop("counts.tsv contains NA/Inf values.")
  }
  if (any(counts_mat < 0)) {
    stop("counts.tsv contains negative values.")
  }

  lib_sizes <- colSums(counts_mat)
  if (any(lib_sizes <= 0)) {
    stop("counts.tsv has non-positive library size in paired samples.")
  }

  cpm_mat <- sweep(counts_mat, 2, lib_sizes, "/") * 1e6
  expr_mat <- log2(cpm_mat + 1)
}

genes_requested <- unique(top_de$gene_id_clean)
genes_present <- genes_requested[genes_requested %in% rownames(expr_mat)]
if (length(genes_present) < 20) {
  stop(
    paste0(
      "Too few top genes found in expression matrix (",
      length(genes_present),
      "). Need at least 20."
    )
  )
}

gene_order_de <- top_de$gene_id_clean[top_de$gene_id_clean %in% genes_present]
gene_order_de <- unique(gene_order_de)
expr_top <- expr_mat[gene_order_de, sample_order, drop = FALSE]

z_mat <- t(apply(expr_top, 1, zscore_row))
rownames(z_mat) <- rownames(expr_top)
colnames(z_mat) <- colnames(expr_top)
z_plot_mat <- pmax(pmin(z_mat, 2.5), -2.5)

row_hc <- hclust(dist(z_mat))
row_order <- rownames(z_plot_mat)[row_hc$order]
z_plot_mat <- z_plot_mat[row_order, sample_order, drop = FALSE]

label_tbl <- top_de %>%
  distinct(gene_id_clean, .keep_all = TRUE) %>%
  select(gene_id_clean, label)

row_labels <- label_tbl$label[match(row_order, label_tbl$gene_id_clean)]
row_labels <- ifelse(is.na(row_labels) | row_labels == "", row_order, row_labels)
label_map <- setNames(row_labels, row_order)

heat_df <- as.data.frame(z_plot_mat) %>%
  mutate(gene_id_clean = rownames(.)) %>%
  tidyr::pivot_longer(
    cols = -gene_id_clean,
    names_to = "sample_id",
    values_to = "z_plot"
  ) %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    gene_id_clean = factor(gene_id_clean, levels = rev(row_order))
  )

anno_df <- paired %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    condition_main = factor(condition_main, levels = c("Normal", "Tumor")),
    y = 1
  )

pair_breaks <- seq(2.5, length(sample_order) - 0.5, by = 2)
if (length(pair_breaks) == 1 && pair_breaks < 2.5) {
  pair_breaks <- numeric(0)
}

cond_colors <- c("Normal" = "#1B9E77", "Tumor" = "#D95F02")

anno_plot <- ggplot(anno_df, aes(x = sample_id, y = y, fill = condition_main)) +
  geom_tile(height = 0.9) +
  {
    if (length(pair_breaks) > 0) geom_vline(xintercept = pair_breaks, color = "#FFFFFF", linewidth = 0.5)
  } +
  scale_fill_manual(values = cond_colors, name = "Condition") +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 6, 0, 10)
  )

heat_plot <- ggplot(heat_df, aes(x = sample_id, y = gene_id_clean, fill = z_plot)) +
  geom_tile() +
  {
    if (length(pair_breaks) > 0) geom_vline(xintercept = pair_breaks, color = "#FFFFFF", linewidth = 0.35)
  } +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Row z-score"
  ) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) label_map[x], expand = c(0, 0)) +
  labs(
    title = "Top 40 DE genes: paired heatmap (row z-scored)",
    subtitle = "Columns are paired samples ordered as Normal then Tumor per patient",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(0, 6, 10, 10)
  ) +
  coord_cartesian(clip = "off")

help_items <- c(
  "Row = gene; column = sample",
  "Color = z-score (scaled expression) [relative level]",
  "Red/high = above gene average; Blue/low = below gene average",
  "Columns are paired: Normal then Tumor per patient",
  "Use this to see patterns, not absolute counts",
  "Some rows show Ensembl IDs (ENSG...) when gene symbol is unavailable."
)

help_intro <- paste(
  strwrap("Heatmap = color map of gene expression patterns", width = 34),
  collapse = "\n"
)

help_bullets <- paste(
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

help_text <- paste(help_intro, "", help_bullets, sep = "\n")

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
      size = 2.8,
      lineheight = 1.1
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      plot.margin = margin(8, 8, 8, 0)
    )

  left_block <- anno_plot / heat_plot + patchwork::plot_layout(heights = c(0.12, 1))
  final_plot <- (left_block | help_panel) + patchwork::plot_layout(widths = c(5.2, 1.3))
} else {
  final_plot <- heat_plot +
    annotate(
      "label",
      x = sample_order[[length(sample_order)]],
      y = row_order[[1]],
      label = paste("How to read this figure\n", help_text),
      hjust = 1,
      vjust = 1,
      size = 2.8,
      lineheight = 1.1,
      label.size = 0.2,
      fill = "white",
      color = "#222222"
    ) +
    theme(plot.margin = margin(10, 90, 10, 10))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 11, height = 7, dpi = 320)

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F05",
  filename = "F05_de_heatmap_top40_paired_v2.png",
  purpose = "DE summary: heatmap of top DE genes (row z-scored expression) in paired cohort",
  inputs = "sample_manifest.tsv; deseq2_paired_v2_results.tsv; expression source used",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F05")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("expression source used:", expr_source_used, "\n")
cat("genes plotted:", nrow(z_plot_mat), "\n")
cat("paired patients used:", n_patients, "\n")
cat("paired samples used:", n_samples, "\n")
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

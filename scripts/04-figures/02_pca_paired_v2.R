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
figure_dir <- "figures_v2/final"
results_dir <- "results_v2"
output_path <- file.path(figure_dir, "F02_qc_pca_pairs.png")
fig_manifest_path <- file.path(results_dir, "fig_manifest.tsv")
counts_path <- "data/processed/counts.tsv"

expr_candidates <- c(
  "results_v2/deseq2/deseq2_paired_v2_vst.tsv",
  "results_v2/deseq2/deseq2_paired_v2_vst_matrix.tsv",
  "results_v2/deseq2/deseq2_paired_v2_log2cpm.tsv"
)

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

manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
assert_columns(
  manifest,
  c("sample_id", "patient_id", "condition_main", "include_paired"),
  "sample_manifest.tsv"
)

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

bad_pairs <- pair_check %>%
  filter(!(n_samples == 2 & n_normal == 1 & n_tumor == 1))

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
  stop(paste0("Invalid paired structure. Expected one Normal and one Tumor per patient: ", bad_msg))
}

patient_order_tbl <- paired %>%
  distinct(patient_id) %>%
  mutate(patient_num = suppressWarnings(as.numeric(patient_id)))

if (all(!is.na(patient_order_tbl$patient_num))) {
  patient_levels <- patient_order_tbl %>%
    arrange(patient_num) %>%
    pull(patient_id)
} else {
  patient_levels <- patient_order_tbl %>%
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

selected_expr <- expr_candidates[file.exists(expr_candidates)]
expr_source_type <- "counts_fallback"
expr_source_path <- counts_path

if (length(selected_expr) > 0) {
  expr_source_type <- "precomputed_log_matrix"
  expr_source_path <- selected_expr[[1]]
} else if (!file.exists(counts_path)) {
  stop(
    paste0(
      "No expression input found. Checked: ",
      paste(c(expr_candidates, counts_path), collapse = ", ")
    )
  )
}

if (expr_source_type == "precomputed_log_matrix") {
  expr_df <- readr::read_tsv(expr_source_path, show_col_types = FALSE)
  assert_columns(expr_df, c("gene_id"), paste0(basename(expr_source_path)))

  missing_samples <- setdiff(sample_order, colnames(expr_df))
  if (length(missing_samples) > 0) {
    stop(
      paste0(
        "Expression matrix missing paired sample columns: ",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  expr_mat <- as.matrix(expr_df[, sample_order, drop = FALSE])
  rownames(expr_mat) <- expr_df$gene_id
  storage.mode(expr_mat) <- "numeric"

  if (any(!is.finite(expr_mat))) {
    stop("Selected expression matrix contains NA/Inf values after numeric conversion.")
  }

  val <- as.numeric(expr_mat)
  frac_integer <- mean(abs(val - round(val)) < 1e-8)
  max_val <- max(val, na.rm = TRUE)
  if (max_val > 1000 || (frac_integer > 0.99 && max_val > 50)) {
    stop(
      paste0(
        "Selected matrix appears non-log (raw-count-like). ",
        "Please provide VST/log2-transformed matrix or use counts fallback."
      )
    )
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

  counts_mat <- as.matrix(counts_df[, sample_order, drop = FALSE])
  rownames(counts_mat) <- counts_df$gene_id
  storage.mode(counts_mat) <- "numeric"

  if (any(!is.finite(counts_mat))) {
    stop("counts.tsv contains NA/Inf values in paired sample columns.")
  }
  if (any(counts_mat < 0)) {
    stop("counts.tsv contains negative values.")
  }

  library_sizes <- colSums(counts_mat)
  if (any(library_sizes <= 0)) {
    stop("At least one paired sample has non-positive library size in counts.tsv.")
  }

  cpm_mat <- sweep(counts_mat, 2, library_sizes, "/") * 1e6
  expr_mat <- log2(cpm_mat + 1)
}

gene_var <- apply(expr_mat, 1, var, na.rm = TRUE)
keep <- is.finite(gene_var) & gene_var > 0
if (sum(keep) < 2) {
  stop("Insufficient non-zero-variance genes for PCA after filtering.")
}

expr_var <- expr_mat[keep, , drop = FALSE]
var_keep <- gene_var[keep]
top_n <- min(2000, nrow(expr_var))
top_idx <- order(var_keep, decreasing = TRUE)[seq_len(top_n)]
expr_top <- expr_var[top_idx, , drop = FALSE]

pca <- prcomp(t(expr_top), center = TRUE, scale. = FALSE)
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
var1 <- var_explained[1]
var2 <- var_explained[2]

scores <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(paired %>% mutate(patient_id = as.character(patient_id)), by = "sample_id")

seg_df <- scores %>%
  select(patient_id, condition_main, PC1, PC2) %>%
  tidyr::pivot_wider(
    names_from = condition_main,
    values_from = c(PC1, PC2),
    names_sep = "_"
  )

if (any(!is.finite(seg_df$PC1_Normal)) || any(!is.finite(seg_df$PC1_Tumor)) ||
    any(!is.finite(seg_df$PC2_Normal)) || any(!is.finite(seg_df$PC2_Tumor))) {
  stop("Could not build complete Normal/Tumor pair segments for PCA plot.")
}

# QC: per-patient paired distance in PCA space.
qc_dir <- file.path(results_dir, "qc")
pca_pair_dist_path <- file.path(qc_dir, "pca_pair_distances_v2.tsv")

pca_pair_dist <- seg_df %>%
  mutate(
    dist = sqrt((PC1_Tumor - PC1_Normal)^2 + (PC2_Tumor - PC2_Normal)^2)
  ) %>%
  select(patient_id, dist) %>%
  arrange(desc(dist))

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(pca_pair_dist, pca_pair_dist_path)

cat("Top 5 patients by Normal-Tumor PCA distance:\n")
print(pca_pair_dist %>% slice_head(n = 5))

palette_condition <- c("Normal" = "#1B9E77", "Tumor" = "#D95F02")

subtitle_text <- "n=21 patients (42 samples); include_paired==TRUE\nTop 2,000 variable genes"

help_items <- c(
  "Each dot = one sample",
  "Line = matched pair (paired design = same patient)",
  "PCA (dimension reduction = compress many genes into 2 axes)",
  "PC1/PC2 = main patterns (largest variance)",
  "If Tumor vs Normal separate -> global expression difference",
  "Outliers (far points) may indicate batch (technical effect)"
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

main_plot <- ggplot() +
  geom_segment(
    data = seg_df,
    aes(x = PC1_Normal, y = PC2_Normal, xend = PC1_Tumor, yend = PC2_Tumor),
    color = "#B8B8B8",
    linewidth = 0.4,
    alpha = 0.9
  ) +
  geom_point(
    data = scores,
    aes(x = PC1, y = PC2, color = condition_main),
    size = 2.8,
    alpha = 0.95
  ) +
  scale_color_manual(values = palette_condition, name = "Condition") +
  labs(
    title = "PCA of paired samples (QC)",
    subtitle = subtitle_text,
    x = sprintf("PC1 (%.1f%% variance)", var1),
    y = sprintf("PC2 (%.1f%% variance)", var2)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
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

  final_plot <- main_plot + help_panel + patchwork::plot_layout(widths = c(5.0, 1.3))
} else {
  xr <- range(scores$PC1, na.rm = TRUE)
  yr <- range(scores$PC2, na.rm = TRUE)
  x_annot <- xr[2] + 0.32 * (xr[2] - xr[1])
  y_annot <- yr[2]

  final_plot <- main_plot +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.45))) +
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
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = final_plot, width = 10, height = 6, dpi = 320)

manifest_cols <- c("figure_id", "filename", "purpose", "inputs")
new_row <- data.frame(
  figure_id = "F02",
  filename = "F02_qc_pca_pairs.png",
  purpose = "QC: PCA of paired cohort with patient-matched pair lines (Normal vs Tumor)",
  inputs = "data/metadata/sample_manifest.tsv; (expression matrix or counts.tsv as used)",
  stringsAsFactors = FALSE
)

if (file.exists(fig_manifest_path)) {
  fig_manifest <- readr::read_tsv(fig_manifest_path, show_col_types = FALSE)
  assert_columns(fig_manifest, manifest_cols, "fig_manifest.tsv")
  fig_manifest <- fig_manifest %>%
    select(all_of(manifest_cols)) %>%
    filter(figure_id != "F02")
  fig_manifest <- bind_rows(fig_manifest, new_row)
} else {
  fig_manifest <- new_row
}

readr::write_tsv(fig_manifest, fig_manifest_path)

out_exists <- file.exists(output_path)
out_size <- if (out_exists) file.info(output_path)$size else NA_integer_

cat("paired patients used:", n_patients, "\n")
cat("paired samples used:", n_samples, "\n")
cat(sprintf("PCA variance explained PC1, PC2: %.1f%%, %.1f%%\n", var1, var2))
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

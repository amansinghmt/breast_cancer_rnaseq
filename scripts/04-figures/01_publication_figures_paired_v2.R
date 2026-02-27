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
qc_path <- "results/qc/qc_summary.tsv"
figure_dir <- "figures_v2/final"
results_dir <- "results_v2"
figure_path <- file.path(figure_dir, "F01_qc_library_size_pairs.png")
fig_manifest_path <- file.path(results_dir, "fig_manifest.tsv")

qc_threshold <- 1e6
subtitle_text <- "n=21 patients (42 samples); include_paired==TRUE; threshold 1e6"

required_manifest_cols <- c("sample_id", "patient_id", "condition_main", "include_paired")
preferred_library_cols <- c(
  "library_size",
  "lib_size",
  "total_reads",
  "n_reads",
  "total_counts",
  "sum_counts"
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

detect_library_column <- function(qc_df) {
  qc_cols <- colnames(qc_df)
  exact_hits <- preferred_library_cols[preferred_library_cols %in% qc_cols]
  if (length(exact_hits) > 0) {
    return(exact_hits[[1]])
  }

  grep_hits <- qc_cols[grepl("lib|read|count", qc_cols, ignore.case = TRUE)]
  if (length(grep_hits) > 0) {
    return(grep_hits[[1]])
  }

  stop(
    paste0(
      "Could not detect library-size column in qc_summary.tsv. Available columns: ",
      paste(qc_cols, collapse = ", ")
    )
  )
}

manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
qc_summary <- readr::read_tsv(qc_path, show_col_types = FALSE)

assert_columns(manifest, required_manifest_cols, "sample_manifest.tsv")
assert_columns(qc_summary, c("sample_id"), "qc_summary.tsv")

library_col <- detect_library_column(qc_summary)

paired_manifest <- manifest %>%
  mutate(include_paired_chr = toupper(as.character(include_paired))) %>%
  filter(include_paired_chr == "TRUE", condition_main %in% c("Normal", "Tumor")) %>%
  dplyr::select(sample_id, patient_id, condition_main)

plot_df <- paired_manifest %>%
  left_join(
    qc_summary %>% dplyr::select(sample_id, library_size_raw = all_of(library_col)),
    by = "sample_id"
  ) %>%
  mutate(
    patient_id = as.character(patient_id),
    condition_main = as.character(condition_main),
    library_size = suppressWarnings(as.numeric(library_size_raw))
  )

if (nrow(plot_df) == 0) {
  stop("No rows after filtering include_paired == TRUE with Normal/Tumor conditions.")
}

if (any(is.na(plot_df$library_size))) {
  bad_ids <- plot_df %>% filter(is.na(library_size)) %>% pull(sample_id)
  stop(
    paste0(
      "Non-numeric or missing library_size values after join for sample_id(s): ",
      paste(bad_ids, collapse = ", ")
    )
  )
}

pair_check <- plot_df %>%
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
  stop(paste0("Invalid paired structure. Expected exactly one Normal and one Tumor per patient: ", bad_msg))
}

patient_sort_df <- plot_df %>%
  group_by(patient_id) %>%
  summarize(pair_median = median(library_size), .groups = "drop") %>%
  mutate(patient_id_num = suppressWarnings(as.numeric(patient_id)))

if (all(!is.na(patient_sort_df$patient_id_num))) {
  patient_levels <- patient_sort_df %>%
    arrange(patient_id_num) %>%
    pull(patient_id)
} else {
  patient_levels <- patient_sort_df %>%
    arrange(pair_median, patient_id) %>%
    pull(patient_id)
}

plot_df <- plot_df %>%
  mutate(
    patient_id = factor(patient_id, levels = patient_levels),
    condition_main = factor(condition_main, levels = c("Normal", "Tumor"))
  ) %>%
  arrange(patient_id, condition_main)

segment_df <- plot_df %>%
  dplyr::select(patient_id, condition_main, library_size) %>%
  tidyr::pivot_wider(names_from = condition_main, values_from = library_size)

palette_condition <- c("Normal" = "#1B9E77", "Tumor" = "#D95F02")

base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

label_df <- data.frame(
  x = qc_threshold * 1.05,
  y = patient_levels[[length(patient_levels)]],
  label = "QC min = 1e6",
  stringsAsFactors = FALSE
)

help_items <- c(
  "Dot = one sample",
  "Line = matched pair (paired design = same patient)",
  "Library size (sequencing depth = total reads)",
  "Log10 scale (log scale)",
  "Dashed line = QC minimum (minimum acceptable reads) = 1e6",
  "If a dot is left of dashed line -> low-depth sample (risky)"
)

wrap_width <- 30
help_text <- paste(
  vapply(
    help_items,
    function(item) {
      wrapped <- strwrap(item, width = wrap_width, initial = "- ", exdent = 2)
      paste(wrapped, collapse = "\n")
    },
    character(1)
  ),
  collapse = "\n"
)

main_plot <- ggplot() +
  geom_segment(
    data = segment_df,
    aes(
      x = Normal,
      xend = Tumor,
      y = patient_id,
      yend = patient_id
    ),
    color = "#B8B8B8",
    linewidth = 0.45
  ) +
  geom_point(
    data = plot_df,
    aes(x = library_size, y = patient_id, color = condition_main),
    size = 2.8,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = qc_threshold,
    linetype = "dashed",
    color = "#444444",
    linewidth = 0.6
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    hjust = 0,
    vjust = -0.7,
    size = 4,
    color = "#333333"
  ) +
  scale_color_manual(values = palette_condition, name = "Condition") +
  scale_x_log10(
    breaks = c(1e6, 3e6, 1e7, 3e7, 1e8),
    labels = c("1e6", "3e6", "1e7", "3e7", "1e8"),
    limits = c(NA, 1e8),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Library size per sample (paired cohort)",
    subtitle = subtitle_text,
    x = "Library size (reads, log10 scale)",
    y = "Patient ID"
  ) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(plot.margin = margin(10, 6, 10, 10))

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
      size = 4.2
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
    theme(
      plot.margin = margin(10, 10, 10, 0)
    )

  f01_plot <- main_plot + help_panel + patchwork::plot_layout(widths = c(4.8, 1.2))
} else {
  help_x <- max(plot_df$library_size, na.rm = TRUE) * 1.35
  help_y <- patient_levels[[length(patient_levels)]]

  f01_plot <- main_plot +
    annotate(
      "label",
      x = help_x,
      y = help_y,
      label = paste("How to read this figure\n", help_text),
      hjust = 0,
      vjust = 1,
      size = 3.0,
      label.size = 0.2,
      fill = "white",
      color = "#222222",
      lineheight = 1.1
    )
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(figure_path, plot = f01_plot, width = 10, height = 6, dpi = 320)

fig_manifest <- data.frame(
  figure_id = "F01",
  filename = "F01_qc_library_size_pairs.png",
  purpose = "QC: paired cohort library size dumbbell plot with 1e6 threshold",
  inputs = "data/metadata/sample_manifest.tsv; results/qc/qc_summary.tsv",
  stringsAsFactors = FALSE
)
readr::write_tsv(fig_manifest, fig_manifest_path)

paired_patients <- n_distinct(plot_df$patient_id)
paired_samples <- nrow(plot_df)
lib_min <- min(plot_df$library_size, na.rm = TRUE)
lib_median <- median(plot_df$library_size, na.rm = TRUE)
lib_max <- max(plot_df$library_size, na.rm = TRUE)
below_thresh <- sum(plot_df$library_size < qc_threshold, na.rm = TRUE)
out_exists <- file.exists(figure_path)
out_size <- if (out_exists) file.info(figure_path)$size else NA_integer_

cat("paired patients used:", paired_patients, "\n")
cat("paired samples used:", paired_samples, "\n")
cat("library_size min/median/max:", lib_min, "/", lib_median, "/", lib_max, "\n")
cat("below 1e6:", below_thresh, "\n")
cat("output file exists:", out_exists, "\n")
cat("output file size:", out_size, "\n")

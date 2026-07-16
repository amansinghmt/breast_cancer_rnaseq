#!/usr/bin/env Rscript

# Purpose:
#   Run pathway-level interpretation for paired DE results using:
#   - Hallmark GSEA (program-level directional enrichment),
#   - GO Biological Process ORA (over-representation of significant genes).
#
# Inputs:
#   - results_v2/deseq2/deseq2_paired_v2_results.tsv
#
# Outputs:
#   - results_v2/enrichment/hallmark_gsea_paired_v2.tsv
#   - results_v2/enrichment/go_bp_ora_paired_v2.tsv
#   - results_v2/enrichment/go_bp_ora_representative_v2.tsv
#   - results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv
#   - results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv
#   - results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv
#   - results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv
#   - results_v2/enrichment/enrichment_diagnostics_v2.tsv
#   - results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt
#
# Determinism / reproducibility:
#   - Seed set from PIPELINE_SEED (default 20260227 from run_v2.sh).
#   - Fixed GSEA implementation: fgseaMultilevel + BiocParallel::SerialParam().
#   - Stable ranking tie handling: deterministic secondary sort + tiny epsilon.
#   - Fixed significance thresholds and universe definitions.
#
# How to run:
#   Rscript --vanilla scripts/03-pathways/01_enrichment_paired_v2.R
#
# Assumptions:
#   - DE table contains unique gene_id rows from paired Tumor-vs-Normal DESeq2.
#   - log2FoldChange_shrunk, padj, and pvalue columns are present and meaningful.

#### 1) Load packages + seed ####
required_cran <- c("readr", "dplyr", "msigdbr")
required_bioc <- c(
  "fgsea",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "BiocParallel"
)

missing_cran <- required_cran[
  !vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)
]
missing_bioc <- required_bioc[
  !vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)
]

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
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(BiocParallel)
})

seed_env <- Sys.getenv("PIPELINE_SEED", "20260227")
seed_value <- suppressWarnings(as.integer(seed_env))
if (is.na(seed_value)) {
  stop(paste0("Invalid PIPELINE_SEED value: ", seed_env))
}
set.seed(seed_value)

de_results_path <- "results_v2/deseq2/deseq2_paired_v2_results.tsv"
go_out <- "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
go_representative_out <- "results_v2/enrichment/go_bp_ora_representative_v2.tsv"
go_tumor_out <- "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv"
go_normal_out <- "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv"
go_tumor_representative_out <- "results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv"
go_normal_representative_out <- "results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv"
gsea_out <- "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
diagnostics_out <- "results_v2/enrichment/enrichment_diagnostics_v2.tsv"
session_info_out <- "results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt"

#### 2) Load DE results + prepare gene ranks ####
de <- readr::read_tsv(de_results_path, show_col_types = FALSE)

# DE input contract:
#   gene_id                = unique gene identifier (Ensembl IDs expected)
#   pvalue                 = nominal p-value from DE test
#   padj                   = multiple-testing adjusted p-value (FDR)
#   log2FoldChange_shrunk  = stable effect size used for thresholding/sign direction
required_cols <- c("gene_id", "padj", "pvalue", "log2FoldChange_shrunk")
missing_cols <- setdiff(required_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required DE columns: ", paste(missing_cols, collapse = ", ")))
}

if (any(duplicated(de$gene_id))) {
  dup_ids <- unique(de$gene_id[duplicated(de$gene_id)])
  stop(
    paste0(
      "Duplicate gene_id values in DE results. Example IDs: ",
      paste(head(dup_ids, 10), collapse = ", ")
    )
  )
}

de <- de %>%
  mutate(lfc_for_selection = log2FoldChange_shrunk)

# Universe for ORA: all tested genes with non-NA padj.
universe_genes <- de %>%
  filter(!is.na(padj)) %>%
  pull(gene_id) %>%
  unique()

# Significant genes for ORA: padj < 0.05 and absolute shrunk effect >= 1.
sig_genes <- de %>%
  filter(
    !is.na(padj),
    !is.na(lfc_for_selection),
    padj < 0.05,
    abs(lfc_for_selection) >= 1
  ) %>%
  pull(gene_id) %>%
  unique()

tumor_sig_genes <- de %>%
  filter(
    !is.na(padj),
    !is.na(lfc_for_selection),
    padj < 0.05,
    lfc_for_selection >= 1
  ) %>%
  pull(gene_id) %>%
  unique()

normal_sig_genes <- de %>%
  filter(
    !is.na(padj),
    !is.na(lfc_for_selection),
    padj < 0.05,
    lfc_for_selection <= -1
  ) %>%
  pull(gene_id) %>%
  unique()

if (length(universe_genes) == 0) {
  stop("Universe is empty: no genes with non-NA padj.")
}
if (length(sig_genes) == 0) {
  warning("No significant genes meet padj < 0.05 and abs(LFC) >= 1; GO ORA will be empty.")
}

# GSEA ranking preparation:
#   Prefer DESeq2 Wald statistic when available;
#   otherwise use sign(LFC) * -log10(pvalue) to encode direction + strength.
de_rankable <- de %>%
  filter(!is.na(pvalue)) %>%
  mutate(
    rank_metric = if ("stat" %in% colnames(de) && !all(is.na(de$stat))) {
      stat
    } else {
      sign(ifelse(is.na(lfc_for_selection), 0, lfc_for_selection)) *
        -log10(pmax(pvalue, .Machine$double.xmin))
    }
  ) %>%
  filter(!is.na(rank_metric))

#### 3) Hallmark pathways retrieval + GSEA ####
# Hallmark GSEA summarizes coordinated pathway-level shifts:
#   NES > 0 indicates enrichment toward the Tumor-higher side of the ranking,
#   NES < 0 indicates enrichment toward the Normal-higher side.
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de_rankable$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, SYMBOL)

de_symbol_rank <- de_rankable %>%
  left_join(symbol_map, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup()

hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

hallmark_pathways <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

gsea_df <- data.frame()
gsea_method_used <- "none"
if (nrow(de_symbol_rank) > 0) {
  # Deterministic tie-breaking:
  #   1) stable pre-sort by descending rank_metric, then SYMBOL
  #   2) add tiny monotonic epsilon so exact ties resolve reproducibly
  # This avoids arbitrary ordering effects in preranked GSEA.
  de_symbol_rank <- de_symbol_rank %>%
    arrange(desc(rank_metric), SYMBOL)

  stats <- de_symbol_rank$rank_metric
  names(stats) <- de_symbol_rank$SYMBOL
  stats <- stats + seq_along(stats) * 1e-12
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]

  # SerialParam disables parallel worker variability and port issues, producing stable
  # behavior across constrained environments while keeping method fixed.
  gsea_method_used <- "fgseaMultilevel (SerialParam)"
  gsea_raw <- fgsea::fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = stats,
    minSize = 15,
    maxSize = 500,
    BPPARAM = BiocParallel::SerialParam()
  )

  # Output columns:
  #   pathway          = Hallmark pathway ID
  #   NES              = normalized enrichment score (direction + magnitude)
  #   pval / padj      = nominal and adjusted enrichment significance
  #   size             = number of genes from the pathway used in analysis
  #   ES               = raw enrichment score
  #   nMoreExtreme     = permutation-based tail count (if returned)
  #   leadingEdge_genes= semicolon-separated core driver genes (if returned)
  gsea_tbl <- as.data.frame(gsea_raw) %>%
    mutate(
      leadingEdge_genes = vapply(
        leadingEdge,
        function(x) paste(x, collapse = ";"),
        character(1)
      )
    )
  gsea_cols <- intersect(
    c(
      "pathway",
      "NES",
      "padj",
      "size",
      "pval",
      "ES",
      "nMoreExtreme",
      "leadingEdge_genes"
    ),
    colnames(gsea_tbl)
  )
  gsea_df <- gsea_tbl %>%
    dplyr::select(all_of(gsea_cols)) %>%
    arrange(padj)
}

#### 4) GO BP ORA ####
# The combined analysis is retained as a non-directional supplementary result.
# Directional analyses use the identical tested-gene universe, identifier rules,
# enrichment parameters and semantic-reduction method.
go_gene_lists_ensembl <- list(
  combined = sig_genes,
  tumor_higher = tumor_sig_genes,
  normal_higher = normal_sig_genes
)
go_keytype_used <- "ENSEMBL"
go_universe_input <- universe_genes
go_gene_lists_input <- go_gene_lists_ensembl

run_go <- function(genes, universe, keytype) {
  if (length(genes) == 0 || length(universe) == 0) return(NULL)
  clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = keytype,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

combined_try <- tryCatch(
  run_go(go_gene_lists_input$combined, go_universe_input, go_keytype_used),
  error = function(e) e
)

if (inherits(combined_try, "error")) {
  go_keytype_used <- "ENTREZID"
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(universe_genes),
    keytype = "ENSEMBL",
    columns = "ENTREZID"
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENSEMBL, ENTREZID)

  go_universe_input <- unique(map_df$ENTREZID)
  go_gene_lists_input <- lapply(go_gene_lists_ensembl, function(ids) {
    map_df %>%
      filter(ENSEMBL %in% ids) %>%
      pull(ENTREZID) %>%
      unique()
  })
  combined_try <- run_go(
    go_gene_lists_input$combined,
    go_universe_input,
    go_keytype_used
  )
}

go_results <- list(
  combined = combined_try,
  tumor_higher = run_go(
    go_gene_lists_input$tumor_higher,
    go_universe_input,
    go_keytype_used
  ),
  normal_higher = run_go(
    go_gene_lists_input$normal_higher,
    go_universe_input,
    go_keytype_used
  )
)

direction_labels <- c(
  combined = "Combined Tumor-higher and Normal-higher",
  tumor_higher = "Tumor-higher",
  normal_higher = "Normal-higher"
)

as_go_table <- function(result, analysis_name) {
  if (is.null(result)) return(data.frame())
  as.data.frame(result) %>%
    mutate(
      analysis_direction = unname(direction_labels[[analysis_name]]),
      selection_rule = ifelse(
        analysis_name == "tumor_higher",
        "padj < 0.05 and shrunken log2FC >= 1",
        ifelse(
          analysis_name == "normal_higher",
          "padj < 0.05 and shrunken log2FC <= -1",
          "padj < 0.05 and abs(shrunken log2FC) >= 1"
        )
      ),
      tested_universe_rule = "DESeq2 genes with non-NA padj",
      p_adjust_method = "Benjamini-Hochberg"
    )
}

simplify_go <- function(result, analysis_name) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) return(data.frame())
  significant_result <- result
  significant_result@result <- significant_result@result %>%
    filter(!is.na(p.adjust), p.adjust < 0.05)
  if (nrow(significant_result@result) == 0) return(data.frame())

  simplified <- clusterProfiler::simplify(
    significant_result,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )
  as.data.frame(simplified) %>%
    arrange(p.adjust) %>%
    slice_head(n = 30) %>%
    mutate(
      analysis_direction = unname(direction_labels[[analysis_name]]),
      simplification_method = "clusterProfiler::simplify",
      semantic_measure = "Wang",
      similarity_cutoff = 0.7,
      representative_selection = "lowest p.adjust within semantic group; first 30 by p.adjust"
    )
}

go_tables <- Map(as_go_table, go_results, names(go_results))
go_representative_tables <- Map(simplify_go, go_results, names(go_results))
go_df <- go_tables$combined
go_representative_df <- go_representative_tables$combined
go_tumor_df <- go_tables$tumor_higher
go_normal_df <- go_tables$normal_higher
go_tumor_representative_df <- go_representative_tables$tumor_higher
go_normal_representative_df <- go_representative_tables$normal_higher

go_sig_term_count <- sum(go_df$p.adjust < 0.05, na.rm = TRUE)

parse_ratio_denominator <- function(x) {
  if (length(x) == 0 || is.na(x[[1]]) || !grepl("/", x[[1]], fixed = TRUE)) {
    return(NA_integer_)
  }
  as.integer(strsplit(as.character(x[[1]]), "/", fixed = TRUE)[[1]][[2]])
}

go_annotated_count <- function(tbl, ratio_column) {
  if (nrow(tbl) == 0) return(NA_integer_)
  parse_ratio_denominator(tbl[[ratio_column]])
}

go_annotated_sig_n <- go_annotated_count(go_df, "GeneRatio")
go_annotated_universe_n <- go_annotated_count(go_df, "BgRatio")

diagnostics_df <- data.frame(
  metric = c(
    "gsea_rankable_ensembl_ids",
    "gsea_ensembl_ids_with_symbol",
    "gsea_unique_symbols_after_duplicate_resolution",
    "hallmark_gene_sets_tested",
    "hallmark_gene_sets_padj_lt_0.05",
    "ora_input_significant_genes",
    "ora_input_universe_genes",
    "ora_go_annotated_significant_genes",
    "ora_go_annotated_universe_genes",
    "go_terms_tested",
    "go_terms_padj_lt_0.05",
    "go_representative_terms_written",
    "ora_tumor_higher_input_genes",
    "ora_tumor_higher_mapped_genes",
    "ora_tumor_higher_terms_tested",
    "ora_tumor_higher_terms_padj_lt_0.05",
    "ora_tumor_higher_representative_terms_written",
    "ora_normal_higher_input_genes",
    "ora_normal_higher_mapped_genes",
    "ora_normal_higher_terms_tested",
    "ora_normal_higher_terms_padj_lt_0.05",
    "ora_normal_higher_representative_terms_written",
    "ora_keytype",
    "ora_multiple_testing_method",
    "ora_semantic_reduction"
  ),
  value = c(
    nrow(de_rankable),
    n_distinct(symbol_map$ENSEMBL),
    nrow(de_symbol_rank),
    nrow(gsea_df),
    sum(gsea_df$padj < 0.05, na.rm = TRUE),
    length(sig_genes),
    length(go_universe_input),
    go_annotated_sig_n,
    go_annotated_universe_n,
    nrow(go_df),
    go_sig_term_count,
    nrow(go_representative_df),
    length(tumor_sig_genes),
    go_annotated_count(go_tumor_df, "GeneRatio"),
    nrow(go_tumor_df),
    sum(go_tumor_df$p.adjust < 0.05, na.rm = TRUE),
    nrow(go_tumor_representative_df),
    length(normal_sig_genes),
    go_annotated_count(go_normal_df, "GeneRatio"),
    nrow(go_normal_df),
    sum(go_normal_df$p.adjust < 0.05, na.rm = TRUE),
    nrow(go_normal_representative_df),
    go_keytype_used,
    "Benjamini-Hochberg",
    "clusterProfiler::simplify; Wang similarity; cutoff 0.7; presentation only"
  ),
  stringsAsFactors = FALSE
)

#### 5) Write outputs + sessionInfo ####
dir.create(dirname(go_out), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(go_df, go_out)
readr::write_tsv(go_representative_df, go_representative_out)
readr::write_tsv(go_tumor_df, go_tumor_out)
readr::write_tsv(go_normal_df, go_normal_out)
readr::write_tsv(go_tumor_representative_df, go_tumor_representative_out)
readr::write_tsv(go_normal_representative_df, go_normal_representative_out)
readr::write_tsv(gsea_df, gsea_out)
readr::write_tsv(diagnostics_df, diagnostics_out)
capture.output(sessionInfo(), file = session_info_out)

cat("GO keyType used:", go_keytype_used, "\n")
cat("GO universe input size:", length(go_universe_input), "\n")
cat("Combined GO input/significant terms:", length(sig_genes), "/", go_sig_term_count, "\n")
cat("Tumor-higher GO input/significant terms:", length(tumor_sig_genes), "/", sum(go_tumor_df$p.adjust < 0.05, na.rm = TRUE), "\n")
cat("Normal-higher GO input/significant terms:", length(normal_sig_genes), "/", sum(go_normal_df$p.adjust < 0.05, na.rm = TRUE), "\n")
cat("GO representative terms written:", nrow(go_representative_df), "\n")
cat("Hallmark GSEA method used:", gsea_method_used, "\n")

go_top10 <- if (nrow(go_df) > 0) {
  go_df %>%
    arrange(p.adjust) %>%
    dplyr::select(Description, p.adjust) %>%
    slice_head(n = 10)
} else {
  data.frame()
}

gsea_top10 <- if (nrow(gsea_df) > 0) {
  gsea_df %>%
    arrange(padj) %>%
    dplyr::select(pathway, NES, padj) %>%
    slice_head(n = 10)
} else {
  data.frame()
}

nonempty_leading_edge <- FALSE
if (nrow(gsea_df) > 0) {
  nonempty_leading_edge <- any(
    !is.na(gsea_df$leadingEdge_genes) & nzchar(gsea_df$leadingEdge_genes)
  )
}

cat("GO ORA top 10 terms (Description + p.adjust):\n")
print(go_top10)
cat("Hallmark GSEA top 10 pathways by padj (pathway + NES + padj):\n")
print(gsea_top10)
cat("leadingEdge_genes_non_empty:", nonempty_leading_edge, "\n")

#!/usr/bin/env Rscript

# Purpose:
#   Run pathway-level interpretation for paired DE results using:
#   - Hallmark GSEA (program-level directional enrichment),
#   - GO Biological Process ORA (over-representation of significant genes).
#
# Inputs:
#   - results_v2/differential_expression/deseq2_paired_v2_results.tsv
#
# Outputs:
#   - results_v2/enrichment/hallmark_gsea_paired_v2.tsv
#   - results_v2/enrichment/go_bp_ora_paired_v2.tsv
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

de_results_path <- "results_v2/differential_expression/deseq2_paired_v2_results.tsv"
go_out <- "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
gsea_out <- "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
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
#   NES > 0 indicates Tumor-enriched pathways (given Tumor-vs-Normal contrast),
#   NES < 0 indicates Normal-enriched pathways.
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
# GO BP ORA asks which biological processes are over-represented among significant
# genes relative to the tested universe.
go_result <- NULL
go_keytype_used <- "ENSEMBL"
go_sig_input_n <- 0L
go_universe_n <- 0L

try_go_ensembl <- tryCatch(
  {
    if (length(sig_genes) == 0) {
      NULL
    } else {
      go_sig_input_n <- length(sig_genes)
      go_universe_n <- length(universe_genes)
      clusterProfiler::enrichGO(
        gene = sig_genes,
        universe = universe_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      )
    }
  },
  error = function(e) e
)

if (inherits(try_go_ensembl, "error")) {
  go_keytype_used <- "ENTREZID"
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(universe_genes),
    keytype = "ENSEMBL",
    columns = c("ENTREZID")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENSEMBL, ENTREZID)

  universe_entrez <- map_df %>%
    pull(ENTREZID) %>%
    unique()
  sig_entrez <- map_df %>%
    filter(ENSEMBL %in% sig_genes) %>%
    pull(ENTREZID) %>%
    unique()

  go_universe_n <- length(universe_entrez)
  go_sig_input_n <- length(sig_entrez)

  if (length(sig_entrez) > 0 && length(universe_entrez) > 0) {
    go_result <- clusterProfiler::enrichGO(
      gene = sig_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
  } else {
    go_result <- NULL
  }
} else {
  go_result <- try_go_ensembl
}

go_df <- if (!is.null(go_result)) as.data.frame(go_result) else data.frame()

go_sig_term_count <- 0L
if (nrow(go_df) > 0) {
  go_sig_term_count <- sum(go_df$p.adjust < 0.05, na.rm = TRUE)
}

#### 5) Write outputs + sessionInfo ####
dir.create(dirname(go_out), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(go_df, go_out)
readr::write_tsv(gsea_df, gsea_out)
capture.output(sessionInfo(), file = session_info_out)

cat("GO keyType used:", go_keytype_used, "\n")
cat("GO universe size:", go_universe_n, "\n")
cat("GO sig gene count used in ORA:", go_sig_input_n, "\n")
cat("GO enriched terms (p.adjust < 0.05):", go_sig_term_count, "\n")
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

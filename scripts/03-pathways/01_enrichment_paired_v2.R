#!/usr/bin/env Rscript

required_cran <- c("readr", "dplyr", "msigdbr")
required_bioc <- c("fgsea", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]

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
})

input_path <- "results_v2/differential_expression/deseq2_paired_v2_results.tsv"
go_out <- "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
gsea_out <- "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
session_info_out <- "results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt"

de <- readr::read_tsv(input_path, show_col_types = FALSE)
required_cols <- c("gene_id", "padj", "pvalue")
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

lfc_col <- "log2FoldChange"
if ("log2FoldChange_shrunk" %in% colnames(de) && !all(is.na(de$log2FoldChange_shrunk))) {
  lfc_col <- "log2FoldChange_shrunk"
}

de <- de %>%
  mutate(lfc_for_selection = .data[[lfc_col]])

# Task A: define tested universe and significant genes.
universe_genes <- de %>%
  filter(!is.na(padj)) %>%
  pull(gene_id) %>%
  unique()

sig_genes <- de %>%
  filter(!is.na(padj), !is.na(lfc_for_selection), padj < 0.05, abs(lfc_for_selection) >= 1) %>%
  pull(gene_id) %>%
  unique()

if (length(universe_genes) == 0) {
  stop("Universe is empty: no genes with non-NA padj.")
}
if (length(sig_genes) == 0) {
  warning("No significant genes meet padj < 0.05 and abs(LFC) >= 1; GO ORA will be empty.")
}

# Task B: GO ORA with explicit universe.
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

dir.create(dirname(go_out), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(go_df, go_out)

# Task C: Hallmark GSEA with leading edge preserved.
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
if (nrow(de_symbol_rank) > 0) {
  stats <- de_symbol_rank$rank_metric
  names(stats) <- de_symbol_rank$SYMBOL
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]

  gsea_raw <- fgsea::fgsea(
    pathways = hallmark_pathways,
    stats = stats,
    minSize = 15,
    maxSize = 500,
    nperm = 10000
  )

  gsea_df <- as.data.frame(gsea_raw) %>%
    mutate(
      leadingEdge_genes = vapply(
        leadingEdge,
        function(x) paste(x, collapse = ";"),
        character(1)
      )
    ) %>%
    dplyr::select(pathway, NES, padj, size, pval, ES, nMoreExtreme, leadingEdge_genes) %>%
    arrange(padj)
}

readr::write_tsv(gsea_df, gsea_out)

capture.output(sessionInfo(), file = session_info_out)

cat("GO keyType used:", go_keytype_used, "\n")
cat("GO universe size:", go_universe_n, "\n")
cat("GO sig gene count used in ORA:", go_sig_input_n, "\n")
cat("GO enriched terms (p.adjust < 0.05):", go_sig_term_count, "\n")

go_top10 <- if (nrow(go_df) > 0) {
  go_df %>% arrange(p.adjust) %>% dplyr::select(Description, p.adjust) %>% slice_head(n = 10)
} else {
  data.frame()
}

gsea_top10 <- if (nrow(gsea_df) > 0) {
  gsea_df %>% arrange(padj) %>% dplyr::select(pathway, NES, padj) %>% slice_head(n = 10)
} else {
  data.frame()
}

nonempty_leading_edge <- FALSE
if (nrow(gsea_df) > 0) {
  nonempty_leading_edge <- any(!is.na(gsea_df$leadingEdge_genes) & nzchar(gsea_df$leadingEdge_genes))
}

cat("GO ORA top 10 terms (Description + p.adjust):\n")
print(go_top10)
cat("Hallmark GSEA top 10 pathways by padj (pathway + NES + padj):\n")
print(gsea_top10)
cat("leadingEdge_genes_non_empty:", nonempty_leading_edge, "\n")

#!/usr/bin/env bash
set -euo pipefail

report_dir="report"
scripts_dir="scripts"
mkdir -p "$report_dir" "$scripts_dir"

log_path="$report_dir/verification_log.txt"
if [ -f "$log_path" ]; then
  ts=$(date +"%Y%m%d_%H%M%S")
  mv "$log_path" "$report_dir/verification_log_backup_${ts}.txt"
fi

get_col_index() {
  local header="$1"
  local target="$2"
  local IFS=$'\t'
  local -a cols
  read -r -a cols <<< "$header"
  local i=1
  for col in "${cols[@]}"; do
    if [ "$col" = "$target" ]; then
      echo "$i"
      return 0
    fi
    i=$((i + 1))
  done
  echo ""
  return 1
}

check_file() {
  local path="$1"
  if [ -f "$path" ]; then
    local size
    size=$(stat -f%z "$path")
    echo "OK\t$path\t${size} bytes"
  else
    echo "MISSING\t$path"
  fi
}

{
  echo "Dataset: breast_cancer_rnaseq (GSE306117)"
  echo "Date: $(date)"
  echo "System: $(uname -a)"
  echo ""
  echo "Sample verification (metadata)"

  metadata_path="data/metadata/metadata.tsv"
  if [ -f "$metadata_path" ]; then
    echo "use_main_tumor_vs_normal (column 10) counts:"
    tail -n +2 "$metadata_path" | tr -d '\r' | awk -F '\t' '{print $10}' | sort | uniq -c | awk '{printf "  %s\t%s\n", $2, $1}'

    yes_count=$(tail -n +2 "$metadata_path" | tr -d '\r' | awk -F '\t' '$10=="yes"{c++} END{print c+0}')
    echo "rows with use_main_tumor_vs_normal == yes: $yes_count"

    echo "condition counts among use_main_tumor_vs_normal == yes:"
    tail -n +2 "$metadata_path" | tr -d '\r' | awk -F '\t' '$10=="yes"{counts[$5]++} END{for (k in counts) printf "  %s\t%d\n", k, counts[k]}' | sort
  else
    echo "metadata file missing: $metadata_path"
  fi

  echo ""
  echo "DESeq2 verification"
  de_path="results/differential_expression/deseq2_results.tsv"
  if [ -f "$de_path" ]; then
    gene_count=$(tail -n +2 "$de_path" | tr -d '\r' | wc -l | awk '{print $1}')
    echo "genes tested: $gene_count"

    header=$(head -n 1 "$de_path" | tr -d '\r')
    padj_idx=$(get_col_index "$header" "padj" || true)
    lfc_idx=$(get_col_index "$header" "log2FoldChange" || true)
    if [ -z "$padj_idx" ] || [ -z "$lfc_idx" ]; then
      echo "missing required columns in DESeq2 results (padj/log2FoldChange)"
    else
      sig_counts=$(tail -n +2 "$de_path" | tr -d '\r' | awk -F '\t' -v padj="$padj_idx" -v lfc="$lfc_idx" '
        ($padj != "" && $padj != "NA" && $padj < 0.05){
          sig++
          if ($lfc > 0) up++
          else if ($lfc < 0) down++
        }
        END{
          printf "%d\t%d\t%d\n", sig+0, up+0, down+0
        }')
      sig=$(echo "$sig_counts" | awk -F '\t' '{print $1}')
      up=$(echo "$sig_counts" | awk -F '\t' '{print $2}')
      down=$(echo "$sig_counts" | awk -F '\t' '{print $3}')
      echo "significant genes (padj < 0.05): $sig"
      echo "upregulated in Tumor: $up"
      echo "downregulated in Tumor: $down"
    fi
  else
    echo "DESeq2 results file missing: $de_path"
  fi

  echo ""
  echo "Enrichment verification"
  hallmark_path="results/enrichment/hallmark_gsea.tsv"
  if [ -f "$hallmark_path" ]; then
    header=$(head -n 1 "$hallmark_path" | tr -d '\r')
    padj_idx=$(get_col_index "$header" "padj" || true)
    if [ -z "$padj_idx" ]; then
      echo "Hallmark GSEA: missing padj column"
    else
      sig_count=$(tail -n +2 "$hallmark_path" | tr -d '\r' | awk -F '\t' -v padj="$padj_idx" '$padj != "" && $padj != "NA" && $padj < 0.05 {c++} END{print c+0}')
      echo "Hallmark significant pathways (padj < 0.05): $sig_count"
    fi
  else
    echo "Hallmark GSEA file missing: $hallmark_path"
  fi

  go_path="results/enrichment/go_enrich.tsv"
  if [ -f "$go_path" ]; then
    header=$(head -n 1 "$go_path" | tr -d '\r')
    padj_idx=$(get_col_index "$header" "p.adjust" || true)
    if [ -z "$padj_idx" ]; then
      echo "GO enrichment: missing p.adjust column"
    else
      sig_count=$(tail -n +2 "$go_path" | tr -d '\r' | awk -F '\t' -v padj="$padj_idx" '$padj != "" && $padj != "NA" && $padj < 0.05 {c++} END{print c+0}')
      echo "GO BP significant terms (p.adjust < 0.05): $sig_count"
    fi
  else
    echo "GO enrichment file missing: $go_path"
  fi

  echo ""
  echo "Figure existence check"
  check_file "figures/qc/library_size.png"
  check_file "figures/qc/pca.png"
  check_file "figures/de/ma_plot.png"
  check_file "figures/de/volcano.png"
  check_file "figures/de/heatmap_top50.png"
  check_file "figures/enrichment/hallmark_gsea_top10.png"
  check_file "figures/enrichment/go_barplot_top15.png"
} > "$log_path"

echo "Wrote: $log_path"

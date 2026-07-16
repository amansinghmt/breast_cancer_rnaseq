#!/usr/bin/env bash
# Smoke test for the maintained v2 paired RNA-seq pipeline.
# Checks that scripts/run_v2.sh runs end-to-end and required outputs exist.
# Run with: bash tests/smoke_test_v2.sh
# Expected success message: SMOKE TEST PASSED

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"


die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

step() {
  printf '\n[%s/7] %s\n' "$1" "$2"
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || die "Missing required file: $path"
}

require_nonempty_file() {
  local path="$1"
  require_file "$path"
  [[ -s "$path" ]] || die "Required file is empty: $path"
}

require_output_manifest_entry() {
  local entry="$1"
  grep -Eq "^${entry}[[:space:]]" results_v2/output_manifest.tsv || \
    die "Missing output manifest entry: ${entry}"
}

latest_run_log() {
  local latest
  latest="$(find results_v2/logs -maxdepth 1 -type f -name 'run_v2_*.log' -print 2>/dev/null | sort | tail -n 1)"
  [[ -n "$latest" ]] || die "No run log found matching results_v2/logs/run_v2_*.log"
  printf '%s\n' "$latest"
}

require_exact_figure_set() {
  local figure_dir="figures_v2/final"
  local expected actual

  [[ -d "$figure_dir" ]] || die "Missing figure directory: $figure_dir"

  expected="$(printf '%s\n' F01.png F02.png F03.png F04.png F05.png F06.png F07.png)"
  actual="$(find "$figure_dir" -mindepth 1 -maxdepth 1 -exec basename {} \; | sort)"

  if [[ "$actual" != "$expected" ]]; then
    printf 'Expected figures:\n%s\n\nActual figures:\n%s\n' "$expected" "$actual" >&2
    die "Figure contract mismatch in $figure_dir"
  fi

  for figure in F01.png F02.png F03.png F04.png F05.png F06.png F07.png; do
    require_nonempty_file "$figure_dir/$figure"
  done
}

step 1 "Running maintained v2 pipeline"
bash scripts/run_v2.sh || die "Pipeline execution failed: bash scripts/run_v2.sh"

step 2 "Checking canonical DE outputs"
require_nonempty_file results_v2/deseq2/deseq2_paired_v2_results.tsv
require_nonempty_file results_v2/deseq2/deseq2_paired_v2_samples_used.tsv
require_nonempty_file results_v2/deseq2/deseq2_paired_v2_vst.tsv
require_nonempty_file results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv
require_nonempty_file results_v2/deseq2/deseq2_paired_v2_size_factors.tsv
require_nonempty_file results_v2/deseq2/sessionInfo_paired_v2.txt

step 3 "Checking enrichment outputs"
require_nonempty_file results_v2/enrichment/hallmark_gsea_paired_v2.tsv
require_nonempty_file results_v2/enrichment/go_bp_ora_paired_v2.tsv
require_nonempty_file results_v2/enrichment/go_bp_ora_representative_v2.tsv
require_nonempty_file results_v2/enrichment/enrichment_diagnostics_v2.tsv
require_nonempty_file results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt

for robustness_file in \
  analysis_metrics_v2.tsv \
  de_threshold_sensitivity_v2.tsv \
  low_count_prefilter_sensitivity_v2.tsv \
  lfc_agreement_v2.tsv \
  pca_outlier_summary_v2.tsv \
  top_de_genes_v2.tsv; do
  require_nonempty_file "results_v2/robustness/${robustness_file}"
done

step 4 "Checking metadata/session/manifests"
require_nonempty_file results_v2/metadata/paired_manifest.tsv
require_nonempty_file results_v2/fig_manifest.tsv
require_nonempty_file results_v2/output_manifest.tsv
require_nonempty_file results_v2/sessionInfo.txt

step 5 "Checking final figure contract"
require_exact_figure_set
for vector_figure in figures_v2/vector/F0*.pdf; do
  require_nonempty_file "$vector_figure"
done

step 6 "Checking run logs"
latest_log="$(latest_run_log)"
require_nonempty_file "$latest_log"
printf 'Latest run log: %s\n' "$latest_log"

step 7 "Checking manifest consistency"
for entry in de_results de_vst de_diagnostics enrichment_hallmark enrichment_go_bp enrichment_go_bp_representative robustness_metrics figure_F01 figure_F02 figure_F03 figure_F04 figure_F05 figure_F06 figure_F07 figure_F01_pdf figure_F07_pdf; do
  require_output_manifest_entry "$entry"
done

fig_rows="$(awk 'NR>1 && NF>0 {count++} END {print count+0}' results_v2/fig_manifest.tsv)"
[[ "$fig_rows" -eq 7 ]] || die "Expected 7 data rows in results_v2/fig_manifest.tsv, found $fig_rows"

printf '\nSMOKE TEST PASSED\n'

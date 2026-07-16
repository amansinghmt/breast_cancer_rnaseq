#!/usr/bin/env bash
# Purpose:
#   Regenerate all v2 analysis outputs from scratch (metadata, paired DE, enrichment,
#   figures, manifests, and logs) with strict validation.
#
# Inputs:
#   - data/raw/GSE306117_series_matrix.txt
#   - data/raw/GSE306117_RAW/
#   - renv.lock
#   - Existing pipeline scripts under scripts/00-metadata, 01-qc, 02-de, 03-pathways,
#     and 04-figures
#
# Outputs:
#   - results_v2/metadata/paired_manifest.tsv
#   - results_v2/deseq2/*
#   - results_v2/enrichment/*
#   - figures_v2/final/F01.png ... F07.png
#   - results_v2/fig_manifest.tsv
#   - results_v2/output_manifest.tsv
#   - results_v2/sessionInfo.txt
#   - results_v2/logs/run_v2_YYYYmmdd_HHMMSS.log
#
# Determinism / reproducibility:
#   - Strict shell mode and explicit step ordering.
#   - A single pipeline seed (PIPELINE_SEED, default 20260227) exported to R.
#   - R scripts are executed with --vanilla.
#   - Final figure contract enforces exactly F01..F07 in figures_v2/final.
#
# How to run:
#   bash scripts/run_v2.sh
#
# Assumptions:
#   - The paired cohort contains one Tumor and one Normal sample per patient
#     (21 matched pairs expected from this dataset).
#   - Required Python and R dependencies are available (renv restore + Python libs).
set -eEuo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

# PIPELINE_SEED is consumed by downstream R scripts to keep seeded operations stable.
RUN_TS="$(date '+%Y%m%d_%H%M%S')"
PIPELINE_SEED="${PIPELINE_SEED:-20260227}"
export PIPELINE_SEED

mkdir -p results_v2
# Temporary writable Matplotlib/font caches prevent permissions issues on macOS.
# These locations are used only for Python subprocesses and are removed on exit.
PYTHON_MPLCONFIGDIR="$(mktemp -d "${TMPDIR:-/tmp}/run_v2_mplconfig.XXXXXX")"
PYTHON_CACHE_HOME="$(mktemp -d "${TMPDIR:-/tmp}/run_v2_cache.XXXXXX")"
mkdir -p "${PYTHON_MPLCONFIGDIR}" "${PYTHON_CACHE_HOME}/fontconfig"

CURRENT_STEP="initialization"

# Print a fatal error and exit.
die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Print a timestamped step banner and track current step for error reporting.
step() {
  CURRENT_STEP="$1"
  printf "\n[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "${CURRENT_STEP}"
}

# Trap handler that reports where and why the run failed.
on_error() {
  local line_no="$1"
  local last_cmd="$2"
  echo
  echo "FAILED at step ${CURRENT_STEP}"
  echo "Line: ${line_no}"
  echo "Last command: ${last_cmd}"
  exit 1
}

trap 'on_error "${LINENO}" "${BASH_COMMAND}"' ERR

# Run an R script with --vanilla and an explicit seed propagated from PIPELINE_SEED.
run_r_script() {
  local script_path="$1"
  [[ -f "${script_path}" ]] || die "Missing R script: ${script_path}"
  Rscript --vanilla -e '
    seed_env <- Sys.getenv("PIPELINE_SEED", "20260227")
    seed_val <- suppressWarnings(as.integer(seed_env))
    if (is.na(seed_val)) {
      stop(paste0("Invalid PIPELINE_SEED value: ", seed_env))
    }
    set.seed(seed_val)
    source(commandArgs(trailingOnly = TRUE)[1], chdir = FALSE)
  ' "${script_path}"
}

# Run Python with isolated cache/config directories for deterministic, writable behavior.
run_python() {
  XDG_CACHE_HOME="${PYTHON_CACHE_HOME}" \
  MPLCONFIGDIR="${PYTHON_MPLCONFIGDIR}" \
  MPLBACKEND=Agg \
  "${PYTHON_BIN}" "$@"
}

# Map canonical figure IDs (F01..F07) to source filenames emitted by figure scripts.
# This keeps figure generation logic unchanged while enforcing the final naming contract.
figure_source_name() {
  local figure_id="$1"
  case "${figure_id}" in
    F01) echo "F01_qc_library_size_pairs.png" ;;
    F02) echo "F02_qc_pca_pairs.png" ;;
    F03) echo "F03_de_ma_paired_v2.png" ;;
    F04) echo "F04_de_volcano_paired_v2.png" ;;
    F05) echo "F05_de_heatmap_top40_paired_v2.png" ;;
    F06) echo "F06_bio_hallmark_nes_paired_v2.png" ;;
    F07) echo "F07_bio_go_bp_dotplot_paired_v2.png" ;;
    *) return 1 ;;
  esac
}

print_help() {
  cat <<'EOF'
Usage:
  bash scripts/run_v2.sh
  bash scripts/run_v2.sh -h|--help

Description:
  Regenerates all v2 pipeline outputs from scratch with strict validation:
    A) metadata + paired manifest
    B) paired DESeq2 differential expression
    C) enrichment (Hallmark GSEA + GO BP ORA)
    D) figures F01..F07
    E) manifests + final output checks

Environment:
  PIPELINE_SEED   Seed propagated to R scripts (default: 20260227)
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  print_help
  exit 0
fi

if [[ $# -gt 0 ]]; then
  die "Unknown argument(s): $* (use -h or --help)"
fi

mkdir -p results_v2/logs
LOG_FILE="results_v2/logs/run_v2_${RUN_TS}.log"
touch "${LOG_FILE}"

LOG_PIPE="results_v2/logs/run_v2_${RUN_TS}.pipe"
rm -f "${LOG_PIPE}"
mkfifo "${LOG_PIPE}"
tee -a "${LOG_FILE}" < "${LOG_PIPE}" &
TEE_PID=$!
exec > "${LOG_PIPE}" 2>&1

# Ensure log pipe + temporary Python cache dirs are cleaned up on exit.
cleanup_logging() {
  rm -f "${LOG_PIPE}" || true
  rm -rf "${PYTHON_MPLCONFIGDIR}" "${PYTHON_CACHE_HOME}" || true
  if [[ -n "${TEE_PID:-}" ]]; then
    kill "${TEE_PID}" 2>/dev/null || true
  fi
}
trap cleanup_logging EXIT

# Pipeline overview:
#   Step A: Build metadata/QC and derive paired manifest for analysis cohort.
#   Step B: Run paired DESeq2 model (~ patient_id + condition_main) and validate outputs.
#   Step C: Run enrichment and robustness summaries; confirm required tables.
#   Step D: Render all figures, then enforce exactly F01..F07 in final output folder.
#   Step E: Validate final artifacts and write figure/output manifests with checksums.

step "Preflight checks"
command -v Rscript >/dev/null 2>&1 || die "Rscript not found. Install R and retry."
if [[ -x "${REPO_ROOT}/.venv/bin/python3" ]]; then
  PYTHON_BIN="${REPO_ROOT}/.venv/bin/python3"
else
  PYTHON_BIN="$(command -v python3 || true)"
fi
[[ -n "${PYTHON_BIN}" ]] || die "python3 not found. Install Python 3 and retry."
[[ -f "renv.lock" ]] || die "renv.lock is missing at repo root."
[[ -f "data/raw/GSE306117_series_matrix.txt" ]] || die "Missing required input: data/raw/GSE306117_series_matrix.txt"
[[ -d "data/raw/GSE306117_RAW" ]] || die "Missing required input directory: data/raw/GSE306117_RAW"

run_python - <<'PY'
import importlib.util
required = ["pandas", "numpy", "matplotlib", "sklearn"]
missing = [m for m in required if importlib.util.find_spec(m) is None]
if missing:
    raise SystemExit("Missing required Python packages: " + ", ".join(missing))
PY

step "Initialize output directories"
mkdir -p results_v2 results_v2/logs results_v2/metadata results_v2/deseq2 results_v2/enrichment
mkdir -p results_v2/robustness figures_v2/final figures_v2/archive figures_v2/vector

ARCHIVE_RUN_DIR="figures_v2/archive/run_${RUN_TS}"
mkdir -p "${ARCHIVE_RUN_DIR}"

shopt -s dotglob nullglob
existing_final=(figures_v2/final/*)
if ((${#existing_final[@]} > 0)); then
  mv "${existing_final[@]}" "${ARCHIVE_RUN_DIR}/"
fi
shopt -u dotglob nullglob

if [[ -d "figures_v2/vector" ]]; then
  mkdir -p "${ARCHIVE_RUN_DIR}/vector"
  shopt -s nullglob
  existing_vector=(figures_v2/vector/*)
  if ((${#existing_vector[@]} > 0)); then
    mv "${existing_vector[@]}" "${ARCHIVE_RUN_DIR}/vector/"
  fi
  shopt -u nullglob
fi

rm -rf results_v2/metadata
rm -rf results_v2/deseq2
rm -rf results_v2/enrichment
rm -rf results_v2/robustness
rm -rf results_v2/differential_expression
rm -rf results_v2/qc
rm -rf results_v2/figures
rm -rf figures_v2/de
rm -rf figures_v2/vector
rm -f results_v2/sessionInfo.txt
rm -f results_v2/fig_manifest.tsv
rm -f results_v2/output_manifest.tsv

mkdir -p results_v2/metadata results_v2/deseq2 results_v2/enrichment results_v2/robustness
mkdir -p figures_v2/final figures_v2/vector figures_v2/archive "${ARCHIVE_RUN_DIR}"

step "Restore R environment with renv"
Rscript --vanilla -e 'if(!requireNamespace("renv", quietly=TRUE)) { cat("ERROR: Package \"renv\" is not installed. Install with install.packages(\"renv\").\n", file=stderr()); quit(status=1) }'
Rscript --vanilla -e 'if(!requireNamespace("renv", quietly=TRUE)) quit(status=1); renv::restore(prompt=FALSE)'

Rscript --vanilla -e '
  pkgs <- c(
    "DESeq2", "readr", "dplyr", "tidyr", "ggplot2", "scales",
    "msigdbr", "fgsea", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "BiocParallel",
    "patchwork", "ggrepel"
  )
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat(
      "Missing required R packages after renv::restore: ",
      paste(missing, collapse = ", "),
      "\n",
      sep = "",
      file = stderr()
    )
    quit(status = 1)
  }
'

Rscript --vanilla -e 'sessionInfo()' > results_v2/sessionInfo.txt
[[ -s "results_v2/sessionInfo.txt" ]] || die "results_v2/sessionInfo.txt is missing or empty."

step "Step A: metadata and paired manifest"
run_python scripts/00-metadata/01_make_metadata.py
run_python scripts/01-qc/02_merge_htseq_counts.py
run_python scripts/01-qc/03_qc_plots.py
run_python scripts/00-metadata/02_make_sample_manifest_v2.py

[[ -s "data/metadata/sample_manifest.tsv" ]] || die "data/metadata/sample_manifest.tsv is missing or empty."
cp -f data/metadata/sample_manifest.tsv results_v2/metadata/paired_manifest.tsv

# Validation block: paired manifest must exist, contain required columns, and encode
# exactly one Tumor + one Normal sample per included patient.
run_python - "${REPO_ROOT}" <<'PY'
import csv
import os
import sys

root = sys.argv[1]
manifest_path = os.path.join(root, "results_v2", "metadata", "paired_manifest.tsv")

if not os.path.exists(manifest_path) or os.path.getsize(manifest_path) == 0:
    raise SystemExit(f"Missing or empty paired manifest: {manifest_path}")

with open(manifest_path, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    cols = reader.fieldnames or []
    rows = list(reader)

required_cols = {"sample_id", "patient_id", "condition_main", "include_paired"}
missing_cols = sorted(required_cols - set(cols))
if missing_cols:
    raise SystemExit(
        "paired_manifest.tsv missing required columns: " + ", ".join(missing_cols)
    )
if not rows:
    raise SystemExit("paired_manifest.tsv has no data rows.")

paired_rows = [
    r for r in rows if str(r.get("include_paired", "")).strip().upper() == "TRUE"
]
if not paired_rows:
    raise SystemExit("paired_manifest.tsv has zero include_paired == TRUE rows.")

by_patient = {}
for row in paired_rows:
    pid = str(row["patient_id"]).strip()
    by_patient.setdefault(pid, []).append(str(row["condition_main"]).strip())

bad = []
for pid, conds in sorted(by_patient.items()):
    n_tumor = sum(c == "Tumor" for c in conds)
    n_normal = sum(c == "Normal" for c in conds)
    if len(conds) != 2 or n_tumor != 1 or n_normal != 1:
        bad.append(f"{pid}(n={len(conds)},tumor={n_tumor},normal={n_normal})")

if bad:
    raise SystemExit("Invalid paired structure in manifest: " + "; ".join(bad))
PY

step "Step B: DESeq2 paired differential expression"
run_r_script scripts/02-de/01_deseq2_paired_v2.R

# Validation block: DE table must be non-empty with required columns, and samples_used
# must match the paired manifest cohort exactly (IDs + patient/condition labels).
run_python - "${REPO_ROOT}" <<'PY'
import csv
import os
import sys

root = sys.argv[1]
de_path = os.path.join(root, "results_v2", "deseq2", "deseq2_paired_v2_results.tsv")
samples_used_path = os.path.join(root, "results_v2", "deseq2", "deseq2_paired_v2_samples_used.tsv")
manifest_path = os.path.join(root, "results_v2", "metadata", "paired_manifest.tsv")

for path in (de_path, samples_used_path, manifest_path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"Missing or empty file: {path}")

with open(de_path, newline="") as fh:
    de_reader = csv.DictReader(fh, delimiter="\t")
    de_cols = de_reader.fieldnames or []
    de_rows = list(de_reader)

required_de_cols = {
    "gene_id", "log2FoldChange", "pvalue", "padj", "log2FoldChange_shrunk"
}
missing_de_cols = sorted(required_de_cols - set(de_cols))
if missing_de_cols:
    raise SystemExit(
        "deseq2_paired_v2_results.tsv missing required columns: "
        + ", ".join(missing_de_cols)
    )
if not de_rows:
    raise SystemExit("deseq2_paired_v2_results.tsv has no data rows.")

with open(manifest_path, newline="") as fh:
    manifest_reader = csv.DictReader(fh, delimiter="\t")
    manifest_rows = list(manifest_reader)

manifest_paired = {}
for row in manifest_rows:
    if str(row.get("include_paired", "")).strip().upper() != "TRUE":
        continue
    sid = str(row["sample_id"]).strip()
    manifest_paired[sid] = (
        str(row["patient_id"]).strip(),
        str(row["condition_main"]).strip(),
    )

if not manifest_paired:
    raise SystemExit("No include_paired == TRUE rows found in paired_manifest.tsv.")

with open(samples_used_path, newline="") as fh:
    used_reader = csv.DictReader(fh, delimiter="\t")
    used_cols = set(used_reader.fieldnames or [])
    required_used_cols = {"sample_id", "patient_id", "condition_main"}
    missing_used_cols = sorted(required_used_cols - used_cols)
    if missing_used_cols:
        raise SystemExit(
            "deseq2_paired_v2_samples_used.tsv missing required columns: "
            + ", ".join(missing_used_cols)
        )
    used_rows = list(used_reader)

used_map = {}
for row in used_rows:
    sid = str(row["sample_id"]).strip()
    used_map[sid] = (
        str(row["patient_id"]).strip(),
        str(row["condition_main"]).strip(),
    )

missing_in_used = sorted(set(manifest_paired) - set(used_map))
extra_in_used = sorted(set(used_map) - set(manifest_paired))
if missing_in_used or extra_in_used:
    parts = []
    if missing_in_used:
        parts.append("missing in samples_used: " + ", ".join(missing_in_used))
    if extra_in_used:
        parts.append("extra in samples_used: " + ", ".join(extra_in_used))
    raise SystemExit("Manifest/sample mismatch: " + "; ".join(parts))

mismatch = []
for sid in sorted(manifest_paired):
    if manifest_paired[sid] != used_map[sid]:
        mismatch.append(sid)

if mismatch:
    raise SystemExit(
        "Manifest/sample mismatch on patient_id/condition_main for sample_id(s): "
        + ", ".join(mismatch)
    )
PY

step "Step C: enrichment (GSEA + ORA)"
run_r_script scripts/03-pathways/01_enrichment_paired_v2.R

# Validation block: enrichment outputs must exist and at least one table row must be
# produced across Hallmark GSEA and GO BP ORA results.
run_python - "${REPO_ROOT}" <<'PY'
import csv
import os
import sys

root = sys.argv[1]
hallmark_path = os.path.join(root, "results_v2", "enrichment", "hallmark_gsea_paired_v2.tsv")
go_path = os.path.join(root, "results_v2", "enrichment", "go_bp_ora_paired_v2.tsv")
go_representative_path = os.path.join(
    root, "results_v2", "enrichment", "go_bp_ora_representative_v2.tsv"
)
go_tumor_path = os.path.join(
    root, "results_v2", "enrichment", "go_bp_ora_tumor_higher_paired_v2.tsv"
)
go_normal_path = os.path.join(
    root, "results_v2", "enrichment", "go_bp_ora_normal_higher_paired_v2.tsv"
)
go_tumor_rep_path = os.path.join(
    root, "results_v2", "enrichment", "go_bp_ora_tumor_higher_representative_v2.tsv"
)
go_normal_rep_path = os.path.join(
    root, "results_v2", "enrichment", "go_bp_ora_normal_higher_representative_v2.tsv"
)
diagnostics_path = os.path.join(
    root, "results_v2", "enrichment", "enrichment_diagnostics_v2.tsv"
)

for path in (
    hallmark_path, go_path, go_representative_path, go_tumor_path,
    go_normal_path, go_tumor_rep_path, go_normal_rep_path, diagnostics_path
):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"Missing or empty enrichment output: {path}")

data_rows = 0
for path in (hallmark_path, go_path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    data_rows += len(rows)

if data_rows == 0:
    raise SystemExit("All enrichment result tables are empty.")
PY

step "Step C2: robustness and presentation summaries"
run_r_script scripts/05-reporting/01_build_robustness_summaries.R

for summary_path in \
  results_v2/robustness/analysis_metrics_v2.tsv \
  results_v2/robustness/de_threshold_sensitivity_v2.tsv \
  results_v2/robustness/low_count_prefilter_sensitivity_v2.tsv \
  results_v2/robustness/lfc_agreement_v2.tsv \
  results_v2/robustness/pca_outlier_summary_v2.tsv \
  results_v2/robustness/top_de_genes_v2.tsv; do
  [[ -s "${summary_path}" ]] || die "Missing robustness output: ${summary_path}"
done

step "Step D: figures F01-F07"
run_r_script scripts/04-figures/01_publication_figures_paired_v2.R
run_r_script scripts/04-figures/02_pca_paired_v2.R
run_r_script scripts/04-figures/03_ma_paired_v2.R
run_r_script scripts/04-figures/04_volcano_paired_v2.R
run_r_script scripts/04-figures/05_heatmap_paired_v2.R
run_r_script scripts/04-figures/06_hallmark_bar_paired_v2.R
run_r_script scripts/04-figures/07_go_bp_dotplot_paired_v2.R

# Figure scripts emit descriptive source names. Copy those exact files to canonical
# F01..F07 names required by downstream manifests and final contract checks.
for figure_id in F01 F02 F03 F04 F05 F06 F07; do
  src_name="$(figure_source_name "${figure_id}")"
  src_path="figures_v2/final/${src_name}"
  dst_path="figures_v2/final/${figure_id}.png"
  [[ -s "${src_path}" ]] || die "Missing expected figure source file: ${src_path}"
  cp -f "${src_path}" "${dst_path}"
done

EXTRA_FIG_ARCHIVE_DIR="${ARCHIVE_RUN_DIR}/final_extra_files"
mkdir -p "${EXTRA_FIG_ARCHIVE_DIR}"

# Archive any non-canonical figure files so figures_v2/final contains only F01..F07.
shopt -s dotglob nullglob
final_items=(figures_v2/final/*)
for item in "${final_items[@]}"; do
  base_name="$(basename "${item}")"
  case "${base_name}" in
    F01.png|F02.png|F03.png|F04.png|F05.png|F06.png|F07.png) ;;
    *) mv "${item}" "${EXTRA_FIG_ARCHIVE_DIR}/${base_name}" ;;
  esac
done
shopt -u dotglob nullglob

step "Step E: strict final validation + manifests"
# Validation block:
#   - Enforce exact figure folder contract (only F01..F07, all non-empty)
#   - Build results_v2/fig_manifest.tsv (figure_id, path, bytes, md5)
#   - Build results_v2/output_manifest.tsv including run_log entry
#   - Hard-fail on any missing/empty required output
run_python - "${REPO_ROOT}" "${LOG_FILE}" <<'PY'
import csv
import hashlib
import os
import sys

root = sys.argv[1]
log_rel = sys.argv[2]

def abs_path(rel_path):
    return os.path.join(root, rel_path)

def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def ensure_nonempty(path):
    if not os.path.exists(path):
        return "missing"
    if os.path.getsize(path) == 0:
        return "empty"
    return ""

expected_figure_files = [f"F{i:02d}.png" for i in range(1, 8)]
expected_figure_ids = [f"F{i:02d}" for i in range(1, 8)]

figure_dir = abs_path("figures_v2/final")
actual_figure_files = sorted(
    name for name in os.listdir(figure_dir)
    if os.path.isfile(os.path.join(figure_dir, name))
)

if actual_figure_files != expected_figure_files:
    raise SystemExit(
        "figures_v2/final must contain exactly: "
        + ", ".join(expected_figure_files)
        + " | found: "
        + ", ".join(actual_figure_files)
    )

fig_rows = []
for figure_id, figure_file in zip(expected_figure_ids, expected_figure_files):
    rel = f"figures_v2/final/{figure_file}"
    path = abs_path(rel)
    status = ensure_nonempty(path)
    if status:
        raise SystemExit(f"Figure is {status}: {rel}")
    fig_rows.append((figure_id, rel, os.path.getsize(path), md5sum(path)))

fig_manifest_rel = "results_v2/fig_manifest.tsv"
fig_manifest_path = abs_path(fig_manifest_rel)
with open(fig_manifest_path, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["figure_id", "path", "bytes", "md5"])
    writer.writerows(fig_rows)

key_outputs = [
    ("metadata_manifest", "results_v2/metadata/paired_manifest.tsv"),
    ("de_results", "results_v2/deseq2/deseq2_paired_v2_results.tsv"),
    ("de_samples_used", "results_v2/deseq2/deseq2_paired_v2_samples_used.tsv"),
    ("de_vst", "results_v2/deseq2/deseq2_paired_v2_vst.tsv"),
    ("de_diagnostics", "results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv"),
    ("de_size_factors", "results_v2/deseq2/deseq2_paired_v2_size_factors.tsv"),
    ("de_session_info", "results_v2/deseq2/sessionInfo_paired_v2.txt"),
    ("enrichment_hallmark", "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"),
    ("enrichment_go_bp", "results_v2/enrichment/go_bp_ora_paired_v2.tsv"),
    ("enrichment_go_bp_representative", "results_v2/enrichment/go_bp_ora_representative_v2.tsv"),
    ("enrichment_go_bp_tumor_higher", "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv"),
    ("enrichment_go_bp_normal_higher", "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv"),
    ("enrichment_go_bp_tumor_higher_representative", "results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv"),
    ("enrichment_go_bp_normal_higher_representative", "results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv"),
    ("enrichment_diagnostics", "results_v2/enrichment/enrichment_diagnostics_v2.tsv"),
    ("enrichment_session_info", "results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt"),
    ("robustness_metrics", "results_v2/robustness/analysis_metrics_v2.tsv"),
    ("robustness_thresholds", "results_v2/robustness/de_threshold_sensitivity_v2.tsv"),
    ("robustness_prefilter", "results_v2/robustness/low_count_prefilter_sensitivity_v2.tsv"),
    ("robustness_lfc_agreement", "results_v2/robustness/lfc_agreement_v2.tsv"),
    ("robustness_pca", "results_v2/robustness/pca_outlier_summary_v2.tsv"),
    ("robustness_top_de", "results_v2/robustness/top_de_genes_v2.tsv"),
    ("robustness_cohort", "results_v2/robustness/cohort_inclusion_summary_v2.tsv"),
    ("pca_sample_diagnostics", "results_v2/qc/pca_sample_diagnostics_v2.tsv"),
    ("robustness_session_info", "results_v2/robustness/sessionInfo_robustness_v2.txt"),
    ("session_info", "results_v2/sessionInfo.txt"),
    ("fig_manifest", "results_v2/fig_manifest.tsv"),
    ("figure_F01", "figures_v2/final/F01.png"),
    ("figure_F02", "figures_v2/final/F02.png"),
    ("figure_F03", "figures_v2/final/F03.png"),
    ("figure_F04", "figures_v2/final/F04.png"),
    ("figure_F05", "figures_v2/final/F05.png"),
    ("figure_F06", "figures_v2/final/F06.png"),
    ("figure_F07", "figures_v2/final/F07.png"),
    ("figure_F01_pdf", "figures_v2/vector/F01_qc_library_size_pairs.pdf"),
    ("figure_F02_pdf", "figures_v2/vector/F02_qc_pca_pairs.pdf"),
    ("figure_F03_pdf", "figures_v2/vector/F03_de_ma_paired_v2.pdf"),
    ("figure_F04_pdf", "figures_v2/vector/F04_de_volcano_paired_v2.pdf"),
    ("figure_F05_pdf", "figures_v2/vector/F05_de_heatmap_top40_paired_v2.pdf"),
    ("figure_F06_pdf", "figures_v2/vector/F06_bio_hallmark_nes_paired_v2.pdf"),
    ("figure_F07_pdf", "figures_v2/vector/F07_bio_go_bp_dotplot_paired_v2.pdf"),
    ("run_log", log_rel),
]

missing_or_empty = []
output_rows = []
for output_id, rel in key_outputs:
    path = abs_path(rel)
    status = ensure_nonempty(path)
    if status:
        missing_or_empty.append(f"{rel} ({status})")
        continue
    output_rows.append((output_id, rel, os.path.getsize(path), md5sum(path)))

if missing_or_empty:
    raise SystemExit(
        "Missing/empty required outputs: " + "; ".join(missing_or_empty)
    )

output_manifest_path = abs_path("results_v2/output_manifest.tsv")
with open(output_manifest_path, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["output_id", "path", "bytes", "md5"])
    writer.writerows(output_rows)

if ensure_nonempty(output_manifest_path):
    raise SystemExit("results_v2/output_manifest.tsv is missing or empty.")
PY

# OncoRNA: Paired Breast Tumor-Normal RNA-seq

OncoRNA is a student-led bioinformatics learning and portfolio project based on public GEO
dataset GSE306117. It asks which gene-expression differences and biological programs are
associated with Tumor tissue after controlling patient identity with matched Normal samples.
The project is reproducible and carefully bounded, but it is not a clinical or biomarker study.

Final frozen version: `oncorna-v1.0-final`.

## Primary Workflow (Maintained)

The maintained and recommended analysis path in this repository is the v2 paired pipeline:

```bash
bash scripts/run_v2.sh
```

This command is the default entry point. It regenerates v2 outputs from scratch with strict
validation, including:

- Paired DE design: `~ patient_id + condition_main`
- Deterministic controls (pipeline seed + fixed methods)
- Final figure contract: `figures_v2/final/F01.png` to `figures_v2/final/F07.png`
- Required output manifests and checksums

## Project Summary

This repository reproduces a paired RNA-seq analysis of breast cancer samples from GEO
accession GSE306117. The main question is differential expression between Tumor and matched
Normal tissue, followed by Hallmark GSEA and GO Biological Process ORA. The v2 workflow
operates on 21 matched Tumor/Normal pairs (42 samples) selected through manifest rules.

## Repository Layout

Key directories:

- `data/`: raw inputs, processed counts, metadata
- `scripts/`: pipeline scripts (`00-metadata`, `01-qc`, `02-de`, `03-pathways`, `04-figures`)
- `dashboard/`: read-only Streamlit presentation of canonical outputs
- `results_v2/`: validated v2 analysis outputs, manifests, logs
- `figures_v2/`: final v2 figures (`final/`) and archived prior runs (`archive/`)

Legacy output directories retained for reference:

- `results/`
- `figures/`

Gitignored data directories:

- `data/raw/`
- `data/processed/`

## Requirements

R environment (renv):

```r
install.packages("renv")
renv::restore()
```

Python environment (`.venv` recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.lock
```

Notes:

- `requirements.lock` is the pinned Python dependency set used in this repo.
- `requirements.txt` is a minimal unpinned list.
- `scripts/run_v2.sh` uses `.venv/bin/python3` if present; otherwise it falls back to
  `python3` on your `PATH`.
- `scripts/run_v2.sh` preflight requires `renv.lock`,
  `data/raw/GSE306117_series_matrix.txt`, and `data/raw/GSE306117_RAW/`.

### Quickstart

```bash
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.lock
Rscript --vanilla -e 'if(!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore(prompt=FALSE)'
bash scripts/run_v2.sh
```

### Smoke Test

```bash
bash tests/smoke_test_v2.sh
```

This smoke test runs the maintained v2 workflow and validates required outputs, manifests,
and the exact final figure contract (`F01.png` to `F07.png`).

### Local Study and Presentation Dashboard

The dashboard reads canonical `results_v2/` tables and `figures_v2/final/F01.png` through
`F07.png`. It does not recompute or modify the analysis.

```bash
./.venv/bin/python -m pip install -r dashboard/requirements.txt
bash run_dashboard.sh
```

Open `http://127.0.0.1:8502`. This is a local demonstration address, not a public link.
Override the port without disturbing another local app with:

```bash
ONCORNA_PORT=8510 bash run_dashboard.sh
```

Verify startup with:

```bash
./.venv/bin/python tests/smoke_test_dashboard.py
```

Professor-facing guides:

- `docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md`
- `docs/ONCORNA_MSC_PORTFOLIO_SUMMARY.md`
- `docs/ONCORNA_PROJECT_UNDERSTANDING.md`
- `docs/ONCORNA_VIVA_SHEET.md`
- `docs/ONCORNA_RESULTS_AND_ROBUSTNESS.md`
- `docs/ONCORNA_CODE_LEARNING_MAP.md`
- `docs/ONCORNA_CLAIM_TRACEABILITY.md`
- `docs/ONCORNA_FUTURE_STUDY_GUIDE.md`

### Tested Software Environment

- macOS arm64
- R 4.5.2
- DESeq2 1.50.2
- Python 3.14.2
- Streamlit 1.59.2
- pandas 3.0.0

Exact R and Python dependencies are recorded in `renv.lock`, `requirements.lock`, and the
session-information files under `results_v2/`.

## What `run_v2.sh` Does (A-E)

- **Step A: metadata + paired manifest**
  - Runs metadata parsing, HTSeq count merging, QC summaries/plots, and paired manifest
    construction.
  - Validates paired manifest structure and required columns, then writes
    `results_v2/metadata/paired_manifest.tsv`.
- **Step B: paired DESeq2 differential expression**
  - Runs `scripts/02-de/01_deseq2_paired_v2.R` with design:
    `design = ~ patient_id + condition_main`.
  - Validates DE columns and that `samples_used` matches the paired manifest cohort.
- **Step C: enrichment**
  - Runs Hallmark GSEA and GO BP ORA via
    `scripts/03-pathways/01_enrichment_paired_v2.R`.
  - Determinism controls include:
    - seed from `PIPELINE_SEED`
    - fixed GSEA method: `fgseaMultilevel` with `BiocParallel::SerialParam()`
    - required DE input columns: `log2FoldChange_shrunk`, `padj`, `pvalue`
- **Step D: figures and naming contract**
  - Runs all v2 figure scripts.
  - Copies/standardizes final set to exactly:
    `figures_v2/final/F01.png` to `figures_v2/final/F07.png`.
  - Archives non-canonical files out of `figures_v2/final/`.
- **Step E: strict final validation + manifests + logs**
  - Enforces required outputs as non-empty.
  - Enforces exact final figure contract (F01-F07 only).
  - Writes checksum manifests and a timestamped run log.

## Outputs (Authoritative)

These are required outputs enforced by `scripts/run_v2.sh`.

### Metadata

- `results_v2/metadata/paired_manifest.tsv`

### Differential Expression

- `results_v2/deseq2/deseq2_paired_v2_results.tsv`
- `results_v2/deseq2/deseq2_paired_v2_samples_used.tsv`
- `results_v2/deseq2/deseq2_paired_v2_vst.tsv` (exploratory figures only)
- `results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv`
- `results_v2/deseq2/sessionInfo_paired_v2.txt`

### Enrichment

- `results_v2/enrichment/hallmark_gsea_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_representative_v2.tsv`
- `results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv`
- `results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv`
- `results_v2/enrichment/enrichment_diagnostics_v2.tsv`
- `results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt`

### Robustness

- `results_v2/robustness/analysis_metrics_v2.tsv`
- `results_v2/robustness/de_threshold_sensitivity_v2.tsv`
- `results_v2/robustness/low_count_prefilter_sensitivity_v2.tsv`
- `results_v2/robustness/lfc_agreement_v2.tsv`
- `results_v2/robustness/pca_outlier_summary_v2.tsv`

### Figures (Final Contract)

- `figures_v2/final/F01.png` (paired cohort library size QC)
- `figures_v2/final/F02.png` (paired cohort PCA QC)
- `figures_v2/final/F03.png` (DE MA plot)
- `figures_v2/final/F04.png` (DE volcano plot)
- `figures_v2/final/F05.png` (top DE genes heatmap)
- `figures_v2/final/F06.png` (Hallmark GSEA NES summary)
- `figures_v2/final/F07.png` (directional GO BP ORA summary)

### Manifests, Session, Logs

- `results_v2/fig_manifest.tsv`
- `results_v2/output_manifest.tsv`
- `results_v2/sessionInfo.txt`
- `results_v2/logs/run_v2_*.log`

## Legacy Workflows (Deprecated)

Older non-v2 workflows are retained for historical comparison and reference only.
They are not the recommended analysis path and are not maintained as the primary pipeline.

- Legacy scripts mainly target `results/` and `figures/` output trees.
- Legacy DE and enrichment scripts in `scripts/02-de/` and `scripts/03-pathways/` are
  preserved but are not used by `bash scripts/run_v2.sh`.
- Use `bash scripts/run_v2.sh` for reproducible, validated paired analysis.

## Citation

Data source: NCBI GEO accession GSE306117.

## Scientific Limitations

The analysis uses one bulk-tissue cohort and does not model cell composition directly. Gene
mapping is incomplete, enrichment terms overlap, and the work has no independent cohort,
protein assay, functional experiment, clinical outcome model, or clinical validation. Results
support cohort-specific biological hypotheses, not causation, diagnosis, prognosis, treatment
recommendations, or validated biomarkers.

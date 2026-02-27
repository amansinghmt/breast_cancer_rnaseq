# Breast Cancer RNA-seq (GSE306117)

## Project Summary
This repository reproduces a paired RNA-seq analysis of breast cancer samples from GEO
accession GSE306117. The primary question is differential expression between Tumor and
matched Normal tissue using a paired design, followed by pathway interpretation with
Hallmark GSEA and GO Biological Process ORA. The stable v2 workflow operates on 21 matched
Tumor/Normal pairs (42 samples) selected through the manifest rules in the pipeline.

## Reproducibility (One Command)
Run the full recommended pipeline from scratch:

```bash
bash scripts/run_v2.sh
```

This command regenerates all v2 outputs and performs strict validation. It exits with a
non-zero status if required outputs are missing/empty or manifest/figure contracts fail.

## Repository Layout
Key directories:

- `data/`: raw inputs, processed counts, metadata
- `scripts/`: pipeline scripts (`00-metadata`, `01-qc`, `02-de`, `03-pathways`, `04-figures`)
- `results_v2/`: validated v2 analysis outputs, manifests, logs
- `figures_v2/`: final figures (`final/`) and archived prior runs (`archive/`)

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
- `requirements.txt` also exists as a minimal unpinned list.
- `scripts/run_v2.sh` uses `.venv/bin/python3` if present; otherwise it falls back to
  `python3` on your `PATH`.
- `scripts/run_v2.sh` preflight also requires:
  - `renv.lock`
  - `data/raw/GSE306117_series_matrix.txt`
  - `data/raw/GSE306117_RAW/`

### Quickstart
```bash
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.lock
Rscript --vanilla -e 'if(!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore(prompt=FALSE)'
bash scripts/run_v2.sh
```

## What The Pipeline Does (A–E)
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
  - Determinism controls in code:
    - seed from `PIPELINE_SEED`
    - fixed GSEA method: `fgseaMultilevel` with `BiocParallel::SerialParam()`
    - required DE input columns for enrichment prep include
      `log2FoldChange_shrunk`, `padj`, `pvalue`.
- **Step D: figures and naming contract**
  - Generates figure sources via all v2 figure scripts.
  - Copies/standardizes final set to exactly:
    `figures_v2/final/F01.png` ... `F07.png`.
  - Archives non-canonical files out of `figures_v2/final/`.
- **Step E: strict final validation + manifests + logs**
  - Enforces required outputs as non-empty.
  - Enforces figure folder contract (exactly F01–F07).
  - Writes checksum manifests and timestamped run log.

## Outputs (Authoritative)
These are the required outputs enforced by `scripts/run_v2.sh`.

### Metadata
- `results_v2/metadata/paired_manifest.tsv`

### Differential Expression
- `results_v2/deseq2/deseq2_paired_v2_results.tsv`
- `results_v2/deseq2/deseq2_paired_v2_samples_used.tsv`
- `results_v2/deseq2/sessionInfo_paired_v2.txt`

### Enrichment
- `results_v2/enrichment/hallmark_gsea_paired_v2.tsv`
- `results_v2/enrichment/go_bp_ora_paired_v2.tsv`
- `results_v2/enrichment/sessionInfo_enrichment_paired_v2.txt`

### Figures (Final Publication Contract)
- `figures_v2/final/F01.png` (paired cohort library size QC)
- `figures_v2/final/F02.png` (paired cohort PCA QC)
- `figures_v2/final/F03.png` (DE MA plot)
- `figures_v2/final/F04.png` (DE volcano plot)
- `figures_v2/final/F05.png` (top DE genes heatmap)
- `figures_v2/final/F06.png` (Hallmark GSEA NES summary)
- `figures_v2/final/F07.png` (GO BP ORA summary)

### Manifests, Session, Logs
- `results_v2/fig_manifest.tsv`
- `results_v2/output_manifest.tsv`
- `results_v2/sessionInfo.txt`
- `results_v2/logs/run_v2_*.log`

## Notes / Caveats
- Counts used for analysis are built by merging GEO supplementary HTSeq count files from
  `data/raw/GSE306117_RAW` via `scripts/01-qc/02_merge_htseq_counts.py`.
- Metadata are parsed from `data/raw/GSE306117_series_matrix.txt`.
- Determinism in v2 includes a fixed pipeline seed (`PIPELINE_SEED`, default `20260227`)
  and deterministic enrichment execution (`fgseaMultilevel` with `SerialParam`).
- A legacy non-v2 workflow still exists under older `results/` and non-v2 scripts, but
  `scripts/run_v2.sh` is the default and recommended path.

## Citation
Data source: NCBI GEO accession GSE306117.

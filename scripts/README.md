# Scripts Overview

## Primary Entry Point

```bash
bash scripts/run_v2.sh
```

This is the default and recommended command for reproducible analysis.

## Directory Roles

### `00-metadata/`

Builds metadata and sample manifests used by the v2 paired pipeline.

### `01-qc/`

Produces counts/QC inputs used by v2.

- `02_merge_htseq_counts.py` and `03_qc_plots.py` are used by `run_v2.sh`.
- `01_extract_counts.py` is retained for reference and is not used by `run_v2.sh`.

### `02-de/`

Contains both active v2 and legacy scripts.

- v2 active script: `01_deseq2_paired_v2.R`
- legacy scripts: `01_deseq2_tumor_vs_normal.R`, `02_volcano_heatmap.R`

### `03-pathways/`

Contains both active v2 and legacy scripts.

- v2 active script: `01_enrichment_paired_v2.R`
- legacy script: `01_enrichment_hallmark_go.R`

### `04-figures/`

Contains v2 publication-style figure scripts used by `run_v2.sh`.

- Canonical final figures: `figures_v2/final/F01.png` to `figures_v2/final/F07.png`

### `05-reporting/`

Builds predefined robustness and presentation summaries. These outputs do not replace the
canonical paired DE result.

## Legacy Scripts

Legacy scripts are retained for historical comparison/reference. They are not deleted, but
they are not the recommended or maintained execution path.

Use `bash scripts/run_v2.sh` for the primary paired, deterministic, validated workflow.

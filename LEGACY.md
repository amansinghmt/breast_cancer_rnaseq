# Legacy vs v2 Workflows

## Use this command for analysis

```bash
bash scripts/run_v2.sh
```

This is the maintained workflow for paired analysis with deterministic controls and strict
final validation.

## v2 (maintained)

- Paired DE design: `~ patient_id + condition_main`
- Deterministic seed/method controls
- Final figures validated as: `figures_v2/final/F01.png` to `figures_v2/final/F07.png`
- Outputs under `results_v2/` and `figures_v2/`

## Legacy (deprecated, retained for reference)

Legacy scripts remain in the repository for historical comparison only. They are not the
recommended path and are not used by `scripts/run_v2.sh`.

Examples:

- `scripts/02-de/01_deseq2_tumor_vs_normal.R`
- `scripts/02-de/02_volcano_heatmap.R`
- `scripts/03-pathways/01_enrichment_hallmark_go.R`
- `scripts/01-qc/01_extract_counts.py`
- `scripts/make_verification_log.sh`

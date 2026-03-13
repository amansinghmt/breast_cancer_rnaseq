# Differential expression (DE)

## Maintained workflow (v2 paired)

Run the repository's primary pipeline entry point:

```bash
bash scripts/run_v2.sh
```

This calls:

- `scripts/02-de/01_deseq2_paired_v2.R`

v2 paired design:

- `~ patient_id + condition_main`

v2 DE outputs are written under `results_v2/`.

## Legacy workflows (deprecated, retained for reference)

The scripts below are preserved for historical comparison and are not used by
`scripts/run_v2.sh`:

- `scripts/02-de/01_deseq2_tumor_vs_normal.R`
  - Legacy design: `~ lfs_status + condition`
  - Legacy outputs under `results/` and `figures/`
- `scripts/02-de/02_volcano_heatmap.R`
  - Legacy visualization path for legacy DE outputs

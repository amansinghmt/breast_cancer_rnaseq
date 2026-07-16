# Pathway enrichment

## Maintained workflow (v2 paired)

Run the repository's primary pipeline entry point:

```bash
bash scripts/run_v2.sh
```

This calls:

- `scripts/03-pathways/01_enrichment_paired_v2.R`

v2 enrichment outputs are written under `results_v2/enrichment/`.

## Legacy workflow (deprecated, retained for reference)

The script below is preserved for historical comparison and is not used by
`scripts/run_v2.sh`:

- `scripts/03-pathways/01_enrichment_hallmark_go.R`
  - Legacy outputs under `results/enrichment/` and `figures/enrichment/`

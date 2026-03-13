# Metadata scripts

## Primary Workflow

```bash
bash scripts/run_v2.sh
```

`run_v2.sh` calls both metadata scripts in this directory:

- `01_make_metadata.py` -> `data/metadata/metadata.tsv`
- `02_make_sample_manifest_v2.py` -> `data/metadata/sample_manifest.tsv`

## Inputs

- `data/raw/GSE306117_series_matrix.txt`
- `results/qc/qc_summary.tsv` (for manifest construction)

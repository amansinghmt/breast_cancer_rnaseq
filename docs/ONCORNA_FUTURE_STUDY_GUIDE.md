# OncoRNA Future Study Guide

## What this project does

OncoRNA studies paired bulk RNA-seq data from GEO GSE306117. It compares breast Tumor and matched Normal tissue from 21 patients using a paired DESeq2 model, then connects gene-level results to Hallmark GSEA and directional GO Biological Process ORA. The project is designed for learning, reproducible rerunning, and careful presentation rather than clinical use.

## Seven pipeline stages

1. **Metadata:** parse GEO sample descriptions and patient/tissue labels.
2. **Count matrix:** extract HTSeq gene counts, remove technical rows, and merge samples.
3. **QC and cohort:** calculate retained count depth and select one valid Tumor-Normal pair per patient.
4. **Paired differential expression:** fit `~ patient_id + condition_main` in DESeq2.
5. **Enrichment:** run Hallmark GSEA plus combined-supplementary and directional GO ORA.
6. **Figures and summaries:** produce F01-F07 and robustness tables.
7. **Validation:** check outputs, hashes, software sessions, dashboard behavior, and reproducibility.

## Most important files

| Role | Files |
|---|---|
| Inputs | `data/raw/`, `data/metadata/`, `data/processed/counts.tsv` |
| Cohort scripts | `scripts/00-metadata/`, `scripts/01-qc/` |
| DESeq2 | `scripts/02-de/01_deseq2_paired_v2.R` |
| Enrichment | `scripts/03-pathways/01_enrichment_paired_v2.R` |
| Figures | `scripts/04-figures/` |
| Reporting | `scripts/05-reporting/01_build_robustness_summaries.R` |
| Pipeline entry point | `scripts/run_v2.sh` |
| Main outputs | `results_v2/` |
| Final figures | `figures_v2/final/F01.png` through `F07.png` |
| Dashboard | `dashboard/app.py`, `dashboard/content.py`, `run_dashboard.sh` |
| Scientific explanation | `docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md` |
| Viva preparation | `docs/ONCORNA_VIVA_SHEET.md` |
| Claim evidence | `docs/ONCORNA_CLAIM_TRACEABILITY.md` |

## Concepts to study first

1. RNA-seq count matrices and sample metadata.
2. Matched samples and within-patient comparisons.
3. Sequencing depth and DESeq2 size-factor normalization.
4. Negative-binomial mean and dispersion.
5. The model formula `~ patient_id + condition_main`.
6. Log2 fold change, standard error and the Wald statistic.
7. P-values, Benjamini-Hochberg adjustment and FDR.
8. VST, PCA and row-wise z-scores.
9. GSEA versus ORA and the importance of the background universe.

## Concepts that can wait

- Alternative shrinkage estimators and Bayesian details.
- Single-cell and spatial RNA-seq.
- Cell-type deconvolution.
- Survival models and clinical prediction.
- Machine learning and classifiers.
- Workflow containers or cluster computing.

These may become useful later, but they are not required to explain the current project accurately.

## Code-reading exercises

1. Trace one included sample from `sample_manifest.tsv` into the DESeq2 sample table.
2. In the DESeq2 script, identify the reference level, design formula, contrast, Wald test and shrinkage call.
3. Recalculate the 506 Tumor-higher and 1,130 Normal-higher strict-gene counts from the DE table without changing it.
4. Explain why VST values are used for PCA but raw counts are used for DESeq2.
5. Follow the Wald statistic from the DE table into the Hallmark ranking.
6. Compare combined GO ORA with the two directional GO tables and explain why the combined analysis has no direction.
7. Read one figure script and write down its source table, selection rule and limitation.
8. Add a temporary local dashboard explanation, confirm no result hash changes, then discard the experiment.

## Running the dashboard

```bash
bash run_dashboard.sh
```

The default local address is `http://127.0.0.1:8502`. To use another port:

```bash
ONCORNA_PORT=8510 bash run_dashboard.sh
```

This is a local demonstration address, not a public website.

## Running validation

```bash
./.venv/bin/python tests/validate_oncorna_outputs.py
./.venv/bin/python tests/smoke_test_dashboard.py
```

The full pipeline smoke test reruns the scientific analysis and should be used only when the analysis environment and raw inputs are available:

```bash
bash tests/smoke_test_v2.sh
```

## What must not be changed accidentally

- The 21-patient, 42-sample paired cohort.
- The count matrix or sample ordering.
- `~ patient_id + condition_main`.
- Normal reference and Tumor-versus-Normal contrast.
- Wald testing, BH adjustment and independent filtering.
- `lfcShrink(type="normal")` and VST `blind=FALSE`.
- The reporting thresholds.
- Hallmark ranking direction and GO universe rules.
- F01-F07 or their source tables without a new scientific review.
- Claims that distinguish association from causation or clinical proof.

Before changing any of these, create a new branch, state the scientific reason, preserve the current outputs, and rerun the full validation suite.

## One-to-two-year revisit checklist

- Read the final report and viva sheet before reading code.
- Confirm the final tag `oncorna-v1.0-final` and repository status.
- Restore the pinned Python and R environments.
- Run output and dashboard validation before experiments.
- Re-explain the model, contrast, FDR, shrinkage, GSEA and ORA in your own words.
- Compare current understanding with the “What I am continuing to learn” list.
- Use a new branch for any experiment and keep the frozen outputs unchanged.
- Check whether a suitable independent matched cohort is available.
- Separate exploratory learning changes from claims intended for presentation.
- Record what changed, why it changed, and which tests passed.

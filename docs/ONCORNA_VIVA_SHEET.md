# OncoRNA Viva Sheet

## 30-second explanation

OncoRNA analyzes public breast RNA-seq data from GEO GSE306117. I selected 21 patients with one quality-controlled Tumor and matched Normal sample each. I used a paired DESeq2 model, `~ patient_id + condition_main`, so each patient acts as their own baseline. The analysis found 6,315 genes with adjusted p-value below 0.05 and 1,636 with both adjusted p-value below 0.05 and at least two-fold shrunken change. Hallmark and GO analyses support cell-cycle and mitotic-associated Tumor-side hypotheses. The results are associations in one bulk-tissue cohort, not clinical or causal proof.

## Two-minute explanation

The project asks which genes and biological programs differ between breast Tumor and matched Normal tissue in GSE306117. I parsed 74 GEO samples and merged gene-level HTSeq counts into a 60,617-gene matrix. The cohort rule requires at least one million retained gene counts and exactly one Tumor and one Normal from the same patient. This leaves 21 matched patients, or 42 samples.

DESeq2 models raw counts with a negative-binomial generalized linear model. The design includes patient identity and condition, which controls stable differences between patients before estimating Tumor versus Normal. Normal is the reference; positive log2 fold change means higher in Tumor. P-values are corrected using Benjamini-Hochberg FDR, and log2 fold changes are shrunk using the preserved `normal` method.

There are 6,315 genes at `padj<0.05`; 1,636 also have `|shrunken log2FC|>=1`, including 506 Tumor-higher and 1,130 Normal-higher genes. Hallmark GSEA tests the full ranked list and identifies 35 significant sets. GO ORA tests the 1,636 selected genes against the 30,244-gene statistical universe and identifies 729 significant terms among 5,102 tested. Cell-cycle and mitotic themes are prominent, but bulk-cell composition and other factors may contribute. Independent cohorts and laboratory work are required before mechanistic or clinical claims.

## Pipeline in ten steps

1. Read GEO GSE306117 metadata.
2. Classify Tumor, Normal, surrounding and other tissues.
3. Read 74 HTSeq count files and remove `__*` technical rows.
4. Merge 60,617 gene rows across samples and validate integer counts.
5. Calculate library sizes and apply the 1e6 count minimum.
6. Require one passing Tumor and Normal per patient; select highest-depth duplicates deterministically.
7. Fit paired DESeq2 model `~ patient_id + condition_main`.
8. Apply Wald tests, BH correction and `normal` LFC shrinkage.
9. Run Hallmark GSEA and GO BP ORA with explicit backgrounds/mapping.
10. Generate F01-F07, robustness tables, logs, manifests and checksums.

## Seven figures

| Figure | Meaning | Main caution |
|---|---|---|
| F01 | Paired gene-level count totals and QC minimum | Depth alone is not complete QC |
| F02 | VST PCA and matched-pair movement | Exploratory, not a test |
| F03 | Shrunken effect versus abundance | Association is not mechanism |
| F04 | Effect versus adjusted evidence | Extreme points are not validated biomarkers |
| F05 | Balanced top DE genes as row z-scores | Clustering does not define subtypes |
| F06 | Positive/negative Hallmark NES | Enrichment is not pathway activation proof |
| F07 | Non-redundant GO presentation | Full overlapping table remains authoritative |

## Ten essential statistical concepts

1. Raw counts are discrete and mean-dependent.
2. Size factors normalize sample depth in the model.
3. Negative-binomial dispersion models extra biological variation.
4. Pairing controls patient-specific baseline differences.
5. A Wald statistic compares effect size with uncertainty.
6. A p-value is a null tail probability, not truth probability.
7. BH FDR addresses tens of thousands of simultaneous tests.
8. Log2 fold change describes direction and magnitude.
9. Shrinkage stabilizes uncertain effects.
10. VST/PCA/z-scores are exploratory transformations, not the DE model.

## Main numerical results

- 74 metadata/count samples.
- 21 matched patients; 42 included samples.
- 60,617 gene rows tested.
- 6,315 genes with `padj<0.05`.
- 3,637 genes with `padj<0.01`.
- 1,636 genes with `padj<0.05` and `|shrunk LFC|>=1`.
- 506 Tumor-higher and 1,130 Normal-higher under the primary rule.
- 35/50 significant Hallmark sets.
- 729/5,102 significant GO BP terms.
- 30 semantically reduced representative GO terms.

## Strongest justified claims

- The maintained workflow reproducibly analyzes a 21-patient paired public cohort.
- Tumor condition is statistically associated with widespread expression differences after controlling patient baseline.
- Cell-cycle/mitotic-associated gene sets are enriched toward the Tumor-side ranking.
- Broad findings are stable under one predefined low-count sensitivity analysis.

## Unsupported claims

- “This proves these genes cause breast cancer.”
- “This can diagnose or predict breast cancer.”
- “These are validated biomarkers.”
- “The pathways are activated or inhibited mechanistically.”
- “The findings apply to every patient or population.”
- “This suggests a treatment recommendation.”

## Twenty likely viva questions

1. **What is the project?** A paired bulk RNA-seq Tumor-Normal analysis of GSE306117.
2. **Why use this dataset?** It provides public HTSeq counts and matched tissues.
3. **Why paired samples?** Each patient supplies their own baseline, reducing between-person variation.
4. **What is RNA-seq?** Sequencing-based measurement of RNA abundance.
5. **What is a count matrix?** Genes by samples with assigned integer counts.
6. **Why DESeq2?** It handles overdispersed counts, normalization and multifactor designs.
7. **Why not a t-test?** A t-test does not model raw count mean-variance behavior or sequencing depth well.
8. **What does the formula mean?** Patient effects are controlled before estimating condition.
9. **What is the reference?** Normal; positive effects mean higher in Tumor.
10. **What is log2 fold change?** +1 is about two-fold Tumor-higher; -1 is two-fold Normal-higher.
11. **What is a p-value?** A null-model probability of an equally or more extreme statistic.
12. **Why adjust p-values?** Thousands of tests otherwise create many false positives.
13. **What is FDR?** Expected false-discovery proportion among discoveries under assumptions.
14. **What is shrinkage?** Stabilization of uncertain fold changes toward zero.
15. **What does PCA show?** Major variance patterns and possible outliers, not significance.
16. **What is GSEA?** Rank-based testing of coordinated gene-set shifts.
17. **How is ORA different?** ORA uses a selected list; GSEA uses the full ranking.
18. **How reliable is it?** Reproducible and broadly robust internally, but not externally validated.
19. **What is the next validation?** Independent matched cohort, then predefined experimental validation.
20. **How was AI used?** AI assisted implementation/review/docs; the owner must verify and understand every scientific decision.

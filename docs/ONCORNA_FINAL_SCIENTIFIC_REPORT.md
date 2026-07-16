# 1. OncoRNA: Paired Breast Tumor-Normal RNA-seq Analysis

## 2. Abstract

OncoRNA is a reproducible bulk RNA-seq analysis of GEO GSE306117. From 74 samples, it selects 21 patients with one Tumor and one matched Normal sample each. A paired DESeq2 model identifies 6,315 genes with BH-adjusted p-value below 0.05 and 1,636 that also have absolute shrunken log2 fold change at least 1. Hallmark GSEA and directional GO ORA support cell-division-associated Tumor-higher hypotheses and distinct Normal-higher tissue-context patterns. These are cohort-specific associations, not causal, clinical, or biomarker conclusions.

## 3. Biological question

Which gene-expression levels and predefined biological programs differ between Tumor and matched Normal breast tissue after controlling each patient's baseline?

## 4. Why matched Tumor-Normal analysis matters

Patients differ genetically and physiologically. Including `patient_id` absorbs stable between-person baselines, so the condition coefficient estimates a within-patient Tumor-Normal effect. Pairing improves alignment with the sampling design but does not remove batch, cell-composition, or other unmeasured confounding.

## 5. Dataset and cohort

The source is public GEO accession GSE306117. The manifest has 74 samples. The maintained cohort has 42 samples from 21 patients, with exactly 21 Tumor and 21 Normal samples.

## 6. Sample-selection process

Samples are classified from GEO metadata, checked for retained gene-level count depth, and required to form a valid pair. Where duplicate eligible tissue samples exist, the highest-depth passing sample is selected deterministically. This controls reproducibility, not biological representativeness.

## 7. RNA-seq count matrix

The matrix contains 60,617 unique gene rows and 42 included sample columns. Counts are non-negative integers. HTSeq technical summary rows are excluded. The model uses raw counts, not VST or z-scores.

## 8. Quality control

F01 shows retained gene-level count depth on a log scale against the predefined one-million-count minimum; all included samples pass. The range is 1,581,646 to 62,668,225 and the median is 4,496,314. This is one QC dimension only. F02 uses saved DESeq2 VST values; PC1 and PC2 explain 32.5% and 15.1% of variance. PCA labels follow an objective robust-distance/paired-shift rule and are exploratory.

## 9. Paired DESeq2 model

The design is `~ patient_id + condition_main`. `condition_main` has levels `Normal`, `Tumor`; the contrast is Tumor versus Normal. Positive effects mean Tumor-higher and negative effects mean Normal-higher.

## 10. Why patient identity is controlled

Without the patient term, naturally high expression in one person could be confused with a tissue-condition effect. The patient coefficients account for those baselines before estimating the common condition difference.

## 11. Normalization, dispersion and negative-binomial modelling

DESeq2 estimates ratio-method size factors (observed range 0.317 to 11.316), gene-wise dispersions, and a parametric dispersion trend. Its negative-binomial model allows count variance to exceed the mean. Normalization does not make every sample identical; it adjusts exposure-like count-depth differences under DESeq2's assumptions.

## 12. Wald testing

For each gene, DESeq2 divides the estimated Tumor-Normal coefficient by its standard error. The resulting Wald statistic measures standardized evidence against a zero condition coefficient under the fitted model.

## 13. P-values and Benjamini-Hochberg FDR

A p-value is a null-model tail probability, not the probability that a gene is truly differential. Because thousands of genes are tested, `results()` applies Benjamini-Hochberg adjustment. `padj<0.05` is the statistical reporting threshold; it does not guarantee biological importance.

## 14. Independent filtering

DESeq2 uses `alpha=0.1` to optimize independent filtering while final reporting applies `padj<0.05`. There are 8,199 all-zero rows, 22,174 positive-p-value rows removed from multiple-testing adjustment by independent filtering, and 30,244 genes with non-missing adjusted p-values. The selected baseMean threshold is 2.2750419.

## 15. Log2 fold-change shrinkage

`lfcShrink(type="normal")` pulls uncertain effects toward zero and stabilizes ranking/visualization. The primary effect rule uses the shrunken coefficient. Shrinkage changes effect estimates; it does not replace the Wald p-values.

## 16. Differential-expression results

Among 60,617 rows, 6,315 genes have `padj<0.05`. Of these, 1,636 also have `|shrunken log2FC|>=1`: 506 Tumor-higher and 1,130 Normal-higher. F03 separates FDR-only genes (4,679) from strict-effect genes. F04 labels up to four genes per direction using lowest adjusted p-value among strict genes with baseMean at least 10; labels are not validated biomarkers.

## 17. PCA and heatmap interpretation

F02 shows global variance structure, matched shifts, and heterogeneity; it is not a group-separation test. F05 uses the VST matrix and row-wise z-scores for 20 objectively selected genes per direction. Samples remain Normal then Tumor within each patient. Genes are clustered only within their known direction block. Relative patterns do not define molecular subtypes.

## 18. Hallmark GSEA

Genes are ranked by the DESeq2 Wald statistic after Ensembl-to-symbol mapping and deterministic duplicate-symbol resolution. `fgseaMultilevel` tests 50 Hallmark sets of mapped size 15–500. Thirty-five have `padj<0.05`: 17 positive NES and 18 negative NES. E2F targets, MYC targets V1, G2M checkpoint and mTORC1 signaling are enriched toward the Tumor-higher side. “Enriched toward” does not mean mechanistically activated.

## 19. Directional GO ORA

The old combined ORA used all 1,636 strict genes and was non-directional. It remains supplementary: 729 of 5,102 terms have BH FDR below 0.05. The corrected main interpretation uses two separate lists against the same 30,244-gene statistical universe and 15,233-gene GO-annotated denominator. Tumor-higher ORA maps 324 of 506 genes and finds 259 significant terms among 3,095 tested. Normal-higher ORA maps 712 of 1,130 genes and finds 694 significant terms among 4,686 tested. `clusterProfiler::simplify` uses Wang similarity at cutoff 0.7 for presentation; F07 then applies a documented maximum of two terms per keyword family. Full tables are unchanged.

### Figure-by-figure interpretation

| Figure | Measured | Observed | May mean | Alternative explanations | Cannot prove |
|---|---|---|---|---|---|
| F01 | Retained count depth | All included samples exceed 1e6 | Adequate depth for retained cohort rule | Assignment efficiency and composition differ | Overall sample quality |
| F02 | VST PCA | Broad condition structure plus heterogeneity | Condition contributes to global variance | Batch/cell mixture/unmeasured covariates | Statistical separation |
| F03 | Mean abundance vs shrunken effect | Both directions; 1,636 strict genes | Widespread association | Composition and model assumptions | Mechanism/clinical value |
| F04 | Shrunken effect vs adjusted evidence | Asymmetric directional counts | Candidate genes for follow-up | Annotation and abundance effects | Biomarker validity |
| F05 | Row-z-scored VST for 40 genes | Recurrent paired patterns with heterogeneity | Selected genes track condition | Patient and cell-mixture variation | Subtypes/classifier |
| F06 | Hallmark NES and FDR | Cell-cycle sets toward Tumor-higher ranking | Coordinated program-level association | Gene-set overlap/composition | Pathway activation |
| F07 | Directional GO ratios/count/FDR | Division themes Tumor-higher; tissue-context themes Normal-higher | Different annotation structures by direction | Mapping loss/redundancy/composition | Mechanism or suppression |

## 20. Robustness checks

A predefined low-count filter retains 34,091 genes. It yields 6,290 significant genes versus 6,310 canonical significant genes among retained rows; 6,219 are shared. Strict-effect counts are 1,636 canonical and 1,709 prefiltered; shrunken LFC Spearman correlation is 0.999734. Unshrunk versus shrunk LFC Spearman correlation is 0.812707. No sample exceeds the exploratory robust PC distance threshold of 3.5.

## 21. Results that are stable

The cohort structure, contrast direction, 6,315 FDR count, 1,636 strict count, broad effect directions, and principal cell-cycle-associated Tumor-side enrichment reproduce under the maintained pipeline. The low-count sensitivity model gives very similar effect ordering.

## 22. Results that remain uncertain

Individual gene priority, exact selected-list counts at alternative thresholds, pathway mechanism, and whether Normal-higher themes reflect tissue loss, cell composition, or tumor-cell regulation remain uncertain.

## 23. Limitations

One public bulk-tissue cohort; no independent validation; incomplete gene mapping; overlapping gene sets; selected-list threshold dependence; no explicit batch or cell-composition model; no protein, spatial, single-cell, functional, survival, treatment-response, or clinical-utility evidence.

## 24. External validation needed

Repeat the same paired model in an independent cohort, verify metadata and technical covariates, assess cell composition, and validate selected genes/programs with orthogonal molecular or functional measurements.

## 25. Reproducibility

Run `bash scripts/run_v2.sh`. The pipeline records session information, timestamped logs, canonical output checksums, figure PDFs/PNGs, complete tables, mapping diagnostics, and validation results. The dashboard launches with `bash run_dashboard.sh` and reads outputs without modifying them.

## 26. AI-use disclosure

AI assistance supported code construction, review, testing, visualization, and educational documentation. The project owner remains responsible for understanding the biological question, cohort, model, thresholds, outputs, limitations, and every public claim.

## 27. Justified conclusion

Within the 21-patient paired GSE306117 cohort, Tumor and matched Normal breast tissues have widespread statistically supported expression differences after controlling patient baseline. Hallmark GSEA and directional GO ORA associate Tumor-higher genes with cell-cycle/division programs and Normal-higher genes with several tissue-context programs. These are reproducible cohort-level hypotheses requiring external and experimental validation.

## 28. Glossary

- **baseMean:** DESeq2 mean normalized count for a gene.
- **BH FDR / padj:** multiple-testing-adjusted p-value controlling expected false-discovery proportion under assumptions.
- **Dispersion:** extra count variability beyond a Poisson mean.
- **GSEA:** rank-based gene-set enrichment analysis.
- **GeneRatio:** GO-annotated selected genes in a term divided by all GO-annotated selected genes.
- **NES:** normalized enrichment score; sign follows ranking direction.
- **ORA:** over-representation analysis of a selected gene list against a background.
- **Shrunken log2FC:** stabilized effect estimate used for direction/effect reporting.
- **VST:** variance-stabilizing transformation used only for exploratory plots.
- **Wald statistic:** coefficient divided by standard error.

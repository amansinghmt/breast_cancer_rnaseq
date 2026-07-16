# OncoRNA Viva Sheet

## 30-second explanation

OncoRNA analyzes public bulk RNA-seq counts from GEO GSE306117. I selected 21 patients with one Tumor and one matched Normal sample and fitted the paired DESeq2 model `~ patient_id + condition_main`. I found 6,315 genes at BH-adjusted p-value below 0.05 and 1,636 that also had absolute shrunken log2 fold change at least 1. Hallmark and directional GO analyses associate Tumor-higher genes with cell-cycle/division programs. These are cohort-specific associations, not causal or clinical proof.

## Two-minute explanation

The original dataset has 74 samples. A deterministic manifest applies tissue labels, a one-million retained-gene-count minimum, and a strict one-Tumor/one-Normal pairing rule. This leaves 42 samples from 21 patients. The count matrix has 60,617 genes.

DESeq2 models raw counts with a negative-binomial generalized linear model. The patient term controls stable between-person differences; Normal is the reference, so a positive condition effect means Tumor-higher. Wald p-values are adjusted by Benjamini-Hochberg, and log2 fold changes are shrunk with the `normal` method. At `padj<0.05`, 6,315 genes are significant; at the stricter effect rule, 506 are Tumor-higher and 1,130 Normal-higher.

Hallmark GSEA uses the full Wald-statistic ranking and finds 35 significant sets. GO ORA uses selected strict genes. The earlier combined GO list mixed both directions and was therefore non-directional. The corrected presentation runs separate Tumor-higher and Normal-higher ORAs with the same background. This supports different directional annotation patterns but not mechanism, diagnosis, prognosis, or biomarkers.

## Five-minute explanation

Start with the biological question and matched design. Explain count-matrix validation, sample pairing, DESeq2 normalization and overdispersion, the patient-adjusted model, Wald tests, BH FDR, independent filtering, and shrinkage. Then use F01/F02 for QC and exploratory structure, F03/F04 for gene-level effect/evidence, F05 for selected relative patterns, F06 for rank-based Hallmark results, and F07 for selected-list directional GO results. End with robustness, mapping loss, bulk-cell composition, one-cohort limitation, and the need for independent and experimental validation.

## Pipeline in ten steps

1. Parse GEO sample metadata.
2. Extract and merge HTSeq gene counts; remove technical rows.
3. Validate integer, non-negative, complete counts and unique IDs.
4. Calculate retained gene-level count depth.
5. Apply deterministic tissue, QC and pairing rules.
6. Require one Tumor and one Normal sample per included patient.
7. Fit `~ patient_id + condition_main` in DESeq2.
8. Apply Wald testing, BH FDR, independent filtering and LFC shrinkage.
9. Run Hallmark GSEA plus combined-supplementary and directional GO ORA.
10. Generate F01-F07, diagnostics, robustness tables, logs, checksums and tests.

## Meaning of the seven figures

| Figure | Exact meaning | Boundary |
|---|---|---|
| F01 | Retained gene-level count depth for each matched sample | One QC dimension only |
| F02 | PC1/PC2 of the 2,000 most variable VST genes, with matched lines | Exploratory, not a significance test |
| F03 | Mean abundance versus shrunken effect; FDR-only and strict genes separated | Association is not mechanism |
| F04 | Shrunken effect versus `-log10(padj)`; labels selected symmetrically | Labels are not biomarkers |
| F05 | Twenty Tumor-higher and twenty Normal-higher strict genes as row z-scores | Relative pattern, not subtype discovery |
| F06 | Balanced significant Hallmark NES; positive toward Tumor-higher ranking | Enrichment is not activation |
| F07 | Separate directional GO ORAs using one background | Over-representation is not mechanism |

## Core methods

### Paired model

`patient_id` absorbs patient-specific baselines. `condition_main` estimates the remaining common Tumor-Normal association. Pairing reduces between-person noise but not all confounding.

### DESeq2

DESeq2 estimates size factors, models overdispersed counts with a negative-binomial GLM, estimates/shrinks dispersion, and tests model coefficients. Raw counts enter DE; VST enters PCA/heatmaps only.

### FDR

BH adjustment addresses testing thousands of genes. `padj<0.05` controls an expected false-discovery proportion under assumptions; it does not mean a 95% probability the gene is true.

### Shrinkage

`lfcShrink(type="normal")` pulls noisy effect estimates toward zero. It stabilizes direction/magnitude displays and strict effect filtering; it does not create statistical significance.

### GSEA versus ORA

GSEA uses the complete ranked list and asks whether a set accumulates toward one end. ORA starts with a selected gene list and asks whether a term appears more often than expected against a background. GSEA direction comes from NES; ORA direction requires separate directional input lists.

### Combined versus directional GO

The combined 1,636-gene ORA supports only “processes over-represented among genes differing between Tumor and Normal.” It cannot say Tumor-specific. Separate 506-gene Tumor-higher and 1,130-gene Normal-higher analyses provide directional annotation summaries.

## Main verified results

- 74 manifest samples; 21 matched patients; 42 included samples.
- 60,617 count-matrix genes.
- 6,315 genes at `padj<0.05`.
- 1,636 strict genes: 506 Tumor-higher, 1,130 Normal-higher.
- PCA variance: PC1 32.5%, PC2 15.1%.
- 35/50 significant Hallmark sets: 17 positive and 18 negative NES.
- Combined supplementary GO: 729/5,102 significant/tested.
- Tumor-higher GO: 324/506 mapped; 259/3,095 significant/tested.
- Normal-higher GO: 712/1,130 mapped; 694/4,686 significant/tested.
- Shared GO-annotated background denominator: 15,233 from a 30,244-gene statistical universe.

## Strongest justified claims

- This cohort shows widespread paired expression differences associated with tissue condition after controlling patient baseline.
- E2F, MYC and G2M Hallmark sets are enriched toward the Tumor-higher ranking.
- Tumor-higher strict genes are over-represented in cell-division/chromosome-related GO terms.
- Normal-higher strict genes include circulation, muscle, extracellular-matrix and other tissue-context GO terms.
- Results are reproducible within this pipeline and broadly stable to the predefined low-count sensitivity check.

## Unsupported claims

- Any pathway is mechanistically activated or inhibited.
- A labelled gene is a validated biomarker, driver, diagnostic, prognostic, or treatment target.
- The project predicts patient outcomes or generalizes to all breast cancers.
- Normal-higher annotations prove suppression inside Tumor cells.
- Statistical significance proves biological importance.

## Likely professor questions

| Question | Concise answer |
|---|---|
| Why paired analysis? | It controls each patient's stable baseline and estimates within-patient condition differences. |
| Why not a t-test? | Counts are discrete, depth-dependent and overdispersed; DESeq2 models these properties and multifactor designs. |
| What is the reference? | Normal; positive log2FC is Tumor-higher. |
| Why `alpha=0.1` but report 0.05? | `alpha=0.1` is preserved for DESeq2 independent-filter optimization; reporting separately uses `padj<0.05`. |
| Why shrink LFC? | To stabilize noisy effect magnitudes; p-values remain from the fitted Wald test. |
| Why did combined GO need correction? | It mixed both LFC signs, so its terms had no Tumor/Normal direction. |
| Why different GO tested-term counts by direction? | Different mapped selected genes annotate to different subsets of GO terms, while the universe stays identical. |
| Does PCA prove separation? | No. It is an exploratory variance projection. |
| Is 6,315 too many? | It reflects power, broad tissue differences and FDR; effect-size filtering narrows interpretation to 1,636. |
| What is the biggest limitation? | One bulk-tissue cohort without cell-composition adjustment or external/experimental validation. |

## AI-use explanation

AI assisted code construction, scientific review, figure refinement, tests and learning material. I am responsible for verifying the data, running the workflow, understanding the model and defending every claim and limitation.

## What the owner must personally understand

The cohort rules; why pairing matters; raw counts versus VST/z-scores; size factors and negative-binomial dispersion; model reference/contrast; Wald p-values; BH FDR; independent filtering; LFC shrinkage; GSEA versus ORA; combined versus directional GO; mapping/background denominators; robustness results; and all claim boundaries.

# OncoRNA: Paired Breast Tumor-Normal RNA-seq Analysis

## Abstract

OncoRNA is a reproducible bulk RNA-seq analysis of GEO dataset GSE306117. It asks which genes and biological programs differ between breast Tumor tissue and matched Normal breast tissue from the same patients. The maintained workflow parses GEO metadata, merges 74 HTSeq count files, applies sample-level quality-control and pairing rules, and selects 21 patients with one Tumor and one Normal sample each. Differential expression is estimated with DESeq2 using the paired model `~ patient_id + condition_main`, which controls patient-specific baseline expression while testing the Tumor-versus-Normal contrast. Of 60,617 gene rows, 6,315 have Benjamini-Hochberg adjusted p-value below 0.05; 1,636 also have absolute shrunken log2 fold change at least 1. Hallmark GSEA identifies 35 significant gene sets, while GO Biological Process over-representation analysis identifies 729 significant terms among 5,102 tested terms. Representative pathway results emphasize cell-cycle, mitotic and chromosome-segregation-associated programs. These findings are statistical associations within one bulk-tissue cohort. They do not demonstrate causation, clinical utility, validated biomarkers or patient-level prediction. The repository includes strict cohort checks, deterministic analysis controls, output manifests, checksums, robustness summaries, seven canonical figures and a read-only presentation dashboard.

## Problem and motivation

Cancer tissue can differ from nearby or control tissue in the abundance of thousands of RNA transcripts. Measuring these differences can reveal candidate genes and coordinated biological programs associated with the Tumor state. A long gene list is difficult to interpret, so this project connects gene-level differential expression to Hallmark and Gene Ontology summaries.

The project was built to demonstrate a complete and reproducible bioinformatics workflow: data organization, quality control, paired statistical modeling, multiple-testing correction, pathway analysis, visualization, validation and scientific communication.

## Biological background

RNA-seq estimates transcript abundance by counting sequenced fragments assigned to genes. In bulk RNA-seq, every tissue sample is a mixture of cell types. A measured difference may therefore reflect changes within cells, changes in cell composition, or both.

Matched Tumor and Normal tissue are useful because human patients differ in genetics, physiology, handling and other characteristics. Comparing each patient with themselves reduces between-person variation. It does not eliminate all confounding, but it aligns the model with the matched study design.

## Research question

Which gene-expression levels and predefined biological gene sets differ between breast Tumor and matched Normal breast tissue in the paired GSE306117 cohort?

## Dataset and cohort

- Source: NCBI GEO accession `GSE306117`.
- Raw analysis inputs: GEO series matrix metadata and 74 compressed HTSeq count files.
- Metadata labels: 40 Tumor, 23 Normal and 11 surrounding/contralateral/other samples not assigned to the main contrast.
- Final cohort: 21 patients, 42 samples, exactly one Tumor and one Normal per patient.
- Biological replication: the 21 patients are independent biological replicate pairs. The two tissues within a patient are matched, not independent replicates.
- QC threshold: at least 1,000,000 retained gene-level counts per sample.
- Exclusions: 32 samples are not included. Reasons include non-main tissue, missing library size, low count depth, no QC-passing matched Normal, or a lower-depth duplicate Tumor sample.

Three patients have more than one Tumor sample in the available main-tissue metadata. For patients 126 and 177, the highest-depth passing Tumor is selected because a passing Normal exists. Patient 50 has no passing Normal, so neither Tumor enters the paired cohort. Highest-depth selection is deterministic and reduces low-depth risk, but it discards within-patient tissue heterogeneity. It is a pragmatic cohort rule, not proof that the selected sample is biologically most representative.

`results_v2/metadata/paired_manifest.tsv` contains all 74 rows, including excluded rows. The filename therefore means “manifest used to define the paired cohort,” not “a file containing only included pairs.” Downstream code filters `include_paired == TRUE` and validates the resulting structure.

## Complete workflow

```text
GEO metadata + 74 HTSeq count files
  -> metadata parsing and tissue classification
  -> merged 60,617 x 74 integer count matrix
  -> count-depth QC and paired-cohort rules
  -> 21 matched patients / 42 samples
  -> DESeq2 paired negative-binomial model
  -> gene-level p-values, BH adjusted p-values and shrunken effects
  -> Hallmark GSEA + GO BP over-representation analysis
  -> robustness summaries, F01-F07, PDFs, manifests and checksums
```

| Stage | Input | Operation | Output | Main risk | Maintained file |
|---|---|---|---|---|---|
| Metadata | GEO series matrix | Parse titles and characteristics; classify tissues | `metadata.tsv` | Wrong labels/patient IDs | `scripts/00-metadata/01_make_metadata.py` |
| Counts | 74 HTSeq files | Remove `__*` rows, outer-join by gene, fill absent genes with zero | `counts.tsv` | Duplicate/non-integer/missing counts | `scripts/01-qc/02_merge_htseq_counts.py` |
| QC | Counts + metadata | Count totals, exploratory log2-CPM PCA | QC summary | Low depth/outliers | `scripts/01-qc/03_qc_plots.py` |
| Cohort | Metadata + QC | Apply threshold, matching and duplicate rules | Paired manifest | Unmatched or arbitrary selection | `scripts/00-metadata/02_make_sample_manifest_v2.py` |
| DE | Raw included counts | Paired negative-binomial GLM and Wald test | DE table + VST + diagnostics | Wrong design/reference/contrast | `scripts/02-de/01_deseq2_paired_v2.R` |
| Enrichment | DE table | Hallmark ranked test; GO ORA against tested background | Full and representative pathway tables | Mapping loss/wrong universe/redundancy | `scripts/03-pathways/01_enrichment_paired_v2.R` |
| Presentation | Canonical tables | F01-F07 and PDF equivalents | Figures | Misleading selection or labels | `scripts/04-figures/` |
| Validation | All outputs | Structural checks and MD5 checksums | Output manifests/log | Silent missing or stale outputs | `scripts/run_v2.sh` |

## Methods and why they were chosen

### Raw counts and normalization

The count matrix contains non-negative integers. The included library-size range is 1,581,646 to 62,668,225 retained gene-level counts, with median 4,496,314. DESeq2 estimates sample size factors using its ratio method. This allows samples with different count depth to be compared without converting the inferential model to CPM.

### Paired negative-binomial model

The maintained design is:

```r
design = ~ patient_id + condition_main
```

The patient coefficients absorb baseline differences between patients. `condition_main` estimates the remaining systematic Tumor-versus-Normal difference. `Normal` is the reference level and the requested contrast is `Tumor` versus `Normal`, so a positive log2 fold change means higher expression in Tumor.

DESeq2 estimates size factors, gene-wise dispersions, a mean-dispersion trend and final shrunken dispersions before fitting Wald tests. The model uses raw counts, not VST values.

### P-values, independent filtering and FDR

DESeq2 is run with its established inferential behavior made explicit: Wald test, parametric dispersion fit, ratio size factors, Cook's-distance filtering enabled, Benjamini-Hochberg adjustment and independent filtering. The `results()` optimization alpha remains 0.1 because that is the original DESeq2 default used by the canonical analysis; biological reporting is separately defined at `padj < 0.05`.

- 8,199 genes have zero counts across all 42 included samples and therefore no nominal p-value.
- 22,174 low-abundance genes have a nominal p-value but `NA` adjusted p-value after independent filtering.
- 30,244 genes enter BH adjustment.
- Zero positive-abundance genes have `NA` nominal p-values, so the current output shows no evidence that Cook's-distance filtering suppressed a positive-abundance gene.

### Effect-size shrinkage

`lfcShrink(type="normal")` stabilizes noisy log2 fold-change estimates. The shrinkage method is preserved because changing it would change reported effect sizes. DESeq2 notes that newer alternatives can have less bias; evaluating a replacement belongs in a separately versioned methodological analysis, not a presentation update.

### VST, PCA and heatmap

The workflow now saves `varianceStabilizingTransformation(dds, blind=FALSE)` values. VST is used only for PCA and heatmap presentation. The DE model still uses raw counts. PCA uses the 2,000 most variable VST genes without per-gene unit-variance scaling. The heatmap uses row-wise z-scores, so its colors represent relative expression within each gene.

### Hallmark GSEA

GSEA ranks genes by the DESeq2 Wald statistic. Ensembl IDs are mapped to symbols; when several Ensembl IDs map to one symbol, the row with largest absolute statistic is kept. Exact ties are broken deterministically. `fgseaMultilevel` uses Hallmark sets between 15 and 500 mapped genes. Positive NES points toward Tumor-higher genes; negative NES points toward Normal-higher genes.

Mapping is incomplete: 52,418 Ensembl rows have nominal p-values, 32,720 Ensembl IDs receive at least one symbol, and 33,306 unique symbols are available before duplicate-symbol resolution because some IDs have one-to-many mappings. This loss limits pathway coverage.

### GO Biological Process ORA

The selected list contains 1,636 genes with `padj < 0.05` and `|shrunken log2FC| >= 1`. The input universe contains 30,244 genes with non-NA adjusted p-value. GO annotation reduces these to 1,036 selected genes and 15,233 universe genes in the actual enrichment ratios.

The current raw table writes all 5,102 tested GO terms; 729 have adjusted p-value below 0.05. The earlier 729-row table contained only cutoff-passing terms, which made “729 of 729 significant” appear more surprising than it was. For presentation only, `clusterProfiler::simplify` uses Wang semantic similarity at cutoff 0.7 and retains the best adjusted-p-value representative from redundant groups. Thirty representative terms are saved; no valid raw term is deleted.

## Verified results

| Result | Verified value |
|---|---:|
| Manifest rows | 74 |
| Included Tumor samples | 21 |
| Included Normal samples | 21 |
| Matched patients | 21 |
| Gene rows tested | 60,617 |
| Genes with `padj < 0.05` | 6,315 |
| Genes with `padj < 0.01` | 3,637 |
| `padj < 0.05` and `|shrunken LFC| >= 1` | 1,636 |
| Tumor-higher under primary rule | 506 |
| Normal-higher under primary rule | 1,130 |
| Hallmark sets tested / significant | 50 / 35 |
| GO BP terms tested / significant | 5,102 / 729 |
| Representative GO table | 30 terms |

The top Hallmark results include Tumor-side E2F targets, MYC targets, G2M checkpoint and mTORC1 signaling. Representative GO terms emphasize mitotic nuclear division and chromosome segregation. These are coherent with a cell-cycle/proliferation-associated Tumor-side pattern, but they do not identify a causal driver.

## Interpretation of the seven figures

### F01: paired count-depth QC

Shows included Tumor and Normal gene-level count totals connected by patient and the predefined QC minimum. All included samples exceed the threshold. It does not prove absence of batch effects or sample degradation.

### F02: paired PCA

Shows the first two principal components of VST expression. PC1 and PC2 explain approximately 32.5% and 15.1% of variance. Pair lines show within-patient movement. PCA is exploratory and not a significance test; no sample crosses the predefined robust-distance flag of 3.5.

### F03: MA plot

Shows shrunken effect size against average normalized abundance. Color requires both FDR and effect-size thresholds. It demonstrates broad statistical differences but not biological mechanism.

### F04: volcano plot

Shows shrunken effect direction/magnitude against adjusted evidence. Four labels are selected from each direction by evidence, where available. Labels are not chosen because a gene was already considered biologically interesting.

### F05: heatmap

Shows up to 20 Tumor-higher and 20 Normal-higher genes passing the primary rule. Values are VST-derived row z-scores and samples are ordered as Normal then Tumor within patient. Clustering is descriptive and does not establish molecular subtypes.

### F06: Hallmark GSEA

Shows balanced significant pathways. Bar direction is NES and point size summarizes adjusted evidence. Positive and negative NES describe positions in the ranked list, not proven pathway activation or inhibition.

### F07: GO BP ORA

Shows 15 representative terms from the semantic-similarity-reduced table. Position is GeneRatio, point size is gene count and color is adjusted evidence. The full tested table remains available and should be used for audit.

## Statistical reliability and robustness

Threshold sensitivity is reported rather than optimized:

| Rule | Tumor higher | Normal higher | Total |
|---|---:|---:|---:|
| `padj<0.05`, no effect cutoff | 3,017 | 3,298 | 6,315 |
| `padj<0.05`, `|shrunk LFC|>=1` | 506 | 1,130 | 1,636 |
| `padj<0.05`, `|shrunk LFC|>=1.5` | 35 | 306 | 341 |
| `padj<0.01`, `|shrunk LFC|>=1` | 419 | 1,092 | 1,511 |

A predefined sensitivity model retains genes with count at least 10 in at least two samples. It keeps 34,091 genes. Within retained genes, 6,310 are significant in the canonical result and 6,290 in the prefiltered rerun; 6,219 are significant in both, 91 only in the canonical result and 71 only after prefiltering. Shrunken effects have Spearman correlation 0.9997. This supports broad robustness while showing that borderline calls depend on filtering and multiplicity.

Unshrunk versus shrunken effects have Spearman correlation 0.8127 across 52,418 estimable genes. Only 43 of the top 100 genes by absolute effect overlap, demonstrating why raw extreme effects should not be ranked without shrinkage.

## Limitations

- One public cohort; no external replication.
- Bulk tissue mixes cancer, stromal, immune and normal cell populations.
- Pairing controls patient baseline but not every batch or clinical covariate.
- Highest-depth duplicate selection may discard within-patient heterogeneity.
- Gene and pathway annotation mapping is incomplete and one-to-many.
- GO terms are overlapping and highly redundant.
- Statistical thresholds are conventions, not biological boundaries.
- RNA differences do not establish protein changes, mechanism or clinical usefulness.

## Reproducibility

- R packages are recorded in `renv.lock`; Python dependencies are pinned for the analysis.
- The pipeline uses one seed and serial GSEA execution.
- Sample order, contrast direction, shrinkage method and figure contract are explicit.
- `results_v2/output_manifest.tsv` records output sizes and MD5 checksums.
- `tests/smoke_test_v2.sh` runs the maintained pipeline and validates outputs.
- Legacy `results/`, `figures/` and legacy scripts are not authoritative.

## Justified conclusions

Within the 21-patient paired GSE306117 cohort, Tumor and Normal breast tissues show widespread statistically supported gene-expression differences. Rank-based Hallmark and selected-gene GO analyses support cell-cycle/mitotic-associated Tumor-side hypotheses. These conclusions are cohort-specific associations and require independent and experimental validation.

## Future validation

1. Reproduce the paired model in an independent matched breast cohort.
2. Compare effect directions and pathway NES across cohorts.
3. Check whether cell-composition estimates explain part of the observed pattern.
4. Validate a predefined, limited gene set at RNA and protein level.
5. Use functional experiments before proposing mechanism.

## AI-use disclosure

AI assistance contributed to code construction, scientific-software review, figure presentation and documentation. AI did not generate the source data. The project owner is responsible for verifying the metadata, cohort, model, contrast, thresholds, outputs and claims, and should not present any method they cannot explain.

## Glossary

- **baseMean:** DESeq2 mean normalized count across included samples.
- **Count matrix:** genes by samples table of non-negative integer counts.
- **Dispersion:** extra biological variability beyond Poisson counting variance.
- **FDR / padj:** multiple-testing-adjusted statistical evidence.
- **GSEA:** rank-based gene-set enrichment analysis.
- **GeneRatio:** selected genes in a GO term divided by annotated selected genes.
- **Log2 fold change:** effect direction and magnitude on a base-two scale.
- **Negative binomial:** count distribution allowing variance greater than the mean.
- **NES:** normalized GSEA enrichment score.
- **ORA:** over-representation analysis of a selected gene list.
- **PCA:** dimension reduction summarizing major variance directions.
- **Shrinkage:** regularization of noisy effect estimates toward zero.
- **VST:** variance-stabilizing transformation for exploratory visualization.
- **Wald test:** coefficient estimate divided by its estimated standard error.
- **Z-score:** value centered and scaled relative to one gene's sample distribution.

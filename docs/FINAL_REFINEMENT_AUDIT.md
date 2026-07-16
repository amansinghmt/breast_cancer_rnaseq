# OncoRNA Final Refinement Audit

## Verified current state

- Baseline reviewed: branch `oncorna-scientific-presentation-v1`, commit `048d8ab`.
- Protected pre-existing work: `.gitignore`, `.Rhistory`, `breast_cancer_rnaseq.Rproj`, and `docs/project_explained.md` were not modified or staged by this refinement.
- Cohort: 74 manifest rows; 42 included samples; 21 patients; exactly one Tumor and one Normal sample per included patient.
- Count matrix: 60,617 unique gene rows; raw non-negative integer counts; HTSeq technical rows excluded upstream.
- Primary model: `~ patient_id + condition_main`; Normal reference; Tumor-versus-Normal contrast; Wald test; BH adjustment; independent filtering with `alpha=0.1`; Cook's cutoff enabled; `lfcShrink(type="normal")`.
- Primary results: 6,315 genes at `padj<0.05`; 1,636 also at `|shrunken log2FC|>=1`; 506 Tumor-higher and 1,130 Normal-higher.
- Exploratory transform: DESeq2 VST with `blind=FALSE`; PCA uses the 2,000 most variable VST rows.
- Hallmark: DESeq2 Wald-statistic ranking; positive NES follows the Tumor-higher side; 35/50 sets have `padj<0.05`.

## Discrepancy found

The former GO ORA input was a combined list of 1,636 genes selected by `padj<0.05` and `|shrunken log2FC|>=1`. It therefore mixed 506 Tumor-higher and 1,130 Normal-higher genes. The combined result is non-directional and cannot support Tumor-specific GO claims.

## Mandatory scientific corrections

1. Preserve the 5,102-row combined GO table as supplementary, non-directional evidence.
2. Add separate Tumor-higher and Normal-higher ORAs with the same 30,244-gene tested universe, Ensembl mapping route, GO BP ontology, BH correction, and unrestricted output cutoffs.
3. Present F07 as two directional panels.
4. State mapping loss and background denominators explicitly.
5. Move teaching instructions out of clean figures and into the report/dashboard.
6. Trace every headline result to its source table and script.

## Optional visual improvements accepted

- Objective PCA labels based on robust PCA distance or the two longest matched-pair shifts.
- Explicit FDR-only category in F03/F04.
- Balanced 20-per-direction heatmap with clustering restricted within directional blocks.
- NES bars with compact FDR labels rather than multiple simultaneous evidence encodings.
- Presentation-only GO family cap after semantic reduction to limit repeated near-synonymous themes.

## Changes explicitly rejected

- Replacing the paired DESeq2 model: no implementation error was found.
- Changing DE thresholds to improve appearance.
- Confidence ellipses on PCA: no inferential model justifies them here.
- Unpaired testing, machine learning, survival analysis, biomarker classification, or clinical prediction.
- Calling enrichment “activation,” calling labelled genes validated biomarkers, or treating bulk-tissue differences as causal regulation.
- Reducing the 40-gene heatmap solely for visual cleanliness.

## Expected maintained changes

- `scripts/03-pathways/01_enrichment_paired_v2.R`
- `scripts/04-figures/01_publication_figures_paired_v2.R` through `07_go_bp_dotplot_paired_v2.R`
- `scripts/05-reporting/01_build_robustness_summaries.R`
- `scripts/run_v2.sh`
- `dashboard/app.py`, `dashboard/content.py`
- scientific documentation and targeted validation tests

## Figure source and selection ledger

| Figure | Source | Exact filter / selection | Generator |
|---|---|---|---|
| F01 | paired manifest plus `results/qc/qc_summary.tsv` | included pairs only; predefined 1e6 retained-count line; patient ID order | `01_publication_figures_paired_v2.R` |
| F02 | saved DESeq2 VST | top 2,000 non-zero-variance genes; labels if robust PC distance >3.5 or in two longest pair shifts; maximum six labels | `02_pca_paired_v2.R` |
| F03 | canonical DE table | shrunken LFC; FDR-only separated; up to three labels per strict direction by padj, then absolute LFC, with baseMean >=10 | `03_ma_paired_v2.R` |
| F04 | canonical DE table | shrunken LFC and `-log10(padj)`; up to four labels per strict direction by padj then absolute LFC, preferring baseMean >=10 | `04_volcano_paired_v2.R` |
| F05 | canonical DE, VST, sample manifest | first 20 strict genes per direction by padj then absolute LFC; cluster within, never across, direction blocks | `05_heatmap_paired_v2.R` |
| F06 | complete Hallmark table | eight lowest-padj significant sets per NES direction; fallback only if fewer than 16 significant sets | `06_hallmark_bar_paired_v2.R` |
| F07 | directional representative GO tables | significant terms; Wang cutoff 0.7 upstream; max two terms per documented keyword family; first ten by padj/count/name | `07_go_bp_dotplot_paired_v2.R` |

## Archive

The pre-refinement canonical state was copied to:

- `results_v2/archive/pre_final_refinement_20260717_010933/`
- `figures_v2/archive/pre_final_refinement_20260717_010933/`

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import subprocess

import pandas as pd
import streamlit as st

from content import (
    CANNOT_SUPPORT,
    CAN_SUPPORT,
    FIGURE_CONTENT,
    METHODS,
    PIPELINE_STAGES,
    VIVA_QA,
)

ROOT = Path(__file__).resolve().parents[1]

PATHS = {
    "manifest": ROOT / "results_v2/metadata/paired_manifest.tsv",
    "de": ROOT / "results_v2/deseq2/deseq2_paired_v2_results.tsv",
    "samples": ROOT / "results_v2/deseq2/deseq2_paired_v2_samples_used.tsv",
    "de_diagnostics": ROOT / "results_v2/deseq2/deseq2_paired_v2_diagnostics.tsv",
    "hallmark": ROOT / "results_v2/enrichment/hallmark_gsea_paired_v2.tsv",
    "go": ROOT / "results_v2/enrichment/go_bp_ora_paired_v2.tsv",
    "go_representative": ROOT / "results_v2/enrichment/go_bp_ora_representative_v2.tsv",
    "go_tumor": ROOT / "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv",
    "go_normal": ROOT / "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv",
    "go_tumor_representative": ROOT / "results_v2/enrichment/go_bp_ora_tumor_higher_representative_v2.tsv",
    "go_normal_representative": ROOT / "results_v2/enrichment/go_bp_ora_normal_higher_representative_v2.tsv",
    "enrichment_diagnostics": ROOT / "results_v2/enrichment/enrichment_diagnostics_v2.tsv",
    "metrics": ROOT / "results_v2/robustness/analysis_metrics_v2.tsv",
    "thresholds": ROOT / "results_v2/robustness/de_threshold_sensitivity_v2.tsv",
    "prefilter": ROOT / "results_v2/robustness/low_count_prefilter_sensitivity_v2.tsv",
    "lfc_agreement": ROOT / "results_v2/robustness/lfc_agreement_v2.tsv",
    "pca_outliers": ROOT / "results_v2/robustness/pca_outlier_summary_v2.tsv",
    "cohort_summary": ROOT / "results_v2/robustness/cohort_inclusion_summary_v2.tsv",
    "top_de": ROOT / "results_v2/robustness/top_de_genes_v2.tsv",
    "output_manifest": ROOT / "results_v2/output_manifest.tsv",
}

DOCUMENT_PATHS = {
    "Final scientific report": ROOT / "docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md",
    "MSc portfolio summary": ROOT / "docs/ONCORNA_MSC_PORTFOLIO_SUMMARY.md",
    "Viva preparation sheet": ROOT / "docs/ONCORNA_VIVA_SHEET.md",
}

FIGURE_PATHS = {
    figure_id: ROOT / f"figures_v2/final/{figure_id}.png"
    for figure_id in FIGURE_CONTENT
}

REQUIRED_PATHS = [
    PATHS["manifest"],
    PATHS["de"],
    PATHS["hallmark"],
    PATHS["go"],
    PATHS["go_tumor"],
    PATHS["go_normal"],
    PATHS["metrics"],
    *FIGURE_PATHS.values(),
]


def validate_required_files() -> list[Path]:
    return [path for path in REQUIRED_PATHS if not path.is_file() or path.stat().st_size == 0]


@st.cache_data(show_spinner=False)
def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def metric_map() -> dict[str, float]:
    table = read_tsv(str(PATHS["metrics"]))
    return dict(zip(table["metric"], table["value"], strict=True))


def latest_analysis_date() -> str:
    log_path: Path | None = None
    if PATHS["output_manifest"].is_file():
        manifest = read_tsv(str(PATHS["output_manifest"]))
        run_rows = manifest[manifest["output_id"] == "run_log"]
        if not run_rows.empty:
            log_path = ROOT / str(run_rows.iloc[0]["path"])
    if log_path is None:
        logs = sorted((ROOT / "results_v2/logs").glob("run_v2_*.log"))
        log_path = logs[-1] if logs else None
    if log_path is None:
        return "Unknown"
    timestamp = log_path.stem.removeprefix("run_v2_")
    try:
        return datetime.strptime(timestamp, "%Y%m%d_%H%M%S").strftime("%d %B %Y, %H:%M")
    except ValueError:
        return log_path.name


def current_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return "Unknown"


def render_figure(figure_id: str) -> None:
    content = FIGURE_CONTENT[figure_id]
    st.subheader(f"{figure_id}. {content['title']}")
    st.image(str(FIGURE_PATHS[figure_id]), width="stretch")
    st.markdown(f"**Caption.** {content['caption']}")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**What it shows.** {content['shows']}")
    with col2:
        st.markdown(f"**What it does not prove.** {content['does_not_prove']}")
    st.caption(f"Source: {content['source']} | Script: {content['script']}")


def download_table(label: str, path: Path) -> None:
    st.download_button(
        label,
        data=path.read_bytes(),
        file_name=path.name,
        mime="text/tab-separated-values",
        width="stretch",
    )


def download_document(label: str, path: Path) -> None:
    if not path.is_file():
        st.warning(f"{label} is not available in this checkout.")
        return
    st.download_button(
        label,
        data=path.read_bytes(),
        file_name=path.name,
        mime="text/markdown",
        width="stretch",
    )


st.set_page_config(
    page_title="OncoRNA | Paired breast-cancer RNA-seq",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
      :root { --ink:#16302b; --teal:#187f69; --orange:#d65a1f; --paper:#f7f3ea; --muted:#4a615a; }
      .stApp { background: linear-gradient(145deg, #fbfaf6 0%, #f2f7f3 58%, #fff8ee 100%); color:var(--ink); }
      h1, h2, h3 { font-family: Georgia, 'Times New Roman', serif; color:var(--ink); letter-spacing:-0.02em; }
      [data-testid="stSidebar"] { background:#112e28; }
      [data-testid="stSidebar"] * { color:#f7f3ea; }
      [data-testid="stSidebar"] [role="radiogroup"] label:hover { background:rgba(255,255,255,.10); border-radius:8px; }
      [data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) { background:#f7f3ea; border-radius:8px; }
      [data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) * { color:#112e28 !important; font-weight:700; }
      div[data-testid="stMetric"] { background:#ffffff; border:1px solid #cbdad3; border-radius:14px; padding:14px; box-shadow:0 8px 28px rgba(22,48,43,.06); }
      [data-testid="stMetricLabel"], [data-testid="stMetricLabel"] * { color:#38544c !important; font-weight:650; }
      div[data-testid="stMetricValue"] { color:#102f27 !important; }
      div[data-testid="stExpander"] { background:#ffffff; border-radius:12px; border-color:#bfd2ca; overflow:hidden; }
      div[data-testid="stExpander"] details summary { color:#173c32 !important; background:#eef6f2; font-weight:650; }
      div[data-testid="stExpander"] details[open] summary { color:#ffffff !important; background:#176b59; }
      div[data-testid="stExpander"] details[open] summary svg { fill:#ffffff !important; color:#ffffff !important; }
      div[data-testid="stDataFrame"] { background:#ffffff; border:1px solid #cbdad3; border-radius:10px; }
      [data-testid="stCaptionContainer"] p { color:#4a615a !important; }
      [data-testid="stSidebar"] [data-testid="stCaptionContainer"] p { color:#d9e8e2 !important; }
      a { color:#0f6f5a !important; text-decoration:underline; text-underline-offset:2px; }
      [data-testid="stDownloadButton"] button { background:#176b59; color:#ffffff; border:1px solid #0f5547; font-weight:650; }
      [data-testid="stDownloadButton"] button:hover { background:#0f5547; color:#ffffff; border-color:#0b463b; }
      [data-testid="stDownloadButton"] button:focus { box-shadow:0 0 0 3px rgba(214,90,31,.25); }
      .hero { padding:28px 34px; border-radius:22px; background:linear-gradient(105deg,#133f35,#187f69 68%,#d65a1f); color:white; margin-bottom:22px; }
      .hero h1 { color:white; margin:0; font-size:3.2rem; }
      .hero p { max-width:850px; font-size:1.08rem; margin:.7rem 0 0; color:#eef8f4; }
      .boundary { border-left:5px solid #d65a1f; background:#fff7ee; padding:15px 18px; border-radius:8px; }
      .route { border:1px solid #cbdad3; background:#ffffff; padding:16px 19px; border-radius:12px; }
      .small-source { color:#5b6d67; font-size:.85rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

missing = validate_required_files()
if missing:
    st.error(
        "The OncoRNA dashboard cannot start because required analysis files are missing:\n\n"
        + "\n".join(f"- {path.relative_to(ROOT)}" for path in missing)
        + "\n\nRun `bash scripts/run_v2.sh` from the repository root."
    )
    st.stop()

metrics = metric_map()
manifest = read_tsv(str(PATHS["manifest"]))

sections = [
    "Overview",
    "Dataset and cohort",
    "Pipeline",
    "Quality control",
    "Differential expression",
    "Pathway analysis",
    "Results and interpretation",
    "Methods explained",
    "Limitations and next steps",
    "Viva mode",
]

st.sidebar.title("OncoRNA")
st.sidebar.caption("Paired RNA-seq student project")
section = st.sidebar.radio("Navigate", sections)
st.sidebar.divider()
st.sidebar.caption(f"Latest pipeline log: {latest_analysis_date()}")
st.sidebar.caption(f"Repository commit: {current_commit()}")
st.sidebar.caption("Analysis version: Paired workflow v2")

if section == "Overview":
    st.markdown(
        """
        <div class="hero">
          <h1>OncoRNA</h1>
          <p>A student-led exploration of paired bulk RNA-seq data, asking which genes and biological programs differ between breast Tumor and matched Normal tissue in GEO dataset GSE306117.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    cols = st.columns(5)
    overview_metrics = [
        ("Matched patients", "paired_patients"),
        ("Included samples", "included_samples"),
        ("Genes tested", "genes_tested"),
        ("padj < 0.05", "genes_padj_lt_0.05"),
        ("FDR + effect rule", "genes_padj_lt_0.05_abs_shrunk_lfc_ge_1"),
    ]
    for col, (label, key) in zip(cols, overview_metrics, strict=True):
        col.metric(label, f"{int(metrics[key]):,}")
    st.subheader("Biological question")
    st.write(
        "This analysis asks which gene-expression differences are associated with Tumor tissue when each "
        "sample is compared with matched Normal tissue from the same patient. I used this question to learn "
        "how count data, paired statistical models and pathway analysis fit together in an RNA-seq workflow."
    )
    st.info(
        "Primary model: `~ patient_id + condition_main`. The patient term controls baseline "
        "between-person differences; the condition term estimates the paired Tumor-versus-Normal contrast."
    )
    st.subheader("Verified current result counts")
    a, b, c = st.columns(3)
    a.metric("Hallmark sets with padj < 0.05", f"{int(metrics['hallmark_sets_padj_lt_0.05']):,}")
    b.metric("GO terms with padj < 0.05", f"{int(metrics['go_terms_padj_lt_0.05']):,}")
    c.metric("Combined GO representatives (supplement)", f"{int(metrics['go_representative_terms']):,}")
    st.caption(f"Current presentation generated from the latest successful run: {latest_analysis_date()}.")

    st.subheader("Suggested presentation route")
    st.markdown(
        """
        <div class="route">
        <strong>For a short professor review:</strong> Overview → Dataset and cohort → Pipeline → Quality
        control → Differential expression → Pathway analysis → Results and interpretation → Limitations and
        next steps.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.subheader("How I approached this project")
    st.markdown(
        """
        1. Identified a paired public breast RNA-seq dataset.
        2. Organised the metadata and gene-count files.
        3. Applied cohort and retained-count-depth checks.
        4. Modelled Tumor-Normal differences while controlling patient identity.
        5. Examined gene-level effects.
        6. Connected ranked genes to pathways and biological processes.
        7. Checked robustness, limitations and reproducibility.
        8. Built this dashboard to study and present the work.
        """
    )

    learn_col, continue_col = st.columns(2)
    with learn_col:
        st.subheader("What I learned through this project")
        st.markdown(
            """
            - How RNA-seq count data are organised.
            - Why matched samples can reduce between-patient variation.
            - Why sequencing depth needs to be accounted for.
            - How paired DESeq2 compares Tumor and Normal tissue.
            - Why thousands of p-values need multiple-testing correction.
            - How gene-level differences connect to pathway hypotheses.
            - Why enrichment does not prove a mechanism.
            - Why reproducibility and limitations matter.
            """
        )
    with continue_col:
        st.subheader("What I am continuing to learn")
        st.markdown(
            """
            - Negative-binomial models and dispersion estimation.
            - Paired generalized linear models and independent filtering.
            - Fold-change shrinkage and enrichment statistics.
            - GSEA and GO ORA assumptions.
            - R and Python code structure.
            - Independent validation and stronger biological interpretation.
            """
        )

    st.subheader("Project documents")
    doc_cols = st.columns(3)
    for col, (label, path) in zip(doc_cols, DOCUMENT_PATHS.items(), strict=True):
        with col:
            download_document(label, path)

elif section == "Pipeline":
    st.title("Pipeline")
    st.graphviz_chart(
        """
        digraph oncorna {
          rankdir=LR; graph [bgcolor="transparent", pad=.2];
          node [shape=box, style="rounded,filled", fillcolor="#eef6f2", color="#187f69", fontname="Helvetica"];
          edge [color="#607a72"];
          A [label="GEO metadata +\nHTSeq counts"]; B [label="Metadata +\ncount matrix"];
          C [label="QC + matched\npairs"]; D [label="Paired DESeq2"];
          E [label="Hallmark GSEA +\nGO ORA"]; F [label="F01-F07 +\nvalidation"];
          A -> B -> C -> D -> E -> F;
        }
        """,
        width="stretch",
    )
    for name, operation, output, failure in PIPELINE_STAGES:
        with st.expander(name):
            st.markdown(f"**Operation:** {operation}")
            st.markdown(f"**Output:** `{output}`")
            st.markdown(f"**Point to check:** {failure}")

elif section == "Dataset and cohort":
    st.title("Dataset and cohort")
    cols = st.columns(4)
    cols[0].metric("GEO accession", "GSE306117")
    cols[1].metric("Manifest rows", f"{int(metrics['manifest_rows'])}")
    cols[2].metric("Included samples", f"{int(metrics['included_samples'])}")
    cols[3].metric("Matched patients", f"{int(metrics['paired_patients'])}")
    st.write(
        "The raw metadata contains Tumor, control breast tissue, surrounding tissue and "
        "contralateral/other labels. Only one QC-passing Tumor and one QC-passing Normal sample "
        "per patient enter the paired model."
    )
    cohort_summary = read_tsv(str(PATHS["cohort_summary"]))
    st.subheader("Inclusion and exclusion summary")
    st.dataframe(cohort_summary, width="stretch", hide_index=True)
    included = manifest[manifest["include_paired"].astype(str).str.upper() == "TRUE"]
    st.subheader("Included matched cohort")
    st.dataframe(
        included[["sample_id", "patient_id", "condition_main", "library_size"]]
        .sort_values(["patient_id", "condition_main"]),
        width="stretch",
        hide_index=True,
    )
    st.markdown(
        '<div class="boundary"><strong>Dataset limitation.</strong> This is one bulk-tissue cohort. '
        "Expression includes mixed cell populations and cannot establish patient-level prediction or population-wide generalization.</div>",
        unsafe_allow_html=True,
    )

elif section == "Quality control":
    st.title("Quality control")
    tabs = st.tabs(["F01 Library size", "F02 PCA", "PCA outlier summary"])
    with tabs[0]:
        render_figure("F01")
    with tabs[1]:
        render_figure("F02")
    with tabs[2]:
        pca_outliers = read_tsv(str(PATHS["pca_outliers"]))
        st.write(
            "Samples are ranked by robust distance in PC1-PC2 space. The flag threshold was "
            "defined before inspection at 3.5; it is exploratory and not an exclusion rule."
        )
        st.dataframe(pca_outliers.head(15), width="stretch", hide_index=True)
        st.metric("Flagged samples", f"{int(metrics['pca_exploratory_outliers'])}")

elif section == "Differential expression":
    st.title("Differential expression")
    strict = read_tsv(str(PATHS["thresholds"]))
    primary = strict[(strict["padj_cutoff"] == 0.05) & (strict["abs_shrunken_log2fc_cutoff"] == 1)].iloc[0]
    cols = st.columns(4)
    cols[0].metric("Genes tested", f"{int(metrics['genes_tested']):,}")
    cols[1].metric("padj < 0.05", f"{int(metrics['genes_padj_lt_0.05']):,}")
    cols[2].metric("Tumor-higher", f"{int(primary['tumor_higher']):,}")
    cols[3].metric("Normal-higher", f"{int(primary['normal_higher']):,}")
    with st.expander("Threshold and robustness summaries"):
        st.dataframe(strict, width="stretch", hide_index=True)
        st.markdown("**Predefined low-count sensitivity**")
        st.dataframe(read_tsv(str(PATHS["prefilter"])), width="stretch", hide_index=True)
        st.markdown("**Unshrunk versus shrunken effect agreement**")
        st.dataframe(read_tsv(str(PATHS["lfc_agreement"])), width="stretch", hide_index=True)
    tabs = st.tabs(["F03 MA", "F04 Volcano", "F05 Heatmap"])
    for tab, figure_id in zip(tabs, ["F03", "F04", "F05"], strict=True):
        with tab:
            render_figure(figure_id)
    st.subheader("Searchable top-gene table")
    top_de = read_tsv(str(PATHS["top_de"]))
    query = st.text_input("Filter by Ensembl ID or gene symbol", placeholder="e.g. EPN3 or ENSG...")
    if query:
        mask = top_de["gene_id"].astype(str).str.contains(query, case=False, na=False) | top_de[
            "SYMBOL"
        ].astype(str).str.contains(query, case=False, na=False)
        top_de = top_de[mask]
    st.dataframe(top_de, width="stretch", hide_index=True, height=430)
    download_table("Download DE results table", PATHS["de"])

elif section == "Pathway analysis":
    st.title("Pathway analysis")
    cols = st.columns(3)
    cols[0].metric("Hallmark sets tested", f"{int(metrics['hallmark_sets_tested'])}")
    cols[1].metric("Significant Hallmark sets", f"{int(metrics['hallmark_sets_padj_lt_0.05'])}")
    cols[2].metric("Combined GO supplement terms", f"{int(metrics['go_terms_padj_lt_0.05'])}")
    tabs = st.tabs(["F06 Hallmark", "F07 Directional GO BP", "Combined GO supplement"])
    with tabs[0]:
        render_figure("F06")
        hallmark = read_tsv(str(PATHS["hallmark"]))
        st.dataframe(hallmark, width="stretch", hide_index=True, height=350)
        download_table("Download Hallmark table", PATHS["hallmark"])
    with tabs[1]:
        render_figure("F07")
        diagnostics = dict(
            zip(
                read_tsv(str(PATHS["enrichment_diagnostics"]))["metric"],
                read_tsv(str(PATHS["enrichment_diagnostics"]))["value"],
                strict=True,
            )
        )
        direction_cols = st.columns(2)
        direction_cols[0].metric(
            "Tumor-higher GO terms (BH FDR < 0.05)",
            f"{int(diagnostics['ora_tumor_higher_terms_padj_lt_0.05']):,}",
            help="From 506 strict Tumor-higher genes; 324 map to GO annotations.",
        )
        direction_cols[1].metric(
            "Normal-higher GO terms (BH FDR < 0.05)",
            f"{int(diagnostics['ora_normal_higher_terms_padj_lt_0.05']):,}",
            help="From 1,130 strict Normal-higher genes; 712 map to GO annotations.",
        )
        representative = pd.concat(
            [
                read_tsv(str(PATHS["go_tumor_representative"])),
                read_tsv(str(PATHS["go_normal_representative"])),
            ],
            ignore_index=True,
        )
        st.dataframe(representative, width="stretch", hide_index=True, height=350)
        col1, col2 = st.columns(2)
        with col1:
            download_table("Download full Tumor-higher GO table", PATHS["go_tumor"])
        with col2:
            download_table("Download full Normal-higher GO table", PATHS["go_normal"])
        st.info(
            "Both directional analyses use the same tested-gene background and BH correction. "
            "Semantic and keyword-family reductions affect presentation only, not the full tables."
        )
    with tabs[2]:
        st.warning(
            "This table combines Tumor-higher and Normal-higher genes. It supports only processes "
            "over-represented among genes differing between conditions; it cannot establish direction."
        )
        st.dataframe(read_tsv(str(PATHS["go_representative"])), width="stretch", hide_index=True)
        download_table("Download complete combined GO table", PATHS["go"])

elif section == "Methods explained":
    st.title("Methods explained")
    st.write("Start with the intuition, then connect it to the code and statistical output.")
    for concept, explanation in METHODS.items():
        with st.expander(concept):
            st.write(explanation)

elif section == "Results and interpretation":
    st.title("Results and interpretation")
    result_rows = [
        (
            "Paired differential expression",
            "6,315 genes have padj<0.05; 1,636 also have |shrunken log2FC|>=1.",
            "The result suggests broad transcript-level differences between Tumor and matched Normal tissue in this cohort.",
            "Bulk-cell composition, tissue structure, batch and other unmeasured factors may contribute. The result would need independent replication and experimental follow-up.",
        ),
        (
            "Hallmark enrichment",
            "35 of 50 Hallmark sets have padj<0.05; E2F, MYC and G2M rank toward Tumor.",
            "Coordinated proliferation and cell-cycle programs are enriched toward the Tumor-higher side of the ranked list.",
            "Gene-set overlap and tissue composition can produce coordinated patterns. Enrichment does not prove pathway activation and needs independent confirmation.",
        ),
        (
            "Directional GO BP enrichment",
            "259 Tumor-higher and 694 Normal-higher GO terms have BH FDR<0.05; both use the same tested background.",
            "Tumor-higher genes are over-represented in cell-division/chromosome themes, while Normal-higher genes include circulation, muscle, extracellular-matrix and tissue-context themes.",
            "GO terms overlap, annotation is incomplete and bulk-cell composition may explain part of either direction. These are hypotheses rather than demonstrated mechanisms.",
        ),
    ]
    for title, finding, interpretation, limitation in result_rows:
        st.subheader(title)
        st.markdown(f"**Finding:** {finding}")
        st.markdown(f"**Interpretation:** {interpretation}")
        st.markdown(f"**Important limitation:** {limitation}")
        st.divider()

elif section == "Limitations and next steps":
    st.title("Limitations and next steps")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("What this project can support")
        for item in CAN_SUPPORT:
            st.markdown(f"- {item}")
    with col2:
        st.subheader("What this project cannot support")
        for item in CANNOT_SUPPORT:
            st.markdown(f"- {item}")
    st.subheader("Important methodological limits")
    st.write(
        "Bulk tissue mixes cell populations; only one public dataset is analyzed; identifier and GO annotation "
        "mapping are incomplete; thresholds affect reported counts; and no independent cohort, clinical study, "
        "protein assay or functional experiment is included."
    )
    st.subheader("What would strengthen the analysis")
    st.write(
        "The next scientific step would be to repeat the paired model in an independent cohort, examine technical "
        "and cell-composition effects, and test selected genes or pathways with orthogonal measurements."
    )
    st.subheader("Project development and authorship")
    st.write(
        "OncoRNA is a student-led learning and portfolio project developed with AI-assisted coding, review and "
        "documentation. I selected and directed the project question and workflow, reviewed the generated outputs "
        "and validation reports, and am using the project to build deeper skills in RNA-seq, statistics, R, Python "
        "and reproducible bioinformatics. I remain responsible for explaining the methods, results and limitations."
    )
    st.subheader("References and project resources")
    st.markdown(
        """
        - NCBI Gene Expression Omnibus: [GSE306117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306117)
        - Love MI, Huber W, Anders S. [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8). *Genome Biology* (2014).
        - [fgsea Bioconductor package](https://bioconductor.org/packages/fgsea)
        - [clusterProfiler Bioconductor package](https://bioconductor.org/packages/clusterProfiler)
        - [MSigDB Hallmark gene sets](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H)
        - [Gene Ontology](https://geneontology.org/)
        """
    )

elif section == "Viva mode":
    st.title("Viva mode")
    st.write(
        "These are short explanations I can practise before discussing the project. Each answer can be expanded "
        "by connecting it to a figure, table or script."
    )
    for question, answer in VIVA_QA:
        with st.expander(question):
            st.write(answer)

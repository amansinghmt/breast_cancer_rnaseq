from __future__ import annotations

from datetime import datetime
from pathlib import Path

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
    "metrics": ROOT / "results_v2/robustness/analysis_metrics_v2.tsv",
    "thresholds": ROOT / "results_v2/robustness/de_threshold_sensitivity_v2.tsv",
    "prefilter": ROOT / "results_v2/robustness/low_count_prefilter_sensitivity_v2.tsv",
    "lfc_agreement": ROOT / "results_v2/robustness/lfc_agreement_v2.tsv",
    "pca_outliers": ROOT / "results_v2/robustness/pca_outlier_summary_v2.tsv",
    "cohort_summary": ROOT / "results_v2/robustness/cohort_inclusion_summary_v2.tsv",
    "top_de": ROOT / "results_v2/robustness/top_de_genes_v2.tsv",
    "output_manifest": ROOT / "results_v2/output_manifest.tsv",
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


def render_figure(figure_id: str) -> None:
    content = FIGURE_CONTENT[figure_id]
    st.subheader(f"{figure_id}. {content['title']}")
    st.image(str(FIGURE_PATHS[figure_id]), use_container_width=True)
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
      :root { --ink:#16302b; --teal:#187f69; --orange:#d65a1f; --paper:#f7f3ea; }
      .stApp { background: linear-gradient(145deg, #fbfaf6 0%, #f2f7f3 58%, #fff8ee 100%); color:var(--ink); }
      h1, h2, h3 { font-family: Georgia, 'Times New Roman', serif; color:var(--ink); letter-spacing:-0.02em; }
      [data-testid="stSidebar"] { background:#112e28; }
      [data-testid="stSidebar"] * { color:#f7f3ea; }
      div[data-testid="stMetric"] { background:rgba(255,255,255,.82); border:1px solid #dce6df; border-radius:14px; padding:14px; box-shadow:0 8px 28px rgba(22,48,43,.06); }
      div[data-testid="stExpander"] { background:rgba(255,255,255,.72); border-radius:12px; border-color:#dce6df; }
      .hero { padding:28px 34px; border-radius:22px; background:linear-gradient(105deg,#133f35,#187f69 68%,#d65a1f); color:white; margin-bottom:22px; }
      .hero h1 { color:white; margin:0; font-size:3.2rem; }
      .hero p { max-width:850px; font-size:1.08rem; margin:.7rem 0 0; color:#eef8f4; }
      .boundary { border-left:5px solid #d65a1f; background:#fff7ee; padding:15px 18px; border-radius:8px; }
      .small-source { color:#5b6d67; font-size:.85rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

missing = validate_required_files()
if missing:
    st.error(
        "The OncoRNA presentation cannot start because required canonical outputs are missing:\n\n"
        + "\n".join(f"- {path.relative_to(ROOT)}" for path in missing)
        + "\n\nRun `bash scripts/run_v2.sh` from the repository root."
    )
    st.stop()

metrics = metric_map()
manifest = read_tsv(str(PATHS["manifest"]))

sections = [
    "Overview",
    "Pipeline",
    "Dataset and cohort",
    "Quality control",
    "Differential expression",
    "Pathway analysis",
    "Methods explained",
    "Results and interpretation",
    "Limitations and claim boundaries",
    "Viva mode",
]

st.sidebar.title("OncoRNA")
st.sidebar.caption("Scientific presentation layer")
section = st.sidebar.radio("Navigate", sections)
st.sidebar.divider()
st.sidebar.caption(f"Latest pipeline log: {latest_analysis_date()}")
st.sidebar.caption("Canonical workflow: v2 paired")

if section == "Overview":
    st.markdown(
        """
        <div class="hero">
          <h1>OncoRNA</h1>
          <p>A reproducible paired bulk RNA-seq analysis asking which genes and biological programs differ between breast Tumor and matched Normal tissue in GEO dataset GSE306117.</p>
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
        ("padj + effect rule", "genes_padj_lt_0.05_abs_shrunk_lfc_ge_1"),
    ]
    for col, (label, key) in zip(cols, overview_metrics, strict=True):
        col.metric(label, f"{int(metrics[key]):,}")
    st.subheader("Why this project exists")
    st.write(
        "Gene-expression differences can reveal biological programs associated with tissue state. "
        "Matched samples are especially useful because each patient provides their own Normal baseline. "
        "The project connects gene-level statistical tests to pathway-level hypotheses while preserving "
        "reproducible inputs, outputs, software versions and checksums."
    )
    st.info(
        "Primary model: `~ patient_id + condition_main`. The patient term controls baseline "
        "between-person differences; the condition term estimates the paired Tumor-versus-Normal contrast."
    )
    st.subheader("Verified current result counts")
    a, b, c = st.columns(3)
    a.metric("Hallmark sets with padj < 0.05", f"{int(metrics['hallmark_sets_padj_lt_0.05']):,}")
    b.metric("GO terms with padj < 0.05", f"{int(metrics['go_terms_padj_lt_0.05']):,}")
    c.metric("Representative GO terms", f"{int(metrics['go_representative_terms']):,}")
    st.caption(f"Current presentation generated from the latest successful run: {latest_analysis_date()}.")

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
        use_container_width=True,
    )
    for name, operation, output, failure in PIPELINE_STAGES:
        with st.expander(name):
            st.markdown(f"**Operation:** {operation}")
            st.markdown(f"**Output:** `{output}`")
            st.markdown(f"**Possible failure:** {failure}")

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
    st.dataframe(cohort_summary, use_container_width=True, hide_index=True)
    included = manifest[manifest["include_paired"].astype(str).str.upper() == "TRUE"]
    st.subheader("Included matched cohort")
    st.dataframe(
        included[["sample_id", "patient_id", "condition_main", "library_size"]]
        .sort_values(["patient_id", "condition_main"]),
        use_container_width=True,
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
        st.dataframe(pca_outliers.head(15), use_container_width=True, hide_index=True)
        st.metric("Flagged samples", f"{int(metrics['pca_exploratory_outliers'])}")

elif section == "Differential expression":
    st.title("Differential expression")
    strict = read_tsv(str(PATHS["thresholds"]))
    primary = strict[(strict["padj_cutoff"] == 0.05) & (strict["abs_shrunken_log2fc_cutoff"] == 1)].iloc[0]
    cols = st.columns(4)
    cols[0].metric("Genes tested", f"{int(metrics['genes_tested']):,}")
    cols[1].metric("padj < 0.05", f"{int(metrics['genes_padj_lt_0.05']):,}")
    cols[2].metric("Tumor higher", f"{int(primary['upregulated']):,}")
    cols[3].metric("Normal higher", f"{int(primary['downregulated']):,}")
    with st.expander("Threshold and robustness summaries"):
        st.dataframe(strict, use_container_width=True, hide_index=True)
        st.markdown("**Predefined low-count sensitivity**")
        st.dataframe(read_tsv(str(PATHS["prefilter"])), use_container_width=True, hide_index=True)
        st.markdown("**Unshrunk versus shrunken effect agreement**")
        st.dataframe(read_tsv(str(PATHS["lfc_agreement"])), use_container_width=True, hide_index=True)
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
    st.dataframe(top_de, use_container_width=True, hide_index=True, height=430)
    download_table("Download canonical DE table", PATHS["de"])

elif section == "Pathway analysis":
    st.title("Pathway analysis")
    cols = st.columns(3)
    cols[0].metric("Hallmark sets tested", f"{int(metrics['hallmark_sets_tested'])}")
    cols[1].metric("Significant Hallmark sets", f"{int(metrics['hallmark_sets_padj_lt_0.05'])}")
    cols[2].metric("Significant GO BP terms", f"{int(metrics['go_terms_padj_lt_0.05'])}")
    tabs = st.tabs(["F06 Hallmark", "F07 GO BP"])
    with tabs[0]:
        render_figure("F06")
        hallmark = read_tsv(str(PATHS["hallmark"]))
        st.dataframe(hallmark, use_container_width=True, hide_index=True, height=350)
        download_table("Download Hallmark table", PATHS["hallmark"])
    with tabs[1]:
        render_figure("F07")
        representative = read_tsv(str(PATHS["go_representative"]))
        st.dataframe(representative, use_container_width=True, hide_index=True, height=350)
        col1, col2 = st.columns(2)
        with col1:
            download_table("Download representative GO table", PATHS["go_representative"])
        with col2:
            download_table("Download full GO test table", PATHS["go"])
    st.info(
        "The 729 significant terms are a subset of 5,102 tested GO terms. Earlier output contained "
        "only cutoff-passing terms, which made 729/729 appear significant. The current raw table "
        "preserves all tested terms; semantic reduction is presentation-only."
    )

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
            "Tumor and matched Normal tissues differ broadly at transcript level.",
            "Cell-type composition, tissue structure, batch and unmeasured covariates may contribute.",
            "Replicate in an independent matched cohort and validate selected genes experimentally.",
        ),
        (
            "Hallmark enrichment",
            "35 of 50 Hallmark sets have padj<0.05; E2F, MYC and G2M rank toward Tumor.",
            "Coordinated proliferation/cell-cycle-associated programs are consistent with the Tumor-side ranking.",
            "Gene-set overlap and tissue composition can produce coordinated rank patterns.",
            "Confirm with independent data and orthogonal proliferation measurements.",
        ),
        (
            "GO BP enrichment",
            "729 of 5,102 tested terms have padj<0.05; representative terms emphasize mitosis and chromosome segregation.",
            "The effect-filtered gene list is enriched for cell-division-associated annotations.",
            "GO terms are highly redundant and annotation coverage is incomplete.",
            "Use independent data and targeted functional assays before mechanistic claims.",
        ),
    ]
    for title, evidence, meaning, alternative, validation in result_rows:
        st.subheader(title)
        st.markdown(f"**Observation and evidence:** {evidence}")
        st.markdown(f"**Possible meaning:** {meaning}")
        st.markdown(f"**Alternative explanation:** {alternative}")
        st.markdown(f"**Required validation:** {validation}")
        st.divider()

elif section == "Limitations and claim boundaries":
    st.title("Limitations and claim boundaries")
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
    st.subheader("AI-use disclosure")
    st.write(
        "AI assistance contributed to code construction, software review, visual presentation and learning material. "
        "The analysis owner must personally understand and defend the cohort rules, statistical model, contrast, "
        "thresholds, outputs, limitations and claim boundaries."
    )

elif section == "Viva mode":
    st.title("Viva mode")
    st.write("Short answers first. Expand each answer by connecting it to a figure or code file.")
    for question, answer in VIVA_QA:
        with st.expander(question):
            st.write(answer)

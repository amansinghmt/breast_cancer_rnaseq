#!/usr/bin/env python3

import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing or empty: {path.relative_to(ROOT)}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    manifest_path = ROOT / "results_v2/metadata/paired_manifest.tsv"
    de_path = ROOT / "results_v2/deseq2/deseq2_paired_v2_results.tsv"
    hallmark_path = ROOT / "results_v2/enrichment/hallmark_gsea_paired_v2.tsv"
    go_path = ROOT / "results_v2/enrichment/go_bp_ora_paired_v2.tsv"
    go_tumor_path = ROOT / "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv"
    go_normal_path = ROOT / "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv"
    enrichment_diagnostics_path = ROOT / "results_v2/enrichment/enrichment_diagnostics_v2.tsv"
    metrics_path = ROOT / "results_v2/robustness/analysis_metrics_v2.tsv"

    manifest = read_tsv(manifest_path)
    included = [row for row in manifest if row["include_paired"].upper() == "TRUE"]
    if len(manifest) != 74 or len(included) != 42:
        raise AssertionError("Manifest row/inclusion count mismatch.")

    pairs: dict[str, list[str]] = {}
    for row in included:
        pairs.setdefault(row["patient_id"], []).append(row["condition_main"])
    if len(pairs) != 21:
        raise AssertionError("Expected 21 included patients.")
    if any(sorted(conditions) != ["Normal", "Tumor"] for conditions in pairs.values()):
        raise AssertionError("Every patient must have one Normal and one Tumor sample.")

    de = read_tsv(de_path)
    required_de = {
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "log2FoldChange_shrunk",
    }
    if not de or not required_de.issubset(de[0]):
        raise AssertionError("DE result columns are incomplete.")
    if len(de) != 60617:
        raise AssertionError(f"Expected 60,617 genes, found {len(de):,}.")

    hallmark = read_tsv(hallmark_path)
    go = read_tsv(go_path)
    go_tumor = read_tsv(go_tumor_path)
    go_normal = read_tsv(go_normal_path)
    if not {"pathway", "NES", "padj"}.issubset(hallmark[0]):
        raise AssertionError("Hallmark result columns are incomplete.")
    if not {"ID", "Description", "GeneRatio", "p.adjust"}.issubset(go[0]):
        raise AssertionError("GO result columns are incomplete.")
    for direction, table in (("Tumor-higher", go_tumor), ("Normal-higher", go_normal)):
        if not table or {row["analysis_direction"] for row in table} != {direction}:
            raise AssertionError(f"GO direction labels are invalid for {direction}.")
        bg_denominators = {row["BgRatio"].split("/")[-1] for row in table}
        if bg_denominators != {"15233"}:
            raise AssertionError(f"GO annotated universe mismatch for {direction}: {bg_denominators}")

    diagnostics = {
        row["metric"]: row["value"] for row in read_tsv(enrichment_diagnostics_path)
    }
    expected_go_diagnostics = {
        "ora_input_universe_genes": "30244",
        "ora_go_annotated_universe_genes": "15233",
        "ora_tumor_higher_input_genes": "506",
        "ora_tumor_higher_mapped_genes": "324",
        "ora_tumor_higher_terms_tested": "3095",
        "ora_tumor_higher_terms_padj_lt_0.05": "259",
        "ora_normal_higher_input_genes": "1130",
        "ora_normal_higher_mapped_genes": "712",
        "ora_normal_higher_terms_tested": "4686",
        "ora_normal_higher_terms_padj_lt_0.05": "694",
        "ora_multiple_testing_method": "Benjamini-Hochberg",
    }
    for key, value in expected_go_diagnostics.items():
        if diagnostics.get(key) != value:
            raise AssertionError(f"Directional GO diagnostic {key}: expected {value}, found {diagnostics.get(key)}")

    metrics = {row["metric"]: float(row["value"]) for row in read_tsv(metrics_path)}
    expected = {
        "manifest_rows": 74,
        "included_samples": 42,
        "paired_patients": 21,
        "genes_tested": 60617,
        "genes_padj_lt_0.05": 6315,
        "genes_padj_lt_0.05_abs_shrunk_lfc_ge_1": 1636,
        "hallmark_sets_padj_lt_0.05": 35,
        "go_terms_padj_lt_0.05": 729,
        "go_terms_tested": 5102,
        "go_tumor_higher_terms_tested": 3095,
        "go_tumor_higher_terms_padj_lt_0.05": 259,
        "go_normal_higher_terms_tested": 4686,
        "go_normal_higher_terms_padj_lt_0.05": 694,
    }
    for key, value in expected.items():
        if metrics.get(key) != value:
            raise AssertionError(f"Metric {key}: expected {value}, found {metrics.get(key)}")

    figure_dir = ROOT / "figures_v2/final"
    actual_figures = sorted(path.name for path in figure_dir.iterdir() if path.is_file())
    expected_figures = [f"F{i:02d}.png" for i in range(1, 8)]
    if actual_figures != expected_figures:
        raise AssertionError(f"Final figure contract mismatch: {actual_figures}")
    for figure in expected_figures:
        require(figure_dir / figure)

    vector_figures = sorted((ROOT / "figures_v2/vector").glob("F0*.pdf"))
    if len(vector_figures) != 7:
        raise AssertionError(f"Expected 7 vector PDFs, found {len(vector_figures)}.")
    for figure in vector_figures:
        require(figure)

    output_manifest = read_tsv(ROOT / "results_v2/output_manifest.tsv")
    for row in output_manifest:
        path = ROOT / row["path"]
        require(path)
        if md5(path) != row["md5"]:
            raise AssertionError(f"Checksum mismatch: {row['path']}")

    documentation = [
        ROOT / "docs/ONCORNA_PROJECT_UNDERSTANDING.md",
        ROOT / "docs/ONCORNA_VIVA_SHEET.md",
        ROOT / "docs/ONCORNA_RESULTS_AND_ROBUSTNESS.md",
        ROOT / "docs/ONCORNA_CODE_LEARNING_MAP.md",
        ROOT / "docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md",
        ROOT / "docs/ONCORNA_CLAIM_TRACEABILITY.md",
        ROOT / "docs/ONCORNA_MSC_PORTFOLIO_SUMMARY.md",
        ROOT / "docs/FINAL_REFINEMENT_AUDIT.md",
    ]
    for path in documentation:
        require(path)

    traceability = (ROOT / "docs/ONCORNA_CLAIM_TRACEABILITY.md").read_text()
    for figure_id in [f"F{i:02d}" for i in range(1, 8)]:
        if figure_id not in traceability:
            raise AssertionError(f"Claim traceability is missing {figure_id}.")

    report = (ROOT / "docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md").read_text()
    for number in ("6,315", "1,636", "506", "1,130", "259", "694", "30,244", "15,233"):
        if number not in report:
            raise AssertionError(f"Final report is missing headline number {number}.")
    if "combined" not in report.lower() or "non-directional" not in report.lower():
        raise AssertionError("Final report does not bound the combined GO interpretation.")

    f07_source = (ROOT / "scripts/04-figures/07_go_bp_dotplot_paired_v2.R").read_text()
    for label in ("Tumor-higher gene set", "Normal-higher gene set"):
        if label not in f07_source:
            raise AssertionError(f"F07 source is missing direction label: {label}")
    if 'go_bp_ora_paired_v2.tsv' in f07_source:
        raise AssertionError("F07 must not read the combined GO table.")

    presentation_files = [
        ROOT / "dashboard/content.py",
        ROOT / "docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md",
        ROOT / "docs/ONCORNA_VIVA_SHEET.md",
        *sorted((ROOT / "scripts/04-figures").glob("*.R")),
    ]
    presentation_text = "\n".join(path.read_text() for path in presentation_files)
    if re.search(r"total reads", presentation_text, flags=re.IGNORECASE):
        raise AssertionError("Presentation language must use retained gene-level counts, not total reads.")

    maintained_figure_text = "\n".join(
        path.read_text() for path in sorted((ROOT / "scripts/04-figures").glob("*.R"))
    )
    for color in ('"Normal" = "#1B9E77"', '"Tumor" = "#D95F02"'):
        if color not in maintained_figure_text:
            raise AssertionError(f"Condition color contract missing: {color}")
    if re.search(r'"Up in (tumor|normal)"', maintained_figure_text, flags=re.IGNORECASE):
        raise AssertionError("Maintained figures must use Tumor-higher/Normal-higher terminology.")

    print("ONCORNA OUTPUT VALIDATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

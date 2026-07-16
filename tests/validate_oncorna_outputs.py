#!/usr/bin/env python3

import csv
import hashlib
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
    if not {"pathway", "NES", "padj"}.issubset(hallmark[0]):
        raise AssertionError("Hallmark result columns are incomplete.")
    if not {"ID", "Description", "GeneRatio", "p.adjust"}.issubset(go[0]):
        raise AssertionError("GO result columns are incomplete.")

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
    ]
    for path in documentation:
        require(path)

    print("ONCORNA OUTPUT VALIDATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

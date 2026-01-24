#!/usr/bin/env python3

import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

COUNTS_PATH = Path("data/processed/counts.tsv")
METADATA_PATH = Path("data/metadata/metadata.tsv")
LIBRARY_PLOT_PATH = Path("figures/qc/library_size.png")
PCA_PLOT_PATH = Path("figures/qc/pca.png")
QC_SUMMARY_PATH = Path("results/qc/qc_summary.tsv")

REQUIRED_METADATA_COLUMNS = {
    "sample_id",
    "condition",
    "tp53_status",
    "lfs_status",
    "use_main_tumor_vs_normal",
}

CONDITION_SET = {"Tumor", "Normal"}
MIN_TUMOR = 10
MIN_NORMAL = 5


def load_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing metadata file: {path}")
    meta = pd.read_csv(path, sep="\t")
    missing_cols = REQUIRED_METADATA_COLUMNS - set(meta.columns)
    if missing_cols:
        raise ValueError(
            "Metadata missing required columns: " + ", ".join(sorted(missing_cols))
        )
    meta = meta[meta["use_main_tumor_vs_normal"] == "yes"].copy()
    meta = meta[meta["condition"].isin(CONDITION_SET)].copy()
    return meta


def read_counts_header(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing counts file: {path}")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip("\n")
    if not header:
        raise ValueError(f"Counts file is empty: {path}")
    columns = header.split("\t")
    if columns[0] != "gene_id":
        raise ValueError("Counts file must have 'gene_id' as the first column.")
    return columns[1:]


def load_counts(path: Path, sample_ids):
    usecols = ["gene_id"] + list(sample_ids)
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    df = df.set_index("gene_id")
    df = df.loc[:, sample_ids]
    return df


def compute_log2_cpm(counts: pd.DataFrame, library_sizes: pd.Series) -> pd.DataFrame:
    cpm = counts.divide(library_sizes, axis=1) * 1e6
    return np.log2(cpm + 1)


def plot_library_sizes(summary: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = summary["condition"].map({"Tumor": "#D55E00", "Normal": "#0072B2"})
    ax.bar(summary["sample_id"], summary["library_size"], color=colors)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Library size")
    ax.set_title("Library size per sample")
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_pca(scores: pd.DataFrame, explained: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))
    condition_colors = {"Tumor": "#D55E00", "Normal": "#0072B2"}

    unique_statuses = list(scores["tp53_status"].dropna().unique())
    marker_pool = ["o", "s", "^", "D", "v", "+", "x"]
    marker_map = {}
    for status in unique_statuses:
        if status == "WT":
            marker_map[status] = "o"
        elif status == "Mutant":
            marker_map[status] = "s"
        elif marker_pool:
            marker_map[status] = marker_pool.pop(0)
        else:
            marker_map[status] = "o"

    for (_, row) in scores.iterrows():
        color = condition_colors.get(row["condition"], "#444444")
        marker = marker_map.get(row["tp53_status"], "o")
        ax.scatter(row["PC1"], row["PC2"], color=color, marker=marker, s=50, alpha=0.8)

    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)")
    ax.set_title("PCA of log2(CPM + 1)")

    condition_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=label, markerfacecolor=color, markersize=8)
        for label, color in condition_colors.items()
    ]
    status_handles = [
        plt.Line2D([0], [0], marker=marker, color="#555555", label=status, linestyle="None")
        for status, marker in marker_map.items()
    ]

    legend1 = ax.legend(handles=condition_handles, title="Condition", loc="upper right")
    ax.add_artist(legend1)
    ax.legend(handles=status_handles, title="TP53 status", loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> int:
    try:
        metadata = load_metadata(METADATA_PATH)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    duplicate_metadata_ids = sorted(
        metadata.loc[metadata["sample_id"].duplicated(), "sample_id"].unique().tolist()
    )
    if duplicate_metadata_ids:
        metadata = metadata.drop_duplicates(subset=["sample_id"], keep="first")

    try:
        count_samples = read_counts_header(COUNTS_PATH)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    duplicate_count_ids = sorted(
        [sample_id for sample_id, count in Counter(count_samples).items() if count > 1]
    )

    metadata_samples = metadata["sample_id"].tolist()
    count_sample_set = set(count_samples)

    missing_in_counts = sorted(set(metadata_samples) - count_sample_set)
    extra_in_counts = sorted(count_sample_set - set(metadata_samples))

    sample_ids_used = [sample_id for sample_id in metadata_samples if sample_id in count_sample_set]
    metadata_used = metadata[metadata["sample_id"].isin(sample_ids_used)].copy()

    tumor_count = int((metadata_used["condition"] == "Tumor").sum())
    normal_count = int((metadata_used["condition"] == "Normal").sum())

    if tumor_count < MIN_TUMOR or normal_count < MIN_NORMAL:
        print(
            "Error: insufficient samples after filtering (Tumor: {tumor}, Normal: {normal})".format(
                tumor=tumor_count, normal=normal_count
            ),
            file=sys.stderr,
        )
        return 1

    if not sample_ids_used:
        print("Error: no matching samples between metadata and counts.", file=sys.stderr)
        return 1

    counts = load_counts(COUNTS_PATH, sample_ids_used)

    library_sizes = counts.sum(axis=0)

    qc_summary = metadata_used[["sample_id", "condition", "tp53_status", "lfs_status"]].copy()
    qc_summary["library_size"] = qc_summary["sample_id"].map(library_sizes)
    qc_summary = qc_summary[["sample_id", "library_size", "condition", "tp53_status", "lfs_status"]]

    QC_SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    qc_summary.to_csv(QC_SUMMARY_PATH, sep="\t", index=False)

    plot_library_sizes(qc_summary, LIBRARY_PLOT_PATH)

    log2_cpm = compute_log2_cpm(counts, library_sizes)
    variances = log2_cpm.var(axis=1)
    top_n = min(2000, variances.shape[0])
    top_genes = variances.sort_values(ascending=False).head(top_n).index
    pca_input = log2_cpm.loc[top_genes].T

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(pca_input.values)

    scores = metadata_used.set_index("sample_id").loc[pca_input.index].copy()
    scores["PC1"] = pcs[:, 0]
    scores["PC2"] = pcs[:, 1]

    explained = pca.explained_variance_ratio_ * 100
    plot_pca(scores, explained, PCA_PLOT_PATH)

    print(f"genes in counts matrix: {counts.shape[0]}")
    print(f"samples after filtering: {counts.shape[1]}")
    print(f"Tumor count: {tumor_count}")
    print(f"Normal count: {normal_count}")
    if missing_in_counts:
        print("missing in counts: " + ", ".join(missing_in_counts))
    else:
        print("missing in counts: none")
    if extra_in_counts:
        print("missing in metadata: " + ", ".join(extra_in_counts))
    else:
        print("missing in metadata: none")
    if duplicate_metadata_ids:
        print("duplicated in metadata: " + ", ".join(duplicate_metadata_ids))
    else:
        print("duplicated in metadata: none")
    if duplicate_count_ids:
        print("duplicated in counts: " + ", ".join(duplicate_count_ids))
    else:
        print("duplicated in counts: none")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

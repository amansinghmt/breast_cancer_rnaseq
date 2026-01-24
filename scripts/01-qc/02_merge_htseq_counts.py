#!/usr/bin/env python3

import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

INPUT_DIR = Path("data/raw/GSE306117_RAW")
OUTPUT_PATH = Path("data/processed/counts.tsv")
MIN_SAMPLES = 70
SAMPLE_ID_RE = re.compile(r"(GSM\d+)")


def list_count_files(input_dir: Path):
    if not input_dir.exists():
        raise FileNotFoundError(f"Missing input directory: {input_dir}")
    files = sorted(input_dir.glob("*.txt.gz"))
    if not files:
        raise ValueError(f"No .txt.gz files found in {input_dir}")
    return files


def extract_sample_id(filename: str):
    match = SAMPLE_ID_RE.search(filename)
    if not match:
        return None
    return match.group(1)


def load_counts(path: Path, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["gene_id", sample_id],
        compression="gzip",
        dtype={"gene_id": str},
    )
    df = df[~df["gene_id"].str.startswith("__", na=False)]
    return df.set_index("gene_id")


def main() -> int:
    try:
        files = list_count_files(INPUT_DIR)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    sample_entries = []
    missing_id_files = []
    for path in files:
        sample_id = extract_sample_id(path.name)
        if not sample_id:
            missing_id_files.append(path.name)
            continue
        sample_entries.append((sample_id, path))

    if missing_id_files:
        print(
            "Error: could not extract GSM ID from files: "
            + ", ".join(missing_id_files),
            file=sys.stderr,
        )
        return 1

    sample_ids = [sample_id for sample_id, _ in sample_entries]
    counts = Counter(sample_ids)
    duplicate_ids = sorted([sample_id for sample_id, count in counts.items() if count > 1])
    if duplicate_ids:
        print(
            "Error: duplicate GSM IDs found: " + ", ".join(duplicate_ids),
            file=sys.stderr,
        )
        return 1

    if len(sample_ids) < MIN_SAMPLES:
        print(
            f"Error: expected at least {MIN_SAMPLES} samples, found {len(sample_ids)}",
            file=sys.stderr,
        )
        return 1

    sample_entries.sort(key=lambda item: item[0])

    frames = []
    for sample_id, path in sample_entries:
        frames.append(load_counts(path, sample_id))

    merged = pd.concat(frames, axis=1, join="outer")
    merged = merged.fillna(0)
    merged = merged.astype(int)
    merged = merged.sort_index()

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    merged.reset_index().to_csv(OUTPUT_PATH, sep="\t", index=False)

    print(f"genes: {merged.shape[0]}")
    print(f"samples: {merged.shape[1]}")
    print("missing sample IDs: none")
    print("duplicated sample IDs: none")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

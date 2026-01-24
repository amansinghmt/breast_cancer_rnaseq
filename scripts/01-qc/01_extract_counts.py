#!/usr/bin/env python3

import sys
from io import StringIO
from pathlib import Path

import pandas as pd

BEGIN_MARKER = "!series_matrix_table_begin"
END_MARKER = "!series_matrix_table_end"


def extract_table_lines(path: Path):
    lines = []
    in_table = False
    found_begin = False
    found_end = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == BEGIN_MARKER:
                in_table = True
                found_begin = True
                continue
            if line == END_MARKER:
                found_end = True
                break
            if in_table:
                lines.append(line)

    if not found_begin or not found_end:
        missing = []
        if not found_begin:
            missing.append(BEGIN_MARKER)
        if not found_end:
            missing.append(END_MARKER)
        raise ValueError(f"Missing marker(s): {', '.join(missing)}")
    if not lines:
        raise ValueError("Expression table is empty between markers.")

    return lines


def load_counts(lines) -> pd.DataFrame:
    table_text = "\n".join(lines) + "\n"
    df = pd.read_csv(StringIO(table_text), sep="\t", header=0)
    if "ID_REF" not in df.columns:
        raise ValueError("Expected 'ID_REF' column in expression table.")
    return df.set_index("ID_REF")


def main() -> int:
    input_path = Path("data/raw/GSE306117_series_matrix.txt")
    output_path = Path("data/processed/counts.tsv")

    try:
        lines = extract_table_lines(input_path)
        df = load_counts(lines)
    except FileNotFoundError:
        print(f"Error: missing input file at {input_path}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    duplicate_ids = df.index[df.index.duplicated()].unique()
    duplicate_count = len(duplicate_ids)
    if duplicate_count:
        df = df.groupby(df.index).sum()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t")

    genes = df.shape[0]
    samples = df.shape[1]
    dup_label = "yes" if duplicate_count else "no"
    print(
        "genes: {genes}, samples: {samples}, duplicates handled: {dup_label} ({count})".format(
            genes=genes,
            samples=samples,
            dup_label=dup_label,
            count=duplicate_count,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

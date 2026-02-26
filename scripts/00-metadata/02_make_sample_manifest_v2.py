#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

METADATA_PATH = Path("data/metadata/metadata.tsv")
QC_SUMMARY_PATH = Path("results/qc/qc_summary.tsv")
OUTPUT_PATH = Path("data/metadata/sample_manifest.tsv")
QC_MIN = 1_000_000

REQUIRED_METADATA_COLS = {
    "sample_id",
    "patient_id",
    "tissue",
    "tissue_type",
    "condition_main",
    "lfs_status",
}

OUTPUT_COLUMNS = [
    "sample_id",
    "patient_id",
    "tissue_raw",
    "tissue_type",
    "condition_main",
    "lfs_status",
    "library_size",
    "include_paired",
    "exclude_reason",
]


def read_tsv(path: Path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        return rows, reader.fieldnames or []


def dedupe_reason_parts(reason_str: str):
    if not reason_str:
        return ""
    seen = set()
    parts = []
    for part in reason_str.split(";"):
        part = part.strip()
        if not part or part in seen:
            continue
        seen.add(part)
        parts.append(part)
    return ";".join(parts)


def main() -> int:
    if not METADATA_PATH.exists():
        raise SystemExit(f"ERROR: missing metadata file: {METADATA_PATH}")
    if not QC_SUMMARY_PATH.exists():
        raise SystemExit(f"ERROR: missing QC summary file: {QC_SUMMARY_PATH}")

    metadata_rows, metadata_cols = read_tsv(METADATA_PATH)
    missing_cols = sorted(REQUIRED_METADATA_COLS - set(metadata_cols))
    if missing_cols:
        raise SystemExit(
            "ERROR: metadata.tsv missing required columns: " + ", ".join(missing_cols)
        )

    qc_rows, qc_cols = read_tsv(QC_SUMMARY_PATH)
    if "sample_id" not in qc_cols or "library_size" not in qc_cols:
        raise SystemExit("ERROR: qc_summary.tsv missing sample_id/library_size columns")

    lib_by_sample = {}
    for row in qc_rows:
        sid = row["sample_id"]
        lib = row.get("library_size", "")
        if lib:
            lib_by_sample[sid] = int(lib)

    manifest_rows = []
    for row in metadata_rows:
        manifest_rows.append(
            {
                "sample_id": row["sample_id"],
                "patient_id": row["patient_id"],
                "tissue_raw": row["tissue"],
                "tissue_type": row["tissue_type"],
                "condition_main": row["condition_main"],
                "lfs_status": row["lfs_status"],
                "library_size": lib_by_sample.get(row["sample_id"]),
                "include_paired": "FALSE",
                "exclude_reason": "",
            }
        )

    for row in manifest_rows:
        reasons = []
        if row["tissue_type"] in {"Surrounding", "Unknown"} or row["condition_main"] == "NA":
            reasons.append("tissue_type_not_main")
        if row["library_size"] is None:
            reasons.append("library_size_missing")
        elif row["library_size"] < QC_MIN:
            reasons.append("low_library_size")
        row["exclude_reason"] = ";".join(reasons)

    by_patient = defaultdict(list)
    for row in manifest_rows:
        by_patient[row["patient_id"]].append(row)

    for _, patient_rows in by_patient.items():
        tumor_candidates = [
            r
            for r in patient_rows
            if r["condition_main"] == "Tumor"
            and r["library_size"] is not None
            and r["library_size"] >= QC_MIN
        ]
        normal_candidates = [
            r
            for r in patient_rows
            if r["condition_main"] == "Normal"
            and r["library_size"] is not None
            and r["library_size"] >= QC_MIN
        ]

        selected_tumor = (
            max(tumor_candidates, key=lambda r: r["library_size"])
            if tumor_candidates
            else None
        )
        selected_normal = (
            max(normal_candidates, key=lambda r: r["library_size"])
            if normal_candidates
            else None
        )

        if selected_tumor is not None and selected_normal is not None:
            selected_tumor["include_paired"] = "TRUE"
            selected_normal["include_paired"] = "TRUE"

        for row in patient_rows:
            if row["condition_main"] not in {"Tumor", "Normal"}:
                continue
            if row["library_size"] is None or row["library_size"] < QC_MIN:
                continue
            if row["include_paired"] == "TRUE":
                continue

            extras = []
            if row["condition_main"] == "Tumor":
                if selected_tumor is not None and row is not selected_tumor:
                    extras.append("duplicate_tumor_lower_library")
                if selected_normal is None:
                    extras.append("no_qc_passing_normal")
            elif row["condition_main"] == "Normal":
                if selected_normal is not None and row is not selected_normal:
                    extras.append("duplicate_normal_lower_library")
                if selected_tumor is None:
                    extras.append("no_qc_passing_tumor")

            if extras:
                merged = row["exclude_reason"]
                if merged:
                    merged = merged + ";" + ";".join(extras)
                else:
                    merged = ";".join(extras)
                row["exclude_reason"] = dedupe_reason_parts(merged)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            out = dict(row)
            out["library_size"] = (
                "NA" if out["library_size"] is None else str(out["library_size"])
            )
            out["exclude_reason"] = dedupe_reason_parts(out["exclude_reason"])
            writer.writerow(out)

    include_rows = [r for r in manifest_rows if r["include_paired"] == "TRUE"]
    include_patients = {r["patient_id"] for r in include_rows}
    print(f"Wrote {OUTPUT_PATH}")
    print(f"rows: {len(manifest_rows)}")
    print(f"include_paired rows: {len(include_rows)}")
    print(f"include_paired patients: {len(include_patients)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
import shlex
from pathlib import Path
import csv

INFILE = Path("data/raw/GSE306117_series_matrix.txt")
OUTFILE = Path("data/metadata/metadata.tsv")

def parse_geo_values(line: str):
    parts = shlex.split(line.strip())
    if len(parts) < 2:
        return None, []
    key = parts[0]
    if len(parts) >= 3 and parts[1] == "=":
        vals = parts[2:]
    else:
        vals = parts[1:]
    return key, vals

sample_id = []
title = []
source_name = []
tissue_vals = []
genotype_vals = []
batch_vals = []
char_lines = []

with INFILE.open() as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            _, sample_id = parse_geo_values(line)
        elif line.startswith("!Sample_title"):
            _, title = parse_geo_values(line)
        elif line.startswith("!Sample_source_name_ch1"):
            _, source_name = parse_geo_values(line)
        elif line.startswith("!Sample_characteristics_ch1"):
            char_lines.append(line)

for line in char_lines:
    _, vals = parse_geo_values(line)
    if not vals:
        continue
    first = vals[0].lower()
    if first.startswith("tissue:"):
        tissue_vals = vals
    elif first.startswith("genotype:"):
        genotype_vals = vals
    elif first.startswith("batch:"):
        batch_vals = vals

n = len(sample_id)
if n == 0:
    raise SystemExit("ERROR: Parsed 0 samples. '!Sample_geo_accession' line not captured.")

def require_len(name, arr):
    if len(arr) != n:
        raise SystemExit(f"ERROR: {name} has {len(arr)} values; expected {n}.")

require_len("title", title)
require_len("source_name", source_name)
require_len("tissue", tissue_vals)
require_len("genotype", genotype_vals)
require_len("batch", batch_vals)

def strip_prefix(x: str):
    return x.split(":", 1)[1].strip() if ":" in x else x.strip()

def tissue_type_from_tissue(t: str):
    tl = t.lower()
    # Priority matters: surrounding/adjacent must be checked before tumor.
    if "surrounding" in tl or "adjacent" in tl:
        return "Surrounding"
    if "control breast tissue" in tl:
        return "Normal"
    if "breast tumor" in tl:
        return "Tumor"
    return "Unknown"

def condition_from_tissue_type(tissue_type: str):
    if tissue_type in {"Tumor", "Normal", "Surrounding"}:
        return tissue_type
    return "Other"

def tp53_status_from_genotype(g: str):
    gl = g.lower()
    if "wt" in gl:
        return "WT"
    if "mutant" in gl:
        return "Mutant"
    return "Other"

def lfs_status_from_batch(b: str):
    bl = b.lower()
    if "non lfs" in bl:
        return "Non_LFS"
    if "lfs" in bl:
        return "LFS"
    return "Other"

OUTFILE.parent.mkdir(parents=True, exist_ok=True)

with OUTFILE.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "sample_id","title","source_name",
        "tissue","condition",
        "genotype","tp53_status",
        "batch","lfs_status",
        "use_main_tumor_vs_normal",
        "patient_id","tissue_type","condition_main"
    ])

    for i in range(n):
        tissue = strip_prefix(tissue_vals[i])
        genotype = strip_prefix(genotype_vals[i])
        batch = strip_prefix(batch_vals[i])
        patient_id = title[i].split("-", 1)[0].strip()

        tissue_type = tissue_type_from_tissue(tissue)
        condition = condition_from_tissue_type(tissue_type)
        condition_main = condition if condition in {"Tumor", "Normal"} else "NA"
        if tissue_type == "Surrounding" and condition_main != "NA":
            raise SystemExit(
                f"ERROR: Surrounding sample {sample_id[i]} has invalid condition_main={condition_main}"
            )
        tp53_status = tp53_status_from_genotype(genotype)
        lfs_status = lfs_status_from_batch(batch)
        use_main = "yes" if condition_main in {"Tumor","Normal"} else "no"

        w.writerow([
            sample_id[i], title[i], source_name[i],
            tissue, condition,
            genotype, tp53_status,
            batch, lfs_status,
            use_main,
            patient_id, tissue_type, condition_main
        ])

print(f"Wrote {OUTFILE} with {n} samples.")

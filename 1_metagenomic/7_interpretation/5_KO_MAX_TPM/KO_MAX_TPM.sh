#!/bin/bash
# Requires: Python 3 (pandas)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

export FILE="/path/to/base/1_metagenomic/7_interpretation/4_MAG_TPM/MAG_KO_TPM.tsv"
export OUT="/path/to/base/1_metagenomic/7_interpretation/5_KO_MAX_TPM"

mkdir -p "${OUT}"

python - <<'PY'
import pandas as pd
import os

infile = os.environ["FILE"]
outdir = os.environ["OUT"]

df = pd.read_csv(infile, sep="\t")

sample_cols = [
    "zr15207_9.sorted TPM",
    "zr15207_10.sorted TPM",
    "zr15207_11.sorted TPM",
    "zr15207_13.sorted TPM",
    "zr15207_14.sorted TPM"
]

# check that the columns exist (fail early if not)
missing = [c for c in sample_cols if c not in df.columns]
if missing:
    raise SystemExit(f"ERROR: Missing TPM columns: {missing}")

meta_cols = ["MAG", "KO", "Gene", "Description", "Pathway", "Pathway_Name"]
missing_meta = [c for c in meta_cols if c not in df.columns]
if missing_meta:
    raise SystemExit(f"ERROR: Missing meta columns: {missing_meta}")

os.makedirs(outdir, exist_ok=True)

for col in sample_cols:
    sample_name = col.split(".sorted")[0]  # zr15207_9
    tmp = df.sort_values(["KO", col], ascending=[True, False])
    top3 = tmp.groupby("KO", as_index=False).head(3)

    out = top3[meta_cols + [col]].rename(columns={col: "TPM"})

    outfile = os.path.join(outdir, f"{sample_name}_KO_Top3MAG.tsv")
    out.to_csv(outfile, sep="\t", index=False)

print(f"Done. Files saved in: {outdir}")
PY

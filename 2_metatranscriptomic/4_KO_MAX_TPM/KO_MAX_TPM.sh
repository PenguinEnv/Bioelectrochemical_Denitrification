#!/bin/bash
# Requires: Python 3 (pandas)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

export FILE="/path/to/base/2_metatranscriptomic/3_TPM/MAG_KO_TPM_with_annotation.csv"
export OUT="/path/to/base/2_metatranscriptomic/4_KO_MAX_TPM"

mkdir -p "${OUT}"

python - <<'PY'
import pandas as pd
import os

infile = os.environ["FILE"]
outdir = os.environ["OUT"]

df = pd.read_csv(infile)

sample_cols = ["zr15207_2T","zr15207_3T","zr15207_4T","zr15207_6T","zr15207_7T"]

missing = [c for c in sample_cols if c not in df.columns]
if missing:
    raise SystemExit(f"ERROR: Missing TPM columns: {missing}\nColumns seen: {list(df.columns)}")

meta_cols = ["MAG", "KO", "KO_Description", "KO_Pathways"]
missing_meta = [c for c in meta_cols if c not in df.columns]
if missing_meta:
    raise SystemExit(f"ERROR: Missing meta columns: {missing_meta}\nColumns seen: {list(df.columns)}")

os.makedirs(outdir, exist_ok=True)

for col in sample_cols:
    # ✅ sample name (zr15207_2T -> zr15207_2)
    sample_name = col.rstrip("T")

    tmp = df.sort_values(["KO", col], ascending=[True, False])
    top3 = tmp.groupby("KO", as_index=False).head(3)

    out = top3[meta_cols + [col]].rename(columns={col: "TPM"})

    outfile = os.path.join(outdir, f"{sample_name}_KO_Top3MAG.tsv")
    out.to_csv(outfile, sep="\t", index=False)

print(f"Done. Files saved in: {outdir}")
PY

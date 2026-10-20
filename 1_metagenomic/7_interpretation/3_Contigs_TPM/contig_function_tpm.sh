#!/bin/bash
# Requires: Python 3 (pandas, numpy); KEGG ko_list / ko_to_pathway / pathway_list files
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# === Paths ===
DEPTH_TSV="/path/to/base/1_metagenomic/4_coverage/coverm_contig_tpm.tsv"
DIAMOND="/path/to/base/1_metagenomic/7_interpretation/2_enrichm/annot_output/annotations_ko/DIAMOND_search.tsv"
KO_INFO="/path/to/your_db/KEGG"
OUTDIR="/path/to/base/1_metagenomic/7_interpretation/3_Contigs_TPM"

KO_TPM_CLEAN="$OUTDIR/KO_TPM.clean.tsv"
LOGFILE="$OUTDIR/pipeline.log"

mkdir -p "$OUTDIR"
echo "[$(date)] === Pipeline started ===" | tee -a "$LOGFILE"

export DEPTH_TSV DIAMOND KO_INFO OUTDIR

python3 - <<'PY' 2>&1 | tee -a "$LOGFILE"
import os
import pandas as pd, numpy as np, re
from pathlib import Path

def log(m): print(m, flush=True)

DEPTH_TSV = os.environ["DEPTH_TSV"]
DIAMOND   = os.environ["DIAMOND"]
KO_INFO   = os.environ["KO_INFO"]
OUTDIR    = os.environ["OUTDIR"]
KO_TPM_CLEAN = str(Path(OUTDIR, "KO_TPM.clean.tsv"))

log("=== Python step started ===")

# 1) Contig TPM (TPM columns only)
depth = pd.read_csv(DEPTH_TSV, sep="\t")
assert depth.columns[0] == "Contig", f"First column is not Contig: {depth.columns[0]}"
depth = depth.rename(columns={"Contig":"contig"})
tpm_cols = [c for c in depth.columns if c.endswith(" TPM")]
if not tpm_cols:
    raise RuntimeError("No TPM columns found (column names ending in ' TPM').")
contig_tpm = depth[["contig"] + tpm_cols].copy()

# 2) Load DIAMOND (as is)
dia = pd.read_csv(
    DIAMOND, sep="\t", header=None,
    names=["query","subject","identity","align_length","mismatch","gap",
           "q_start","q_end","s_start","s_end","evalue","bitscore"]
)
dia["align_length"] = pd.to_numeric(dia["align_length"], errors="coerce")
dia["bitscore"] = pd.to_numeric(dia["bitscore"], errors="coerce")

# 2-1) query -> contig / query_id (split on tilde)
# e.g. final.contigs_1000.part_031~k141_1043518_1
q = dia["query"].astype(str)
has_tilde = q.str.contains("~", regex=False)
token = q.copy()
token[has_tilde] = q[has_tilde].str.split("~").str[-1]  # k141_1043518_1

dia["query_id"] = token  # for best hit per query
# contig = k141_####
dia["contig"] = token.str.extract(r'^(k\d+_\d+)')

# 2-2) subject -> KO (K + 5 digits)
s = dia["subject"].astype(str)
ko1 = s.str.extract(r'~(K\d{5})')[0]     # UniRef...~Kxxxxx
ko2 = s.str.extract(r'ko:(K\d{5})')[0]   # ko:Kxxxxx
ko3 = s.str.extract(r'(K\d{5})')[0]      # fallback
dia["KO"] = ko1.fillna(ko2).fillna(ko3)

# keep rows with required values
dia = dia.dropna(subset=["contig","query_id","KO","align_length","bitscore"]).copy()

# 3) one best hit per query (max bitscore, then longest alignment)
dia = dia.sort_values(["query_id","bitscore","align_length"], ascending=[True, False, False])
best = dia.drop_duplicates(subset=["query_id"], keep="first").copy()

# 4) distribute 100% of contig TPM over its best-hit queries by alignment length
best = best.merge(contig_tpm, on="contig", how="left")
best = best.dropna(subset=tpm_cols)
best["denom"] = best.groupby("contig")["align_length"].transform("sum")
best = best[best["denom"] > 0].copy()
w = best["align_length"] / best["denom"]
for c in tpm_cols:
    best[c] = best[c] * w

# 5) sum TPM per KO x sample
ko_tpm = best.groupby("KO", as_index=False)[tpm_cols].sum()

# 6) KEGG KO descriptions (from ko_list)
desc = pd.read_csv(f"{KO_INFO}/ko_list.txt", sep="\t", header=None, names=["KO","Description"])
desc["KO"] = desc["KO"].str.replace("ko:","",regex=False)
desc["Gene"] = desc["Description"].str.extract(r'^([A-Za-z0-9_:-]+)')
desc["Description"] = desc["Description"].str.replace(r'^([A-Za-z0-9_:-]+)\s+','',regex=True)

# keep valid KOs only (drop malformed IDs)
valid = set(desc["KO"])
ko_tpm = ko_tpm[ko_tpm["KO"].isin(valid)].merge(desc[["KO","Gene","Description"]], on="KO", how="left")

# 6-1) KO -> pathway mapping
ko2path = pd.read_csv(f"{KO_INFO}/ko_to_pathway.txt", sep="\t", header=None, names=["ko","path"])
ko2path["KO"] = ko2path["ko"].str.replace("ko:","", regex=False)
ko2path["Pathway_ID"] = ko2path["path"].str.replace("path:","", regex=False)

path_df = pd.read_csv(f"{KO_INFO}/pathway_list.txt", sep="\t", header=None, names=["Pathway_ID","Pathway_Name"])
path_map = dict(zip(path_df["Pathway_ID"], path_df["Pathway_Name"]))

# pathway names per KO (joined with ",")
ko_to_pathnames = (
    ko2path
    .groupby("KO")["Pathway_ID"]
    .apply(lambda ids: ", ".join(path_map.get(i, i) for i in ids))
    .to_dict()
)

ko_tpm["Pathways"] = ko_tpm["KO"].map(ko_to_pathnames).fillna("")

# 7) strip sample-name suffix and save
rename = {c: c.replace(" TPM","") for c in tpm_cols}
out = ko_tpm.rename(columns=rename)

# column order: KO, Gene, Description, Pathways, samples
out = out[["KO","Gene","Description","Pathways"] + [rename[c] for c in tpm_cols]]

out.to_csv(KO_TPM_CLEAN, sep="\t", index=False)
log(f"✅ KO TPM saved → {KO_TPM_CLEAN}")

# 8) check totals (KO sum <= contig sum)
tot = contig_tpm[tpm_cols].sum()
kos = ko_tpm[tpm_cols].sum()
excess = (kos - tot).round(6)
if (excess > 1e-3).any():
    log("⚠️  WARNING: Some samples exceed total contig TPM. Check DIAMOND duplication or regex.")
else:
    log("✔︎ Validation passed: KO sums do not exceed total contig TPM (expected < 1e6).")

log("=== Python step completed ===")
PY

echo "[$(date)] ✅ Done: $KO_TPM_CLEAN" | tee -a "$LOGFILE"

#!/bin/bash
# Requires: Python 3 (pandas)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

################################################################################
# USER SETTINGS
# Only edit this section
################################################################################

INDIR="/path/to/base/1_metagenomic/7_interpretation/6_Gene_normal"

# Previous MAG-KO activity output
MAG_ACTIVITY_FILE="$INDIR/MAG_KO_activity_RNA_over_DNAcopy.tsv"

OUTDIR="/path/to/base/1_metagenomic/7_interpretation/6_Gene_normal/pathway_MAG_activity"
mkdir -p "$OUTDIR"

# Minimum DNA copy threshold for a MAG-KO to be included
MIN_DNA_COPY=0.1

# Conditions
CONDITIONS="2T 3T 4T 6T 7T"

# Exclude broad/unhelpful KEGG pathway names? 1=yes, 0=no
EXCLUDE_BROAD_PATHWAYS=1

################################################################################
# DO NOT EDIT BELOW THIS LINE
################################################################################

echo "============================================================"
echo "MAG KEGG pathway activity summary"
echo "============================================================"
echo "Input file   = $MAG_ACTIVITY_FILE"
echo "Output dir   = $OUTDIR"
echo "MIN_DNA_COPY = $MIN_DNA_COPY"
echo "CONDITIONS   = $CONDITIONS"
echo "============================================================"

if [[ ! -s "$MAG_ACTIVITY_FILE" ]]; then
    echo "ERROR: input file not found or empty: $MAG_ACTIVITY_FILE" >&2
    exit 1
fi

python3 - << PY
import os
import re
import pandas as pd

MAG_ACTIVITY_FILE = "$MAG_ACTIVITY_FILE"
OUTDIR = "$OUTDIR"
MIN_DNA_COPY = float("$MIN_DNA_COPY")
CONDITIONS = "$CONDITIONS".split()
EXCLUDE_BROAD_PATHWAYS = int("$EXCLUDE_BROAD_PATHWAYS")

cond_to_dna = {
    "2T": "zr15207_9.sorted_DNA_copy_per_MAG",
    "3T": "zr15207_10.sorted_DNA_copy_per_MAG",
    "4T": "zr15207_11.sorted_DNA_copy_per_MAG",
    "6T": "zr15207_13.sorted_DNA_copy_per_MAG",
    "7T": "zr15207_14.sorted_DNA_copy_per_MAG",
}

# Broad KEGG categories that often dominate but are not mechanistically specific
broad_pathways = {
    "Metabolic pathways",
    "Biosynthesis of secondary metabolites",
    "Microbial metabolism in diverse environments",
    "Biosynthesis of amino acids",
    "Carbon metabolism",
    "Biosynthesis of cofactors",
    "Two-component system",
    "ABC transporters",
}

def clean_pathway_name(x):
    x = str(x).strip()
    x = re.sub(r"^map\d{5}\s*", "", x)
    x = re.sub(r"^ko\d{5}\s*", "", x)
    x = re.sub(r"\s+", " ", x)
    return x.strip()

def split_pathways(x):
    if pd.isna(x):
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "nan":
        return []
    
    # KO_Pathways in your file is comma-separated:
    # "Glycolysis / Gluconeogenesis, ko00010, Fatty acid degradation, ko00071, ..."
    parts = [p.strip() for p in x.split(",") if p.strip()]
    
    # Remove ko/map IDs and keep pathway names only.
    names = []
    for p in parts:
        if re.fullmatch(r"(ko|map)\d{5}", p):
            continue
        if re.fullmatch(r"K\d{5}", p):
            continue
        name = clean_pathway_name(p)
        if name:
            names.append(name)
    
    # unique while preserving order
    seen = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out

print("Reading MAG-KO activity file...")
df = pd.read_csv(MAG_ACTIVITY_FILE, sep="\t")

required = {"MAG", "KO"}
if not required.issubset(df.columns):
    raise SystemExit(f"ERROR: input must contain {required}. Found: {list(df.columns)}")

# Detect pathway column
pathway_col = None
for cand in ["KO_Pathways", "Pathway", "Pathways", "KO_Pathway", "Pathway_Name"]:
    if cand in df.columns:
        pathway_col = cand
        break

if pathway_col is None:
    raise SystemExit(
        "ERROR: no pathway column found. Expected one of KO_Pathways, Pathway, Pathways, KO_Pathway, Pathway_Name."
    )

print("Using pathway column:", pathway_col)

# Explode KO_Pathways to long format
records = []
for idx, row in df.iterrows():
    paths = split_pathways(row[pathway_col])
    if EXCLUDE_BROAD_PATHWAYS:
        paths = [p for p in paths if p not in broad_pathways]
    if not paths:
        continue
    
    base = row.to_dict()
    for p in paths:
        rec = base.copy()
        rec["KEGG_Pathway"] = p
        records.append(rec)

long = pd.DataFrame(records)

if long.empty:
    raise SystemExit("ERROR: no pathway assignments after parsing/filtering.")

# Convert numeric columns
available = []
for cond in CONDITIONS:
    rna_col = f"zr15207_{cond}"
    dna_col = cond_to_dna.get(cond)
    ratio_col = f"{cond}_RNA_over_DNAcopy"
    
    if rna_col in long.columns and dna_col in long.columns:
        long[rna_col] = pd.to_numeric(long[rna_col], errors="coerce").fillna(0)
        long[dna_col] = pd.to_numeric(long[dna_col], errors="coerce").fillna(0)
        if ratio_col in long.columns:
            long[ratio_col] = pd.to_numeric(long[ratio_col], errors="coerce")
        available.append((cond, rna_col, dna_col, ratio_col))
    else:
        print("WARNING: missing columns for condition:", cond, rna_col, dna_col)

if not available:
    raise SystemExit("ERROR: no usable condition columns found.")

# Save exploded KO-pathway mapping
mapping_cols = ["MAG", "KO", "KEGG_Pathway"]
for c in ["KO_Description", "KO_Pathways"]:
    if c in long.columns:
        mapping_cols.append(c)

long[mapping_cols].drop_duplicates().to_csv(
    os.path.join(OUTDIR, "MAG_KO_to_KEGG_pathway.exploded.tsv"),
    sep="\t", index=False
)

# Main summary: MAG x Pathway x Condition
summary_records = []

for cond, rna_col, dna_col, ratio_col in available:
    sub = long[long[dna_col] >= MIN_DNA_COPY].copy()
    
    for (mag, pathway), g in sub.groupby(["MAG", "KEGG_Pathway"], dropna=False):
        # Avoid double-counting same KO assigned multiple times to same pathway within MAG
        g2 = g.drop_duplicates(subset=["MAG", "KO", "KEGG_Pathway"]).copy()
        
        sum_rna = g2[rna_col].sum()
        sum_dna = g2[dna_col].sum()
        weighted_activity = sum_rna / sum_dna if sum_dna > 0 else None
        
        ratio_values = None
        mean_ratio = None
        median_ratio = None
        sum_ratio = None
        
        if ratio_col in g2.columns:
            ratio_values = pd.to_numeric(g2[ratio_col], errors="coerce").dropna()
            if len(ratio_values) > 0:
                mean_ratio = ratio_values.mean()
                median_ratio = ratio_values.median()
                sum_ratio = ratio_values.sum()
        
        summary_records.append({
            "MAG": mag,
            "KEGG_Pathway": pathway,
            "Condition": cond,
            "n_KO_used": g2["KO"].nunique(),
            "KO_list_used": ";".join(sorted(g2["KO"].unique())),
            "sum_RNA_TPM": sum_rna,
            "sum_DNA_copy": sum_dna,
            "weighted_activity_sumRNA_over_sumDNA": weighted_activity,
            "mean_RNA_over_DNAcopy": mean_ratio,
            "median_RNA_over_DNAcopy": median_ratio,
            "sum_RNA_over_DNAcopy": sum_ratio
        })

summary = pd.DataFrame(summary_records)

if summary.empty:
    raise SystemExit("ERROR: summary table is empty after filtering. Try lowering MIN_DNA_COPY.")

out_long = os.path.join(OUTDIR, f"MAG_KEGG_pathway_activity.long_DNAcopyge{MIN_DNA_COPY}.tsv")
summary.to_csv(out_long, sep="\t", index=False)
print("Written:", out_long)

# Wide tables
metrics = [
    "sum_RNA_TPM",
    "sum_DNA_copy",
    "weighted_activity_sumRNA_over_sumDNA",
    "mean_RNA_over_DNAcopy",
    "median_RNA_over_DNAcopy",
    "n_KO_used"
]

for metric in metrics:
    wide = summary.pivot_table(
        index=["MAG", "KEGG_Pathway"],
        columns="Condition",
        values=metric,
        aggfunc="first"
    ).reset_index()
    wide.columns.name = None
    
    out = os.path.join(OUTDIR, f"MAG_KEGG_pathway_{metric}.wide_DNAcopyge{MIN_DNA_COPY}.tsv")
    wide.to_csv(out, sep="\t", index=False)
    print("Written:", out)

# Top pathways per MAG-condition
top_records = []
for (mag, cond), g in summary.groupby(["MAG", "Condition"]):
    gg = g.sort_values("weighted_activity_sumRNA_over_sumDNA", ascending=False).head(20).copy()
    gg["rank_within_MAG_condition"] = range(1, len(gg) + 1)
    top_records.append(gg)

top_mag_pathways = pd.concat(top_records, ignore_index=True)
out_top_mag = os.path.join(OUTDIR, f"Top20_pathways_by_MAG_condition_DNAcopyge{MIN_DNA_COPY}.tsv")
top_mag_pathways.to_csv(out_top_mag, sep="\t", index=False)
print("Written:", out_top_mag)

# Top MAGs per pathway-condition
top_records = []
for (pathway, cond), g in summary.groupby(["KEGG_Pathway", "Condition"]):
    gg = g.sort_values("weighted_activity_sumRNA_over_sumDNA", ascending=False).head(20).copy()
    gg["rank_within_pathway_condition"] = range(1, len(gg) + 1)
    top_records.append(gg)

top_mags = pd.concat(top_records, ignore_index=True)
out_top_mags = os.path.join(OUTDIR, f"Top20_MAGs_by_KEGG_pathway_activity_DNAcopyge{MIN_DNA_COPY}.tsv")
top_mags.to_csv(out_top_mags, sep="\t", index=False)
print("Written:", out_top_mags)

# Community summary from MAGs: all MAGs per pathway-condition
community_records = []
for cond, rna_col, dna_col, ratio_col in available:
    sub = long[long[dna_col] >= MIN_DNA_COPY].copy()
    
    for pathway, g in sub.groupby("KEGG_Pathway"):
        g2 = g.drop_duplicates(subset=["MAG", "KO", "KEGG_Pathway"]).copy()
        sum_rna = g2[rna_col].sum()
        sum_dna = g2[dna_col].sum()
        community_records.append({
            "KEGG_Pathway": pathway,
            "Condition": cond,
            "n_MAG_used": g2["MAG"].nunique(),
            "n_KO_used": g2["KO"].nunique(),
            "sum_RNA_TPM": sum_rna,
            "sum_DNA_copy": sum_dna,
            "weighted_activity_sumRNA_over_sumDNA": sum_rna / sum_dna if sum_dna > 0 else None,
            "KO_list_used": ";".join(sorted(g2["KO"].unique()))
        })

community = pd.DataFrame(community_records)
out_comm = os.path.join(OUTDIR, f"Community_from_MAG_KEGG_pathway_activity_DNAcopyge{MIN_DNA_COPY}.tsv")
community.to_csv(out_comm, sep="\t", index=False)
print("Written:", out_comm)

# Top community pathways per condition
top_records = []
for cond, g in community.groupby("Condition"):
    gg = g.sort_values("weighted_activity_sumRNA_over_sumDNA", ascending=False).head(30).copy()
    gg["rank_within_condition"] = range(1, len(gg) + 1)
    top_records.append(gg)

top_comm = pd.concat(top_records, ignore_index=True)
out_top_comm = os.path.join(OUTDIR, f"Top30_community_KEGG_pathways_by_activity_DNAcopyge{MIN_DNA_COPY}.tsv")
top_comm.to_csv(out_top_comm, sep="\t", index=False)
print("Written:", out_top_comm)

print("\\nDONE.")
print("Main files:")
print(out_long)
print(out_comm)
print(out_top_comm)
PY

echo ""
echo "Finished."
echo "Results are in:"
echo "$OUTDIR"

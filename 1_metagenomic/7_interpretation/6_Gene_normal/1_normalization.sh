#!/bin/bash
# Requires: Python 3 (pandas, numpy)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

################################################################################
# Only edit this section
################################################################################

PROJECT="/path/to/base"

# ---- GTDB-Tk marker summary ----
GTDB_MARKERS="$PROJECT/1_metagenomic/7_interpretation/1_gtdbtk/gtdbtk_out_drep/identify/gtdbtk.bac120.markers_summary.tsv"

# ---- DIAMOND KO annotation from EnrichM ----
DIAMOND="$PROJECT/1_metagenomic/7_interpretation/2_enrichm/annot_output/annotations_ko/DIAMOND_search.tsv"

# ---- Contig to MAG/bin table ----
CONTIG_BIN="$PROJECT/1_metagenomic/7_interpretation/4_MAG_TPM/contig_to_bin.tsv"

# ---- DNA coverage files ----
CONTIG_COV="$PROJECT/1_metagenomic/4_coverage/coverm_contig_mean.tsv"
MAG_COV="$PROJECT/1_metagenomic/4_coverage/coverm_genome_mean.tsv"

# ---- RNA TPM files ----
RNA_KO_TPM="$PROJECT/2_metatranscriptomic/3_TPM/KO_TPM_contig_level_with_annotation.csv"
RNA_MAG_KO_TPM="$PROJECT/2_metatranscriptomic/3_TPM/MAG_KO_TPM_with_annotation.csv"

# ---- Output directory ----
OUTDIR="$PROJECT/1_metagenomic/7_interpretation/6_Gene_normal"

# ---- Sample mapping ----
# DNA sample names should match columns in coverm_contig_mean.tsv / coverm_genome.tsv
# RNA sample names should match columns in KO_TPM_contig_level_with_annotation.csv
SAMPLE_MAP_CONTENT=$(cat << 'EOF'
DNA_sample	RNA_sample	Condition
zr15207_9.sorted	zr15207_2T	2T
zr15207_10.sorted	zr15207_3T	3T
zr15207_11.sorted	zr15207_4T	4T
zr15207_13.sorted	zr15207_6T	6T
zr15207_14.sorted	zr15207_7T	7T
EOF
)

# ---- Optional filters for final activity tables ----
# These do not remove rows from the main output, but create additional filtered tables.
MIN_COMMUNITY_DNA=0.01
MIN_MAG_DNA_COPY=0.10

################################################################################
# DO NOT EDIT BELOW THIS LINE
################################################################################

mkdir -p "$OUTDIR"

cat > "$OUTDIR/sample_map.tsv" <<< "$SAMPLE_MAP_CONTENT"

echo "============================================================"
echo "GTDB-marker DNA/RNA activity pipeline"
echo "============================================================"
echo "PROJECT        = $PROJECT"
echo "GTDB_MARKERS   = $GTDB_MARKERS"
echo "DIAMOND        = $DIAMOND"
echo "CONTIG_BIN     = $CONTIG_BIN"
echo "CONTIG_COV     = $CONTIG_COV"
echo "MAG_COV        = $MAG_COV"
echo "RNA_KO_TPM     = $RNA_KO_TPM"
echo "RNA_MAG_KO_TPM = $RNA_MAG_KO_TPM"
echo "OUTDIR         = $OUTDIR"
echo "============================================================"

# Check input files
for f in "$GTDB_MARKERS" "$DIAMOND" "$CONTIG_BIN" "$CONTIG_COV" "$MAG_COV" "$RNA_KO_TPM" "$RNA_MAG_KO_TPM"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: file not found or empty: $f" >&2
        exit 1
    fi
done

python3 - << PY
import os
import re
import numpy as np
import pandas as pd

################################################################################
# Paths from bash
################################################################################

GTDB_MARKERS = r"$GTDB_MARKERS"
DIAMOND = r"$DIAMOND"
CONTIG_BIN = r"$CONTIG_BIN"
CONTIG_COV = r"$CONTIG_COV"
MAG_COV = r"$MAG_COV"
RNA_KO_TPM = r"$RNA_KO_TPM"
RNA_MAG_KO_TPM = r"$RNA_MAG_KO_TPM"
OUTDIR = r"$OUTDIR"
SAMPLE_MAP = os.path.join(OUTDIR, "sample_map.tsv")

MIN_COMMUNITY_DNA = float("$MIN_COMMUNITY_DNA")
MIN_MAG_DNA_COPY = float("$MIN_MAG_DNA_COPY")

PSEUDO = 1e-9

################################################################################
# Helper functions
################################################################################

def norm_sample(c):
    """
    Normalize sample column names from CoverM outputs.
    Example:
    'zr15207_10.sorted TPM' -> 'zr15207_10.sorted'
    """
    c = str(c).strip()
    c = re.sub(r"\s+TPM$", "", c)
    c = re.sub(r"\s+Mean$", "", c)
    return c

def find_col(df, candidates):
    lower = {str(c).lower(): c for c in df.columns}
    for x in candidates:
        if x.lower() in lower:
            return lower[x.lower()]
    for c in df.columns:
        for x in candidates:
            if x.lower() in str(c).lower():
                return c
    raise ValueError(f"Cannot find any of {candidates} in columns: {list(df.columns)}")

def parse_marker_list(x):
    """
    Parse comma-separated GTDB-Tk marker list.
    """
    if pd.isna(x):
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "nan":
        return []
    parts = re.split(r",", x)
    clean = []
    for p in parts:
        p = p.strip()
        p = re.sub(r"\s+", "", p)
        if p:
            clean.append(p)
    return clean

def write(df, filename):
    path = os.path.join(OUTDIR, filename)
    df.to_csv(path, sep="\t", index=False)
    print(f"Written: {path}")
    return path

################################################################################
# Step 1. Parse GTDB-Tk bac120 marker summary
################################################################################

print("\n[1] Parsing GTDB-Tk bac120 unique single-copy markers")

markers = pd.read_csv(GTDB_MARKERS, sep="\t", dtype=str)
markers.columns = [c.strip() for c in markers.columns]

required_cols = {"name", "list_unique_genes"}
missing = required_cols - set(markers.columns)
if missing:
    raise SystemExit(f"GTDB marker summary missing columns: {missing}. Found: {list(markers.columns)}")

marker_rows = []
for _, row in markers.iterrows():
    mag = row["name"]
    unique_markers = parse_marker_list(row["list_unique_genes"])
    for m in unique_markers:
        marker_rows.append((mag, m, 1))

marker_presence = pd.DataFrame(marker_rows, columns=["MAG", "marker_id", "copy_number"]).drop_duplicates()

marker_count_by_mag = (
    marker_presence.groupby("MAG", as_index=False)
    .agg(n_unique_markers=("marker_id", "nunique"))
)

write(marker_presence, "GTDB_bac120_unique_marker_presence.tsv")
write(marker_count_by_mag, "GTDB_bac120_unique_marker_count_by_MAG.tsv")

print(f"  MAGs in GTDB marker summary: {markers['name'].nunique()}")
print(f"  Unique marker families: {marker_presence['marker_id'].nunique()}")
print(f"  Marker presence rows: {len(marker_presence)}")

################################################################################
# Step 2. Build contig_to_KO_MAG from DIAMOND + contig_to_bin
################################################################################

print("\n[2] Building contig_to_KO_MAG table from DIAMOND + contig_to_bin")

dia = pd.read_csv(
    DIAMOND,
    sep="\t",
    header=None,
    names=[
        "query", "subject", "identity", "align_length", "mismatch", "gap",
        "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"
    ]
)

dia["align_length"] = pd.to_numeric(dia["align_length"], errors="coerce")
dia["bitscore"] = pd.to_numeric(dia["bitscore"], errors="coerce")
dia["evalue"] = pd.to_numeric(dia["evalue"], errors="coerce")

# query -> query_id and contig
q = dia["query"].astype(str)
token = q.copy()
has_tilde = q.str.contains("~", regex=False)
token[has_tilde] = q[has_tilde].str.split("~").str[-1]

dia["query_id"] = token
dia["Contig"] = token.str.extract(r"^(k\d+_\d+)")

# subject -> KO
s = dia["subject"].astype(str)
ko1 = s.str.extract(r"~(K\d{5})")[0]
ko2 = s.str.extract(r"ko:(K\d{5})")[0]
ko3 = s.str.extract(r"(K\d{5})")[0]
dia["KO"] = ko1.fillna(ko2).fillna(ko3)

dia = dia.dropna(subset=["Contig", "KO", "bitscore"])

# Best hit per predicted gene/query
dia = dia.sort_values(
    ["query_id", "bitscore", "align_length"],
    ascending=[True, False, False]
)
best = dia.drop_duplicates("query_id", keep="first").copy()

# Read contig_to_bin.tsv
rows = []
with open(CONTIG_BIN) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith("Processing"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            rows.append((parts[0], parts[1]))

contig_bin = pd.DataFrame(rows, columns=["Contig", "MAG"]).drop_duplicates()

anno = best.merge(contig_bin, on="Contig", how="left")
rows_without_mag = anno["MAG"].isna().sum()
anno = anno.dropna(subset=["MAG"]).copy()

anno = anno[
    [
        "Contig", "MAG", "query_id", "KO",
        "identity", "align_length", "evalue", "bitscore", "subject"
    ]
].copy()

write(anno, "contig_to_KO_MAG.tsv")

print(f"  DIAMOND best-hit genes: {len(best)}")
print(f"  Rows without MAG: {rows_without_mag}")
print(f"  Final annotated MAG-gene rows: {len(anno)}")
print(f"  Unique MAGs: {anno['MAG'].nunique()}")
print(f"  Unique KOs: {anno['KO'].nunique()}")

################################################################################
# Step 3. Read DNA coverage tables
################################################################################

print("\n[3] Reading DNA coverage tables")

contig_cov = pd.read_csv(CONTIG_COV, sep="\t")
mag_cov = pd.read_csv(MAG_COV, sep="\t")

contig_cov.columns = [norm_sample(c) for c in contig_cov.columns]
mag_cov.columns = [norm_sample(c) for c in mag_cov.columns]

contig_col = find_col(contig_cov, ["Contig", "contig"])
mag_col = find_col(mag_cov, ["MAG", "Genome", "genome", "bin"])

contig_cov = contig_cov.rename(columns={contig_col: "Contig"})
mag_cov = mag_cov.rename(columns={mag_col: "MAG"})

for df in [contig_cov, mag_cov]:
    for c in df.columns:
        if c not in ["Contig", "MAG"]:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

smap = pd.read_csv(SAMPLE_MAP, sep="\t")
dna_samples = [x for x in smap["DNA_sample"] if x in contig_cov.columns]

if len(dna_samples) == 0:
    raise SystemExit("No DNA sample columns found in contig coverage table. Check sample_map.tsv.")

missing_mag_cov = [x for x in dna_samples if x not in mag_cov.columns]
if missing_mag_cov:
    raise SystemExit(f"These DNA samples are missing from MAG coverage table: {missing_mag_cov}")

print(f"  DNA samples used: {dna_samples}")

################################################################################
# Step 4. Community-level GTDB marker denominator
################################################################################

print("\n[4] Calculating community-level genome-equivalent coverage using GTDB markers")

# For each marker family, sum coverages of MAGs that contain the marker as unique/single-copy.
marker_cov = marker_presence.merge(mag_cov[["MAG"] + dna_samples], on="MAG", how="left")
marker_cov[dna_samples] = marker_cov[dna_samples].fillna(0)

marker_family_cov = (
    marker_cov.groupby("marker_id", as_index=False)[dna_samples]
    .sum()
)

write(marker_family_cov, "GTDB_marker_family_summed_MAG_coverage.tsv")

genome_equiv_rows = []
for sample in dna_samples:
    vals = marker_family_cov[sample].replace(0, np.nan).dropna()
    genome_equiv_rows.append({
        "DNA_sample": sample,
        "n_marker_families_used": int(vals.shape[0]),
        "marker_family_sum_mean": float(vals.mean()) if vals.shape[0] else np.nan,
        "marker_family_sum_median": float(vals.median()) if vals.shape[0] else np.nan,
        "marker_family_sum_min": float(vals.min()) if vals.shape[0] else np.nan,
        "marker_family_sum_max": float(vals.max()) if vals.shape[0] else np.nan
    })

genome_equiv = pd.DataFrame(genome_equiv_rows)
write(genome_equiv, "GTDB_marker_genome_equivalent_coverage.tsv")

print(genome_equiv.to_string(index=False))

ge_median = dict(zip(genome_equiv["DNA_sample"], genome_equiv["marker_family_sum_median"]))

################################################################################
# Step 5. Community-level KO DNA abundance
################################################################################

print("\n[5] Calculating community-level KO DNA abundance")

ann_cov = anno.merge(contig_cov, on="Contig", how="inner")

# Important:
# Multiple predicted genes on the same contig can have the same KO.
# Deduplicate Contig-KO to avoid counting the same contig coverage repeatedly for the same KO.
ann_cov_dedup = ann_cov.drop_duplicates(subset=["Contig", "KO"]).copy()

community_ko_cov = (
    ann_cov_dedup.groupby("KO", as_index=False)[dna_samples]
    .sum()
)

write(community_ko_cov, "community_KO_DNA_raw_contig_coverage.tsv")

community_ko_dna = community_ko_cov.copy()
for sample in dna_samples:
    community_ko_dna[sample + "_DNA_cov_per_GTDB_marker"] = (
        community_ko_dna[sample] / (ge_median[sample] + PSEUDO)
    )

community_dna_cols = ["KO"] + [s + "_DNA_cov_per_GTDB_marker" for s in dna_samples]
write(community_ko_dna[community_dna_cols], "community_KO_DNA_cov_per_GTDB_marker.tsv")

print(f"  Community KO rows: {len(community_ko_dna)}")

################################################################################
# Step 6. Community-level RNA/DNA activity
################################################################################

print("\n[6] Calculating community-level RNA/DNA activity")

rna_ko = pd.read_csv(RNA_KO_TPM)

if "KO" not in rna_ko.columns:
    raise SystemExit("RNA KO TPM file must contain a 'KO' column.")

community_activity = rna_ko.merge(
    community_ko_dna[community_dna_cols],
    on="KO",
    how="left"
)

activity_cols = []
for _, row in smap.iterrows():
    d = row["DNA_sample"]
    r = row["RNA_sample"]
    cond = row["Condition"]
    dcol = d + "_DNA_cov_per_GTDB_marker"
    outcol = cond + "_RNA_over_DNA_GTDBmarker"
    if r in community_activity.columns and dcol in community_activity.columns:
        community_activity[outcol] = community_activity[r] / (community_activity[dcol] + PSEUDO)
        activity_cols.append(outcol)

write(community_activity, "community_KO_activity_RNA_over_DNA_GTDB_marker.tsv")

# Filtered activity table to avoid extreme ratios from near-zero DNA abundance
filter_mask = np.ones(len(community_activity), dtype=bool)
for _, row in smap.iterrows():
    d = row["DNA_sample"]
    dcol = d + "_DNA_cov_per_GTDB_marker"
    if dcol in community_activity.columns:
        filter_mask = filter_mask & (
            community_activity[dcol].fillna(0).astype(float) >= MIN_COMMUNITY_DNA
        )

community_activity_filtered = community_activity.loc[filter_mask].copy()
write(
    community_activity_filtered,
    f"community_KO_activity_RNA_over_DNA_GTDB_marker.filtered_DNAge{MIN_COMMUNITY_DNA}.tsv"
)

################################################################################
# Step 7. MAG-level KO DNA copy estimate
################################################################################

print("\n[7] Calculating MAG-level KO DNA copy estimates")

# MAG-level:
# KO copy per MAG = coverage of KO-carrying contigs within MAG / MAG coverage.
# Deduplicate MAG-Contig-KO to avoid overcounting same contig for same KO in same MAG.
ann_cov_mag_dedup = ann_cov.drop_duplicates(subset=["MAG", "Contig", "KO"]).copy()

magko_cov = (
    ann_cov_mag_dedup.groupby(["MAG", "KO"], as_index=False)[dna_samples]
    .sum()
)

write(magko_cov, "MAG_KO_DNA_raw_contig_coverage.tsv")

magko_copy = magko_cov.copy()
for sample in dna_samples:
    mc = mag_cov[["MAG", sample]].rename(columns={sample: sample + "_MAG_coverage"})
    magko_copy = magko_copy.merge(mc, on="MAG", how="left")
    magko_copy[sample + "_DNA_copy_per_MAG"] = (
        magko_copy[sample] / (magko_copy[sample + "_MAG_coverage"] + PSEUDO)
    )

write(magko_copy, "MAG_KO_DNA_copy_per_MAG.tsv")

print(f"  MAG-KO rows: {len(magko_copy)}")

################################################################################
# Step 8. MAG-level RNA/DNA activity
################################################################################

print("\n[8] Calculating MAG-level RNA/DNA activity")

rna_magko = pd.read_csv(RNA_MAG_KO_TPM)

if not {"MAG", "KO"}.issubset(set(rna_magko.columns)):
    raise SystemExit("RNA MAG-KO TPM file must contain 'MAG' and 'KO' columns.")

mag_copy_cols = ["MAG", "KO"] + [s + "_DNA_copy_per_MAG" for s in dna_samples]

magko_activity = rna_magko.merge(
    magko_copy[mag_copy_cols],
    on=["MAG", "KO"],
    how="left"
)

mag_activity_cols = []
for _, row in smap.iterrows():
    d = row["DNA_sample"]
    r = row["RNA_sample"]
    cond = row["Condition"]
    dcol = d + "_DNA_copy_per_MAG"
    outcol = cond + "_RNA_over_DNAcopy"
    if r in magko_activity.columns and dcol in magko_activity.columns:
        magko_activity[outcol] = magko_activity[r] / (magko_activity[dcol] + PSEUDO)
        mag_activity_cols.append(outcol)

write(magko_activity, "MAG_KO_activity_RNA_over_DNAcopy.tsv")

# Filtered MAG-level activity table
filter_mask = np.ones(len(magko_activity), dtype=bool)
for _, row in smap.iterrows():
    d = row["DNA_sample"]
    dcol = d + "_DNA_copy_per_MAG"
    if dcol in magko_activity.columns:
        filter_mask = filter_mask & (
            magko_activity[dcol].fillna(0).astype(float) >= MIN_MAG_DNA_COPY
        )

magko_activity_filtered = magko_activity.loc[filter_mask].copy()
write(
    magko_activity_filtered,
    f"MAG_KO_activity_RNA_over_DNAcopy.filtered_DNAcopyge{MIN_MAG_DNA_COPY}.tsv"
)

################################################################################
# Step 9. Summary
################################################################################

summary = {
    "n_MAGs_in_GTDB_marker_summary": int(markers["name"].nunique()),
    "n_GTDB_marker_families_detected": int(marker_presence["marker_id"].nunique()),
    "n_DIAMOND_best_hit_genes": int(len(best)),
    "n_annotated_gene_rows_with_MAG": int(len(anno)),
    "n_unique_MAGs_in_annotation": int(anno["MAG"].nunique()),
    "n_unique_KOs_in_annotation": int(anno["KO"].nunique()),
    "n_community_KO_rows": int(len(community_ko_dna)),
    "n_MAG_KO_rows": int(len(magko_copy)),
    "min_community_DNA_filter": MIN_COMMUNITY_DNA,
    "min_MAG_DNA_copy_filter": MIN_MAG_DNA_COPY,
}

summary_df = pd.DataFrame([summary])
write(summary_df, "pipeline_summary.tsv")

print("\n============================================================")
print("DONE")
print("============================================================")
print(f"Output directory: {OUTDIR}")
print("")
print("Key output files:")
for fname in [
    "GTDB_marker_genome_equivalent_coverage.tsv",
    "community_KO_DNA_cov_per_GTDB_marker.tsv",
    "community_KO_activity_RNA_over_DNA_GTDB_marker.tsv",
    f"community_KO_activity_RNA_over_DNA_GTDB_marker.filtered_DNAge{MIN_COMMUNITY_DNA}.tsv",
    "MAG_KO_DNA_copy_per_MAG.tsv",
    "MAG_KO_activity_RNA_over_DNAcopy.tsv",
    f"MAG_KO_activity_RNA_over_DNAcopy.filtered_DNAcopyge{MIN_MAG_DNA_COPY}.tsv",
    "pipeline_summary.tsv",
]:
    print(" - " + os.path.join(OUTDIR, fname))
PY

echo ""
echo "Pipeline finished."
echo "Results are in:"
echo "$OUTDIR"

#!/bin/bash
# Requires: Python 3 (pandas); KEGG ko_list / ko_to_pathway / pathway_list files
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# === Output folder ===
OUTDIR="/path/to/base/2_metatranscriptomic/3_TPM"
mkdir -p "$OUTDIR"
LOGFILE="$OUTDIR/pipeline.log"

echo "[$(date)] === MAG & Contig KO TPM pipeline started ===" | tee -a "$LOGFILE"

python3 - <<'PY' 2>&1 | tee -a "$LOGFILE"
import os
import pandas as pd
from glob import glob

# ==========================
# 🔧 0. Paths
# ==========================

# (1) Salmon quant results (contig level)
#   -> parent folder of directories like zr15207_2T_salmon_quant/quant.sf
contig_quant_dir = "/path/to/base/2_metatranscriptomic/2_mapping"

# (2) EnrichM annotation results (contigs >= 1000 bp)
base_annot_dir = "/path/to/base/1_metagenomic/7_interpretation/2_enrichm/annot_output"
diamond_tsv    = os.path.join(base_annot_dir, "annotations_ko",  "DIAMOND_search.tsv")
gff_dir        = os.path.join(base_annot_dir, "annotations_gff")

# (3) contig -> MAG mapping file
ctb_path = "/path/to/base/1_metagenomic/7_interpretation/4_MAG_TPM/contig_to_bin.tsv"

# (4) KO / pathway annotation files
ko_desc_file     = "/path/to/your_db/KEGG/ko_list.txt"
ko_pathway_file  = "/path/to/your_db/KEGG/ko_to_pathway.txt"
pathway_desc_file= "/path/to/your_db/KEGG/pathway_list.txt"

# (5) output folder
OUTDIR = "/path/to/base/2_metatranscriptomic/3_TPM"
os.makedirs(OUTDIR, exist_ok=True)

print("▶ DIAMOND:", diamond_tsv)
print("▶ GFF dir:", gff_dir)
print("▶ contig_to_bin:", ctb_path)
print("▶ Salmon quant dir:", contig_quant_dir)
print("▶ Output dir:", OUTDIR)

# ==========================
# 1. DIAMOND: CDS -> KO, contig mapping
# ==========================
print("\n[1] DIAMOND: reading CDS-contig-KO mapping...")

diamond = pd.read_csv(
    diamond_tsv,
    sep="\t",
    header=None,
    usecols=[0, 1],
    names=["CDS_raw", "Target"]
)

# e.g. final.contigs.part_014~k141_103172_1  -> CDS_ID = k141_103172_1
diamond["CDS_ID"] = diamond["CDS_raw"].str.extract(r"~(k141_\d+_\d+)")
# contig: k141_103172
diamond["Contig"] = diamond["CDS_ID"].str.extract(r"^(k141_\d+)")
# KO: Kxxxx
diamond["KO"] = diamond["Target"].str.extract(r"~(K\d+)")

mapping_df = diamond.dropna(subset=["CDS_ID", "Contig", "KO"]).copy()
print(f"  - total DIAMOND hits: {len(diamond):,}")
print(f"  - valid mappings (with CDS_ID+contig+KO): {len(mapping_df):,}")

# ==========================
# 2. GFF: CDS lengths
# ==========================
print("\n[2] reading CDS lengths from GFF...")

gff_files = sorted(glob(os.path.join(gff_dir, "final.contigs.part_*.gff")))
if not gff_files:
    raise FileNotFoundError(f"GFF files not found: {gff_dir}")

cds_rows = []
for gff_path in gff_files:
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            fields = line.strip().split("\t")
            contig = fields[0]
            start = int(fields[3])
            end   = int(fields[4])
            length = abs(end - start) + 1
            attr = fields[8]
            # seq_id=k141_34544_1;...
            if "seq_id=" not in attr:
                continue
            cds_id = attr.split("seq_id=")[1].split(";")[0]
            cds_rows.append((contig, cds_id, length))

cds_len_df = pd.DataFrame(cds_rows, columns=["Contig", "CDS_ID", "CDS_Length"])
print(f"  - CDS read from GFF: {len(cds_len_df):,}")

# join CDS lengths to the DIAMOND mapping
mapping_df = mapping_df.merge(
    cds_len_df,
    on=["Contig", "CDS_ID"],
    how="inner"
)
mapping_df = mapping_df.dropna(subset=["CDS_Length"])

print(f"  - mapping rows with CDS length: {len(mapping_df):,}")
if mapping_df.empty:
    raise RuntimeError("mapping_df is empty. Check that GFF and DIAMOND contig/CDS IDs match.")

# ==========================
# 3. contig_to_bin: contig -> MAG mapping
# ==========================
print("\n[3] reading contig_to_bin (contig -> MAG)...")

ctb = pd.read_csv(
    ctb_path,
    sep="\t",
    header=None,
    names=["Contig", "MAG"]
)

# filter out lines such as "Processing ..."
ctb = ctb[ctb["Contig"].str.startswith("k141_")].copy()
print(f"  - contig_to_bin mappings: {len(ctb):,}")

# add MAG to mapping_df (by contig)
mapping_df = mapping_df.merge(ctb, on="Contig", how="left")
mag_mapped = mapping_df["MAG"].notna().sum()
print(f"  - CDS assigned to a MAG: {mag_mapped:,}")

# ==========================
# 4. Salmon quant: distribute contig TPM to CDS
#     - KO × Sample (community-level)
#     - MAG × KO × Sample (MAG-level)
# ==========================
print("\n[4] reading Salmon quant.sf and distributing TPM...")

quant_files = sorted(glob(os.path.join(contig_quant_dir, "zr15207_*_salmon_quant", "quant.sf")))
if not quant_files:
    raise FileNotFoundError(f"Salmon quant.sf files not found: {contig_quant_dir}")

print(f"  - quant.sf files found: {len(quant_files)}")

sample_ko_tpm = {}       # KO × Sample
sample_mag_ko_tpm = {}   # MAG × KO × Sample

for qf in quant_files:
    sample = os.path.basename(os.path.dirname(qf)).replace("_salmon_quant", "")
    print(f"    ▶ processing sample: {sample} ({qf})")

    quant_df = pd.read_csv(qf, sep="\t")
    # Name -> Contig; Salmon keeps the final.contigs.fa headers (k141_xxx)
    quant_df = quant_df.rename(columns={"Name": "Contig"})

    # join with mapping_df by contig
    merged = mapping_df.merge(quant_df[["Contig", "TPM"]], on="Contig", how="inner")
    if merged.empty:
        print(f"      - warning: no contigs mapped to KO/MAG in {sample} (merged empty)")
        continue

    # distribute TPM by CDS length within each contig
    ko_rows = []
    mag_ko_rows = []

    for contig, g in merged.groupby("Contig"):
        total_len = g["CDS_Length"].sum()
        if total_len <= 0:
            continue
        contig_tpm = g["TPM"].iloc[0]  # same contig, so same TPM
        frac = g["CDS_Length"] / total_len
        tpm_cds = contig_tpm * frac

        for (_, row), tpm_val in zip(g.iterrows(), tpm_cds):
            ko = row["KO"]
            ko_rows.append((ko, tpm_val))

            mag = row["MAG"]
            if isinstance(mag, str):  # only when a MAG is assigned
                mag_ko_rows.append((mag, ko, tpm_val))

    if not ko_rows:
        print(f"      - warning: no KO-level TPM in {sample}")
    else:
        df_ko = pd.DataFrame(ko_rows, columns=["KO", "TPM"])
        sample_ko_tpm[sample] = df_ko.groupby("KO")["TPM"].sum()

    if not mag_ko_rows:
        print(f"      - note: no CDS in MAGs for {sample} (no MAG x KO TPM)")
    else:
        df_mag_ko = pd.DataFrame(mag_ko_rows, columns=["MAG", "KO", "TPM"])
        sample_mag_ko_tpm[sample] = df_mag_ko.groupby(["MAG", "KO"])["TPM"].sum()

print(f"\n  - samples with KO-level TPM: {len(sample_ko_tpm)}")
print(f"  - samples with MAG x KO TPM: {len(sample_mag_ko_tpm)}")
if not sample_ko_tpm:
    raise RuntimeError("No KO-level TPM was computed for any sample.")

# ==========================
# 5. KO x sample TPM matrix (community level)
# ==========================
print("\n[5] building KO x sample TPM matrix...")

combined_ko = pd.concat(sample_ko_tpm, names=["Sample"])
long_ko = combined_ko.reset_index()   # Sample, KO, TPM

ko_matrix = long_ko.pivot_table(
    index="KO",
    columns="Sample",
    values="TPM",
    fill_value=0
).reset_index()

# ==========================
# 6. MAG x KO x sample TPM matrix (MAG level)
# ==========================
if sample_mag_ko_tpm:
    print("\n[6] building MAG x KO x sample TPM matrix...")

    combined_mag = pd.concat(sample_mag_ko_tpm, names=["Sample"])
    long_mag = combined_mag.reset_index()  # Sample, MAG, KO, TPM

    mag_matrix = long_mag.pivot_table(
        index=["MAG", "KO"],
        columns="Sample",
        values="TPM",
        fill_value=0
    ).reset_index()
else:
    print("  - no MAG x KO TPM data; MAG matrix not created.")
    mag_matrix = None

# ==========================
# 7. add KO descriptions and pathway annotation
# ==========================
print("\n[7] adding KO descriptions and pathways...")

ko_desc_df = pd.read_csv(ko_desc_file, sep="\t", header=None, names=["KO", "Description"])
ko_desc_dict = dict(zip(ko_desc_df["KO"], ko_desc_df["Description"]))

ko_pathway_df = pd.read_csv(ko_pathway_file, sep="\t", header=None, names=["ko", "pathway"])
ko_pathway_df["KO"] = ko_pathway_df["ko"].str.replace("ko:", "", regex=False)
ko_pathway_df["Pathway_ID"] = ko_pathway_df["pathway"].str.replace("path:", "", regex=False)
ko_to_pathways = ko_pathway_df.groupby("KO")["Pathway_ID"].apply(list).to_dict()

pathway_df = pd.read_csv(pathway_desc_file, sep="\t", header=None, names=["Pathway_ID", "Pathway_Name"])
pathway_dict = dict(zip(pathway_df["Pathway_ID"], pathway_df["Pathway_Name"]))

def get_pathways(ko):
    ids = ko_to_pathways.get(ko, [])
    return ", ".join([pathway_dict.get(pid, pid) for pid in ids]) if ids else ""

# --- annotate the KO-level matrix ---
ko_matrix["KO_Description"] = ko_matrix["KO"].map(ko_desc_dict).fillna("Unknown KO")
ko_matrix["KO_Pathways"]    = ko_matrix["KO"].map(get_pathways)

# tidy columns
base_cols_ko = ["KO", "KO_Description", "KO_Pathways"]
sample_cols_ko = [c for c in ko_matrix.columns if c not in base_cols_ko]
sample_cols_ko = sorted(sample_cols_ko)

ko_final = ko_matrix[base_cols_ko + sample_cols_ko]

# --- annotate the MAG-level matrix (if present) ---
if mag_matrix is not None:
    mag_matrix["KO_Description"] = mag_matrix["KO"].map(ko_desc_dict).fillna("Unknown KO")
    mag_matrix["KO_Pathways"]    = mag_matrix["KO"].map(get_pathways)

    base_cols_mag = ["MAG", "KO", "KO_Description", "KO_Pathways"]
    sample_cols_mag = [c for c in mag_matrix.columns if c not in base_cols_mag]
    sample_cols_mag = sorted(sample_cols_mag)

    mag_final = mag_matrix[base_cols_mag + sample_cols_mag]
else:
    mag_final = None

# ==========================
# 8. save
# ==========================
print("\n[8] writing CSV...")

ko_out = os.path.join(OUTDIR, "KO_TPM_contig_level_with_annotation.csv")
ko_final.to_csv(ko_out, index=False)
print(f"  - KO TPM (contig level, community) saved -> {ko_out}")

if mag_final is not None:
    mag_out = os.path.join(OUTDIR, "MAG_KO_TPM_with_annotation.csv")
    mag_final.to_csv(mag_out, index=False)
    print(f"  - MAG x KO TPM (MAG level) saved -> {mag_out}")
else:
    print("  - MAG x KO TPM file not created (no MAG mapping).")

print("\n✅ All done!")

PY

echo "[$(date)] ✅ Done: MAG & Contig KO TPM pipeline" | tee -a "$LOGFILE"

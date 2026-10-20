
#!/bin/bash
# Requires: Python 3 (pandas); KEGG ko_list / ko_to_pathway / pathway_list files
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# === 📁 Paths ===
DREP_DIR="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"
DEPTH_TSV="/path/to/base/1_metagenomic/4_coverage/coverm_contig_tpm.tsv"
DIAMOND="/path/to/base/1_metagenomic/7_interpretation/2_enrichm/annot_output/annotations_ko/DIAMOND_search.tsv"
KO_INFO="/path/to/your_db/KEGG"
OUTDIR="/path/to/base/1_metagenomic/7_interpretation/4_MAG_TPM"
CTB="$OUTDIR/contig_to_bin.tsv"
MAG_ABUNDANCE="$OUTDIR/MAG_abundance_by_sample.tsv"
MAG_KO_TPM="$OUTDIR/MAG_KO_TPM.tsv"
LOGFILE="$OUTDIR/pipeline.log"

mkdir -p "$OUTDIR"

echo "=== Pipeline started ===" | tee -a "$LOGFILE"

# 🔄 Step 1: contig-to-bin mapping from the dereplicated MAGs
echo "Generating contig-to-bin mapping..." | tee -a "$LOGFILE"
for f in "$DREP_DIR"/*.fa; do
    echo "Processing $f" | tee -a "$LOGFILE"
    grep ">" "$f" | sed 's/>//' | awk -v bin=$(basename "$f" .fa) '{print $1 "\t" bin}'
done > "$CTB"
echo "Done: $CTB" | tee -a "$LOGFILE"

# 🧮 Step 2: MAG abundance and KO TPM in Python
echo "Starting Python TPM calculation..." | tee -a "$LOGFILE"
python3 - <<EOF 2>&1 | tee -a "$LOGFILE"
import pandas as pd
import sys

def log(msg):
    print(msg, flush=True)

log("=== Python step started ===")

# === Bash variables -> Python variables ===
DEPTH_TSV = "$DEPTH_TSV"
DIAMOND = "$DIAMOND"
KO_INFO = "$KO_INFO"
CTB = "$CTB"
MAG_ABUNDANCE = "$MAG_ABUNDANCE"
MAG_KO_TPM = "$MAG_KO_TPM"
KO_TPM = "$OUTDIR/KO_TPM.tsv"

# 1) MAG abundance
log("Calculating MAG abundance per sample...")
depth = pd.read_csv(DEPTH_TSV, sep="\t")
depth = depth.rename(columns={"Contig": "contig"})
log(f"Depth columns: {list(depth.columns)}")

mapping = pd.read_csv(CTB, sep="\t", names=["contig","MAG"])
log(f"Mapping columns: {list(mapping.columns)}")

merged = pd.merge(mapping, depth, on="contig")
abundance = merged.groupby("MAG").sum()
abundance.to_csv(MAG_ABUNDANCE, sep="\t")
log(f"✅ MAG abundance saved → {MAG_ABUNDANCE}")

# 2) MAG x KO TPM with length-based distribution
log("Parsing DIAMOND results...")
diamond = pd.read_csv(DIAMOND, sep="\t", header=None,
                      names=["query","subject","identity","align_length","mismatch","gap",
                             "q_start","q_end","s_start","s_end","evalue","bitscore"])
diamond["contig"] = diamond["query"].str.extract(r'~(k141_[^_]+)')
diamond["KO"] = diamond["subject"].str.extract(r'~(K\d+)')
diamond = diamond.dropna(subset=["contig","KO"])
log(f"DIAMOND entries after filtering: {len(diamond)}")

# prepare merge
merged = diamond[["contig","KO","align_length"]].merge(depth, on="contig").merge(mapping, on="contig")
samples = [col for col in merged.columns if col not in ["contig","KO","MAG","align_length"]]

# length-based distribution
log("Applying length-based distribution for contig TPM...")
merged["len_sum"] = merged.groupby("contig")["align_length"].transform("sum")
for sample in samples:
    merged[sample] = merged[sample] * (merged["align_length"] / merged["len_sum"])

# === (A) MAG × KO TPM ===
tpm_table = merged.groupby(["MAG","KO"])[samples].sum().reset_index()

# === (B) total KO TPM (all MAGs together) ===
ko_table = merged.groupby(["KO"])[samples].sum().reset_index()

# KO annotation
log("Merging KO annotation and pathway info...")
desc = pd.read_csv(f"{KO_INFO}/ko_list.txt", sep="\t", header=None, names=["KO","Description"])
desc["KO"] = desc["KO"].str.replace("ko:","",regex=False)
desc["Gene"] = desc["Description"].str.extract(r'^([a-zA-Z0-9_-]+)')
desc["Description"] = desc["Description"].str.replace(r'^([a-zA-Z0-9_-]+)\s+','',regex=True)

ko2path = pd.read_csv(f"{KO_INFO}/ko_to_pathway.txt", sep="\t", header=None, names=["KO","Pathway"])
ko2path["KO"] = ko2path["KO"].str.replace("ko:","",regex=False)
ko2path["Pathway"] = ko2path["Pathway"].str.replace("path:","",regex=False)

path_name = pd.read_csv(f"{KO_INFO}/pathway_list.txt", sep="\t", header=None, names=["Pathway","Pathway_Name"])
path_name["Pathway"] = path_name["Pathway"].str.replace("path:","",regex=False)

# (A) MAG × KO TPM annotation
annot = tpm_table.merge(desc[["KO","Gene","Description"]], on="KO", how="left")
annot = annot.merge(ko2path, on="KO", how="left")
annot = annot.merge(path_name, on="Pathway", how="left")

grouped = annot.groupby(["MAG","KO","Gene","Description"], as_index=False).agg({
    "Pathway": lambda x: "; ".join(sorted(set(x.dropna()))),
    "Pathway_Name": lambda x: "; ".join(sorted(set(x.dropna()))),
    **{col:"first" for col in samples}
})

cols = ["MAG","KO","Gene","Description","Pathway","Pathway_Name"] + samples
grouped[cols].to_csv(MAG_KO_TPM, sep="\t", index=False)
log(f"✅ MAG × KO TPM + annotation saved → {MAG_KO_TPM}")

# (B) KO TPM annotation
annot2 = ko_table.merge(desc[["KO","Gene","Description"]], on="KO", how="left")
annot2 = annot2.merge(ko2path, on="KO", how="left")
annot2 = annot2.merge(path_name, on="Pathway", how="left")

grouped2 = annot2.groupby(["KO","Gene","Description"], as_index=False).agg({
    "Pathway": lambda x: "; ".join(sorted(set(x.dropna()))),
    "Pathway_Name": lambda x: "; ".join(sorted(set(x.dropna()))),
    **{col:"first" for col in samples}
})

cols2 = ["KO","Gene","Description","Pathway","Pathway_Name"] + samples
grouped2[cols2].to_csv(KO_TPM, sep="\t", index=False)
log(f"✅ KO TPM + annotation saved → {KO_TPM}")

log("=== Python step completed ===")
EOF

echo "🎉 Pipeline finished!" | tee -a "$LOGFILE"

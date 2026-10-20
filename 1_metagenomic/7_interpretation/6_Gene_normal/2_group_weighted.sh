#!/bin/bash
# Requires: Python 3 (pandas)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

################################################################################
# USER SETTINGS
# Only edit this section
################################################################################

# Directory containing previous output files
INDIR="/path/to/base/1_metagenomic/7_interpretation/6_Gene_normal"

# Main MAG-level activity file from previous pipeline
MAG_ACTIVITY_FILE="$INDIR/MAG_KO_activity_RNA_over_DNAcopy.tsv"

# Output directory
OUTDIR="$INDIR/grouped_MAG_activity"
mkdir -p "$OUTDIR"

# Minimum DNA copy threshold for a MAG-KO to be included in each condition
# Recommended: 0.1 for permissive, 0.5 for conservative
MIN_DNA_COPY=0.1

# Conditions used in the MAG_KO_activity file
CONDITIONS="2T 3T 4T 6T 7T"

################################################################################
# DO NOT EDIT BELOW THIS LINE
################################################################################

echo "============================================================"
echo "MAG functional activity summary"
echo "============================================================"
echo "Input file      = $MAG_ACTIVITY_FILE"
echo "Output dir      = $OUTDIR"
echo "MIN_DNA_COPY    = $MIN_DNA_COPY"
echo "CONDITIONS      = $CONDITIONS"
echo "============================================================"

if [[ ! -s "$MAG_ACTIVITY_FILE" ]]; then
    echo "ERROR: input file not found or empty: $MAG_ACTIVITY_FILE" >&2
    exit 1
fi

################################################################################
# Functional KO group file
################################################################################

KO_GROUP_FILE="$OUTDIR/KO_function_groups.tsv"

cat > "$KO_GROUP_FILE" << 'EOF'
KO	Gene	Function_Group
K02575	NRT2/narK/nrtP/nasA	Nitrate_transport
K15576	nrtA/nasF/cynA	Nitrate_transport
K15577	nrtB/nasE/cynB	Nitrate_transport
K15578	nrtC/nasD	Nitrate_transport

K00370	narG/narZ/nxrA	Respiratory_nitrate_reduction
K00371	narH/narY/nxrB	Respiratory_nitrate_reduction
K00374	narI/narV	Respiratory_nitrate_reduction
K02567	napA	Periplasmic_nitrate_reduction
K02568	napB	Periplasmic_nitrate_reduction

K00360	nasB	Assimilatory_nitrate_reduction
K00372	nasC/nasA	Assimilatory_nitrate_reduction
K00367	narB	Assimilatory_nitrate_reduction

K02586	nifD	Nitrogen_fixation
K02588	nifH	Nitrogen_fixation
K02591	nifK	Nitrogen_fixation

K00362	nirB	DNRA_or_assimilatory_nitrite_reduction
K00363	nirD	DNRA_or_assimilatory_nitrite_reduction
K03385	nrfA	DNRA
K15876	nrfH	DNRA

K15864	nirS	Denitrification_nitrite_reduction
K00368	nirK	Denitrification_nitrite_reduction

K00366	nirA	Assimilatory_nitrite_reduction

K04561	norB	NO_reduction
K02305	norC	NO_reduction
K00376	nosZ	N2O_reduction

K05601	hcp	NO_detox_or_hybrid_cluster_protein
K12264	norV	NO_detox
K12265	norW	NO_detox
K05916	hmp	NO_detox

K00265	gltB	N_assimilation
K00266	gltD	N_assimilation
K01915	glnA/GLUL	N_assimilation

K00436	hoxH	Bidirectional_hydrogenase
K18005	hoxF	Bidirectional_hydrogenase
K18006	hoxU	Bidirectional_hydrogenase
K18007	hoxY	Bidirectional_hydrogenase

K18008	hydA	FeFe_hydrogenase
K00437	hydB	FeFe_hydrogenase

K06281	hyaB	NiFe_uptake_hydrogenase
K06282	hyaA	NiFe_uptake_hydrogenase

K14109	echA	Energy_converting_hydrogenase
K14110	echB	Energy_converting_hydrogenase
K14111	echC	Energy_converting_hydrogenase
K14112	echD	Energy_converting_hydrogenase
K14113	echE	Energy_converting_hydrogenase
K14114	echF	Energy_converting_hydrogenase

K00440	mvhA	Methyl_viologen_reducing_hydrogenase
K00441	mvhG	Methyl_viologen_reducing_hydrogenase
K00442	mvhD	Methyl_viologen_reducing_hydrogenase

K04564	sodA	Oxidative_stress_SOD
K04565	sodB	Oxidative_stress_SOD
K04566	sodC	Oxidative_stress_SOD

K03781	katG	Oxidative_stress_catalase
K03782	katE	Oxidative_stress_catalase

K03386	ahpC	Oxidative_stress_peroxiredoxin
K03387	ahpF	Oxidative_stress_peroxiredoxin

K00383	gor	Redox_buffering_glutathione
K03553	recA	DNA_damage_repair

K04761	oxyR	Oxidative_stress_regulator
K03711	fur	Iron_redox_regulation

K00425	cydA	Terminal_oxidase
K00426	cydB	Terminal_oxidase

K00404	ccoN	Terminal_oxidase
K00405	ccoO	Terminal_oxidase
K00406	ccoP	Terminal_oxidase
EOF

################################################################################
# Python summary
################################################################################

python3 - << PY
import os
import pandas as pd

INDIR = "$INDIR"
OUTDIR = "$OUTDIR"
MAG_ACTIVITY_FILE = "$MAG_ACTIVITY_FILE"
KO_GROUP_FILE = "$KO_GROUP_FILE"
MIN_DNA_COPY = float("$MIN_DNA_COPY")
CONDITIONS = "$CONDITIONS".split()

print("Reading MAG activity file...")
df = pd.read_csv(MAG_ACTIVITY_FILE, sep="\t")

print("Reading KO group file...")
groups = pd.read_csv(KO_GROUP_FILE, sep="\t")
groups = groups.dropna(subset=["KO", "Function_Group"]).drop_duplicates()

required = {"MAG", "KO"}
if not required.issubset(df.columns):
    raise SystemExit(f"ERROR: MAG activity file must contain columns {required}. Found: {list(df.columns)}")

# Merge group labels
dfg = df.merge(groups, on="KO", how="inner")

if dfg.empty:
    raise SystemExit("ERROR: No matching KO found between MAG activity file and KO group file.")

# Check available condition-specific columns
available = []
missing = []

for cond in CONDITIONS:
    rna_col = f"zr15207_{cond}"
    dna_col = None
    ratio_col = f"{cond}_RNA_over_DNAcopy"
    
    # DNA copy columns are named by DNA sample, not condition.
    # Map condition to DNA sample from previous pipeline.
    cond_to_dna = {
        "2T": "zr15207_9.sorted_DNA_copy_per_MAG",
        "3T": "zr15207_10.sorted_DNA_copy_per_MAG",
        "4T": "zr15207_11.sorted_DNA_copy_per_MAG",
        "6T": "zr15207_13.sorted_DNA_copy_per_MAG",
        "7T": "zr15207_14.sorted_DNA_copy_per_MAG",
    }
    dna_col = cond_to_dna.get(cond)

    if rna_col in dfg.columns and dna_col in dfg.columns and ratio_col in dfg.columns:
        available.append((cond, rna_col, dna_col, ratio_col))
    else:
        missing.append((cond, rna_col, dna_col, ratio_col))

if missing:
    print("WARNING: Some condition columns were not found and will be skipped:")
    for x in missing:
        print("  ", x)

if not available:
    raise SystemExit("ERROR: No usable condition columns found.")

print("Available conditions:")
for x in available:
    print("  ", x)

# Convert numeric columns
for _, rna_col, dna_col, ratio_col in available:
    dfg[rna_col] = pd.to_numeric(dfg[rna_col], errors="coerce").fillna(0)
    dfg[dna_col] = pd.to_numeric(dfg[dna_col], errors="coerce").fillna(0)
    dfg[ratio_col] = pd.to_numeric(dfg[ratio_col], errors="coerce")

# Long-format summary
records = []

for cond, rna_col, dna_col, ratio_col in available:
    sub = dfg.copy()
    
    # condition-specific DNA filter
    sub = sub[sub[dna_col] >= MIN_DNA_COPY].copy()
    
    if sub.empty:
        continue
    
    grouped = sub.groupby(["MAG", "Function_Group"], dropna=False)
    
    for (mag, func), g in grouped:
        sum_rna = g[rna_col].sum()
        sum_dna = g[dna_col].sum()
        weighted_activity = sum_rna / sum_dna if sum_dna > 0 else None
        
        ratios = g[ratio_col].dropna()
        mean_ratio = ratios.mean() if len(ratios) else None
        median_ratio = ratios.median() if len(ratios) else None
        sum_ratio = ratios.sum() if len(ratios) else None
        
        records.append({
            "MAG": mag,
            "Function_Group": func,
            "Condition": cond,
            "n_KO_used": g["KO"].nunique(),
            "KO_list_used": ";".join(sorted(g["KO"].unique())),
            "Gene_list_used": ";".join(sorted(g["Gene"].dropna().astype(str).unique())),
            "sum_RNA_TPM": sum_rna,
            "sum_DNA_copy": sum_dna,
            "weighted_activity_sumRNA_over_sumDNA": weighted_activity,
            "mean_RNA_over_DNAcopy": mean_ratio,
            "median_RNA_over_DNAcopy": median_ratio,
            "sum_RNA_over_DNAcopy": sum_ratio
        })

summary_long = pd.DataFrame(records)

if summary_long.empty:
    raise SystemExit("ERROR: Summary table is empty after filtering. Try lowering MIN_DNA_COPY.")

out_long = os.path.join(OUTDIR, f"MAG_function_group_activity.long_DNAcopyge{MIN_DNA_COPY}.tsv")
summary_long.to_csv(out_long, sep="\t", index=False)
print("Written:", out_long)

# Wide tables for plotting
metrics = [
    "sum_RNA_TPM",
    "sum_DNA_copy",
    "weighted_activity_sumRNA_over_sumDNA",
    "mean_RNA_over_DNAcopy",
    "median_RNA_over_DNAcopy",
    "n_KO_used"
]

for metric in metrics:
    wide = summary_long.pivot_table(
        index=["MAG", "Function_Group"],
        columns="Condition",
        values=metric,
        aggfunc="first"
    ).reset_index()
    wide.columns.name = None
    
    out_wide = os.path.join(
        OUTDIR,
        f"MAG_function_group_{metric}.wide_DNAcopyge{MIN_DNA_COPY}.tsv"
    )
    wide.to_csv(out_wide, sep="\t", index=False)
    print("Written:", out_wide)

# Top MAGs per function and condition by weighted activity
top_records = []
for (func, cond), g in summary_long.groupby(["Function_Group", "Condition"]):
    gg = g.sort_values("weighted_activity_sumRNA_over_sumDNA", ascending=False)
    gg = gg.head(20).copy()
    gg["rank_within_function_condition"] = range(1, len(gg) + 1)
    top_records.append(gg)

top = pd.concat(top_records, ignore_index=True)
out_top = os.path.join(OUTDIR, f"Top20_MAGs_by_function_activity_DNAcopyge{MIN_DNA_COPY}.tsv")
top.to_csv(out_top, sep="\t", index=False)
print("Written:", out_top)

# Function-level total across all MAGs
func_records = []
for cond, rna_col, dna_col, ratio_col in available:
    sub = dfg[dfg[dna_col] >= MIN_DNA_COPY].copy()
    for func, g in sub.groupby("Function_Group"):
        sum_rna = g[rna_col].sum()
        sum_dna = g[dna_col].sum()
        func_records.append({
            "Function_Group": func,
            "Condition": cond,
            "n_MAG_used": g["MAG"].nunique(),
            "n_KO_used": g["KO"].nunique(),
            "sum_RNA_TPM": sum_rna,
            "sum_DNA_copy": sum_dna,
            "weighted_activity_sumRNA_over_sumDNA": sum_rna / sum_dna if sum_dna > 0 else None,
            "KO_list_used": ";".join(sorted(g["KO"].unique()))
        })

func_summary = pd.DataFrame(func_records)
out_func = os.path.join(OUTDIR, f"Community_from_MAG_function_group_activity_DNAcopyge{MIN_DNA_COPY}.tsv")
func_summary.to_csv(out_func, sep="\t", index=False)
print("Written:", out_func)

print("\\nDONE.")
print("Main output:")
print(out_long)
PY

echo ""
echo "Finished."
echo "Results are in:"
echo "$OUTDIR"

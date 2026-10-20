#!/bin/bash
# Requires: bash, awk
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

set -euo pipefail

# Master output directory for all bin-to-contig mapping tables
TSV_OUT_DIR="/path/to/base/1_metagenomic/6_bincheck/1_tsv"
mkdir -p "$TSV_OUT_DIR"
echo "✅ Saving all .tsv files to: $TSV_OUT_DIR"

################################
# 1) VAMB
################################
VAMB_DIR="/path/to/base/1_metagenomic/5_binning/1_vamb/vamb_out/bins"
VAMB_TSV="$TSV_OUT_DIR/vamb.tsv"

echo "🔄 Creating vamb.tsv (prefix 'vamb_')..."
: > "$VAMB_TSV"
shopt -s nullglob
for f in "$VAMB_DIR"/*.fna; do
  bin="$(basename "$f" .fna)"
  # extract contig IDs, assign bin with 'vamb_' prefix
  grep "^>" "$f" | sed 's/^>//' | awk -v b="vamb_$bin" '{print $1 "\t" b}' >> "$VAMB_TSV"
done

################################
# 2) MaxBin2
################################
MAXBIN_DIR="/path/to/base/1_metagenomic/5_binning/3_maxbin2/maxbin_out"
MAXBIN_TSV="$TSV_OUT_DIR/maxbin.tsv"

echo "🔄 Creating maxbin.tsv (replace '.' with '_')..."
: > "$MAXBIN_TSV"

# support both *.fasta and *.fa
FOUND=0
for f in "$MAXBIN_DIR"/*.fasta "$MAXBIN_DIR"/*.fa; do
  [[ -e "$f" ]] || continue
  FOUND=1
  rawbin="$(basename "$f")"
  rawbin="${rawbin%.fasta}"
  rawbin="${rawbin%.fa}"
  bin="${rawbin//./_}"
  grep "^>" "$f" | sed 's/^>//' | awk -v b="$bin" '{print $1 "\t" b}' >> "$MAXBIN_TSV"
done
if [[ $FOUND -eq 0 ]]; then
  echo "⚠️  No MaxBin2 *.fasta/*.fa files found in $MAXBIN_DIR"
fi

################################
# 3) SemiBin2
################################
SEMIBIN_IN="/path/to/base/1_metagenomic/5_binning/4_semibin/3_bin/contig_bins.tsv"
SEMIBIN_TSV="$TSV_OUT_DIR/semibin.tsv"

echo "🔄 Creating semibin.tsv (prefix 'semibin_')..."
if [[ -f "$SEMIBIN_IN" ]]; then
  awk 'BEGIN{OFS="\t"} NR>1 {$2="semibin_"$2; print $1,$2}' "$SEMIBIN_IN" > "$SEMIBIN_TSV"
else
  echo "⚠️  SemiBin contig_bins.tsv not found: $SEMIBIN_IN"
fi

################################
# 4) CONCOCT
################################
CLUSTERING="/path/to/base/1_metagenomic/5_binning/2_concoct/concoct_out/clustering_merged.csv"
CONCOCT_TSV="$TSV_OUT_DIR/concoct_fixed.tsv"

awk -F',' 'NR>1{print $1"\tconcoct_bin_"$2}' \
  "$CLUSTERING" \
  > "$CONCOCT_TSV"

echo "🎉 Done. Generated TSV files:"
ls -lh "$TSV_OUT_DIR"/*.tsv || true

#!/bin/bash
# Requires: CONCOCT 1.1.0
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
BAM_DIR="/path/to/base/1_metagenomic/3_mapping"
OUTDIR="/path/to/base/1_metagenomic/5_binning/2_concoct/concoct_out"
MINFASTA=200000
THREADS=8
# 32G memory will be enough

#Step 4: clustering results -> merge accroding to original contig
merge_cutup_clustering.py \
  "$OUTDIR/clustering_gt1000.csv" \
  > "$OUTDIR/clustering_merged.csv"
echo "CONCOCT clustring.py finished"

#Step 5: extract the fasta for each bin from orignial contig FASTA
mkdir -p "$OUTDIR/bins"

extract_fasta_bins.py \
  "$FASTA" \
  "$OUTDIR/clustering_merged.csv" \
  --output_path "$OUTDIR/bins"

echo "✅ CONCOCT fasta_bin_binning finished!"

#Step 6: flitering bin size >200,000bp
mkdir -p "$OUTDIR/bins_200k"

for f in "$OUTDIR"/bins/*.fa; do
  size=$(awk 'BEGIN{sum=0} /^>/ {next} {sum+=length($0)} END{print sum}' "$f")
  if [ "$size" -ge "$MINFASTA" ]; then
    cp "$f" "$OUTDIR/bins_200k/"
  fi
done

echo "flitering bin size"

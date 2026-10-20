#!/bin/bash
# Requires: bash, awk
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Paths and settings
FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
DEPTH_TSV="/path/to/base/1_metagenomic/4_coverage/coverm_contig_mean.tsv"
OUTDIR="/path/to/base/1_metagenomic/5_binning/3_maxbin2/maxbin_out"
THREADS=8
MINFASTA=200000
# 4-8 CPU and 8G memory will be enough

#flitering bin size >200,000bp
mkdir -p "$OUTDIR/bins_200k"

for f in "$OUTDIR"/*.fasta; do
  size=$(awk 'BEGIN{sum=0} /^>/ {next} {sum+=length($0)} END{print sum}' "$f")
  if [ "$size" -ge "$MINFASTA" ]; then
    cp "$f" "$OUTDIR/bins_200k/"
  fi
done

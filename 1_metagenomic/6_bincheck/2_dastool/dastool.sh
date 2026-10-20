#!/bin/bash
#   script to run das_tool
#
# Requires: DAS Tool 1.1.7, DIAMOND 0.9.36
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Input files
TSVDIR="/path/to/base/1_metagenomic/6_bincheck/1_tsv"
FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
OUTDIR="/path/to/base/1_metagenomic/6_bincheck/2_dastool"
THREADS=16

# Run DAS Tool
DAS_Tool \
  -i "$TSVDIR/vamb.tsv","$TSVDIR/maxbin.tsv","$TSVDIR/concoct_fixed.tsv","$TSVDIR/semibin.tsv" \
  -l vamb,maxbin,concoct,semibin \
  -c "$FASTA" \
  -o "$OUTDIR/DAS_Tool" \
  --threads "$THREADS" \
  --write_bins \
  --write_bin_evals

echo "✅ DAS Tool integration complete!"
echo "📂 Results are in: $OUTDIR"

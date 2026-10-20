#!/bin/bash
#
# Requires: VAMB 5.0.3
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
DEPTH_ORG="/path/to/base/1_metagenomic/4_coverage/coverm_contig_mean.tsv"
DEPTH_VAMB="/path/to/base/1_metagenomic/5_binning/1_vamb/coverm_depth_vamb.tsv"
OUTDIR="/path/to/base/1_metagenomic/5_binning/1_vamb/vamb_out"

## 8 CPU, 8G memory will be enough

# ✅ change header 'Contig' → 'contigname'
sed '1s/^Contig/contigname/' "$DEPTH_ORG" > "$DEPTH_VAMB"

# ✅ run VAMB
vamb bin default \
  --fasta "$FASTA" \
  --abundance_tsv "$DEPTH_VAMB" \
  --outdir "$OUTDIR" \
  --minfasta 200000 \

echo "✅ VAMB binning finished!"

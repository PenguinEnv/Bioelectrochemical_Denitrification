#!/bin/bash
# Requires: SemiBin2 2.2.0 (GPU optional)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Base paths
FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
MAPDIR="/path/to/base/1_metagenomic/3_mapping"
OUTBASE="/path/to/base/1_metagenomic/5_binning/4_semibin"
THREADS=32
MINFASTA=200000

mkdir -p "$OUTBASE"

# 1. Feature generation
SemiBin2 generate_sequence_features_single \
  -i "$FASTA" \
  -b \
    "$MAPDIR/zr15207_10.sorted.bam" \
    "$MAPDIR/zr15207_11.sorted.bam" \
    "$MAPDIR/zr15207_13.sorted.bam" \
    "$MAPDIR/zr15207_14.sorted.bam" \
    "$MAPDIR/zr15207_9.sorted.bam" \
  -o "$OUTBASE/1_generate" \
  --threads "$THREADS"

# 2. Model training
SemiBin2 train_self \
  --data "$OUTBASE/1_generate/data.csv" \
  --data-split "$OUTBASE/1_generate/data_split.csv" \
  -o "$OUTBASE/2_train" \
  --threads "$THREADS"

# 3. Binning
SemiBin2 bin_short \
  -i "$FASTA" \
  --model "$OUTBASE/2_train/model.pt" \
  --data "$OUTBASE/1_generate/data.csv" \
  -o "$OUTBASE/3_bin" \
  --threads "$THREADS"

# Filtering bins only >200K  
mkdir -p "$OUTBASE/bins_200k"

for f in "$OUTBASE"/3_bin/output_bins/*.fa.gz; do
  size=$(awk 'BEGIN{sum=0} /^>/ {next} {sum+=length($0)} END{print sum}' "$f")
  if [ "$size" -ge "$MINFASTA" ]; then
    cp "$f" "$OUTBASE/bins_200k/"
  fi
done

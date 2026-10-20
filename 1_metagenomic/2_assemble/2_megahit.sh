#!/bin/bash
# Requires: MEGAHIT 1.2.9
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

THREADS=16
TRIMMED_DIR="/path/to/base/1_metagenomic/2_assemble"
MEGAHIT_OUT="/path/to/base/1_metagenomic/2_assemble/megahit_out"

# Run megahit

megahit \
	-1 "$TRIMMED_DIR/all_R1.fastq.gz" \
	-2 "$TRIMMED_DIR/all_R2.fastq.gz" \
	-o "$MEGAHIT_OUT" \
	-t 16 \
	--min-contig-len 1000

echo "✅ MEGAHIT finished for all samples!"

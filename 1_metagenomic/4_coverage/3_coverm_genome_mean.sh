#!/bin/bash
# Requires: CoverM 0.7.0
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

set -euo pipefail

THREADS=4

# Input directories/files
BAM_DIR="/path/to/base/1_metagenomic/3_mapping"
MAG_DIR="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"

# Output directories/files
WORKDIR="/path/to/base/1_metagenomic/4_coverage"
GENOME_TPM="$WORKDIR/coverm_genome_mean.tsv"

mkdir -p "$WORKDIR"

# ==========================
# 1️⃣ Run CoverM genome → TPM
# ==========================

echo "🔄 Running CoverM genome tpm..."
coverm genome\
  --bam-files "$BAM_DIR"/*.sorted.bam \
  --genome-fasta-directory "$MAG_DIR" \
  --genome-fasta-extension fa \
  --methods mean \
  --min-read-percent-identity 0.95 \
  --min-read-aligned-percent 0.75 \
  --threads $THREADS \
  --output-file "$GENOME_TPM"
echo "✅ CoverM genome TPM done: $GENOME_TPM"

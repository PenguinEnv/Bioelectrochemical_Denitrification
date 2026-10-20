#!/bin/bash
# Requires: CoverM 0.7.0
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

set -euo pipefail

THREADS=8

# Input directories/files
BAM_DIR="/path/to/base/1_metagenomic/3_mapping"

# Output directories/files
WORKDIR="/path/to/base/1_metagenomic/4_coverage"
CONTIG_TPM="$WORKDIR/coverm_contig_tpm.tsv"
CONTIG_MEAN="$WORKDIR/coverm_contig_mean.tsv"

mkdir -p "$WORKDIR"

# ==========================
# 1️⃣ Run CoverM contig 
# ==========================
echo "🔄 Running CoverM contig tpm..."
coverm contig \
  --bam-files "$BAM_DIR"/*.sorted.bam \
  --methods tpm \
  --min-read-percent-identity 0.95 \
  --min-read-aligned-percent 0.75 \
  --threads $THREADS \
  --output-file "$CONTIG_TPM"
echo "✅ CoverM contig TPM done: $CONTIG_TPM"

echo "🔄 Running CoverM contig mean..."
coverm contig \
  --bam-files "$BAM_DIR"/*.sorted.bam \
  --methods mean \
  --min-read-percent-identity 0.95 \
  --min-read-aligned-percent 0.75 \
  --threads $THREADS \
  --output-file "$CONTIG_MEAN"
echo "✅ CoverM contig MEAN done: $CONTIG_MEAN"

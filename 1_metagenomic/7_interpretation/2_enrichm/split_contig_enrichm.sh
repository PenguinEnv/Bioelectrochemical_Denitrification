#!/bin/bash
# --------------------------------------------------------------------
#  EnrichM v0.6.5 : sequential annotation using --genome_directory
#  Author: Hyejeong | Date: 2025-10-xx
# --------------------------------------------------------------------
# Requires: EnrichM 0.6.5, DIAMOND, SeqKit
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# === Paths ===
FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
WORKDIR="/path/to/base/1_metagenomic/7_interpretation/2_enrichm"
CONTIG="$WORKDIR/final.contigs.fna"
SPLIT_DIR="$WORKDIR/split_contigs"
OUTDIR="$WORKDIR/annot_output"

THREADS=8
SPLIT_SIZE=10000   # number of contigs per chunk

# === Preparation ===
mkdir -p "$WORKDIR" "$SPLIT_DIR" 
cd "$WORKDIR" || exit 1

exec > >(tee enrichm_run.log)
exec 2>&1

cp "$FASTA" "$CONTIG"

# === Split contigs ===
echo "✂️  Splitting FASTA into chunks of ~$SPLIT_SIZE contigs..."
seqkit split -s "$SPLIT_SIZE" "$CONTIG" -O "$SPLIT_DIR"

# === Run EnrichM (per directory) ===
echo "🚀 Running EnrichM annotate on directory..."
enrichm annotate \
    --ko \
    --genome_directory "$SPLIT_DIR" \
    --output "$OUTDIR" \
    --threads "$THREADS" 

echo "✅ EnrichM annotation completed."

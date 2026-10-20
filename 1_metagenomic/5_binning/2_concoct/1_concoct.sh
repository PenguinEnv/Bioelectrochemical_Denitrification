#!/bin/bash
# Requires: CONCOCT 1.1.0, SAMtools 1.18
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
BAM_DIR="/path/to/base/1_metagenomic/3_mapping"
OUTDIR="/path/to/base/1_metagenomic/5_binning/2_concoct/concoct_out"
MINFASTA=200000
THREADS=16
# 32G memory will be enough

mkdir -p "$OUTDIR"

# Step 1: contig cut-up and clustering
cut_up_fasta.py "$FASTA" -c 10000 -o 0 --merge_last -b "$OUTDIR/contigs_10K.bed" > "$OUTDIR/contigs_10K.fa"

# Step 2: coverage table from BAM files
concoct_coverage_table.py "$OUTDIR/contigs_10K.bed" "$BAM_DIR"/*.sorted.bam > "$OUTDIR/coverage_table.tsv"

# Step 3: binning
concoct --composition_file "$OUTDIR/contigs_10K.fa" \
        --coverage_file "$OUTDIR/coverage_table.tsv" \
        -b "$OUTDIR/" \
        -c 200 \
        --threads "$THREADS"

#!/usr/bin/env bash
# Requires: bash, gzip
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.
set -euo pipefail

SRC_DIR="/path/to/base/1_metagenomic/1_qc/1_fastp_out"
OUT_DIR="/path/to/base/1_metagenomic/2_assemble"
LOGFILE="$OUT_DIR/combine_reads.log"

mkdir -p "$OUT_DIR"

# remove previous output (avoid appending twice)
rm -f "$OUT_DIR/all_R1.fastq.gz" "$OUT_DIR/all_R2.fastq.gz"

echo "[$(date)] Start: Combining trimmed reads" > "$LOGFILE"

cd "$SRC_DIR"

# combine R1
echo "[$(date)] Combining R1 reads" >> "$LOGFILE"
for R1 in zr15207*_clean_R1.fastq.gz; do
    echo "[$(date)] Added: $R1" >> "$LOGFILE"
    cat "$R1" >> "$OUT_DIR/all_R1.fastq.gz"
done
echo "[$(date)] Finished combining R1 reads" >> "$LOGFILE"

# combine R2
echo "[$(date)] Combining R2 reads" >> "$LOGFILE"
for R2 in zr15207*_clean_R2.fastq.gz; do
    echo "[$(date)] Added: $R2" >> "$LOGFILE"
    cat "$R2" >> "$OUT_DIR/all_R2.fastq.gz"
done
echo "[$(date)] Finished combining R2 reads" >> "$LOGFILE"

# gzip integrity check
gzip -t "$OUT_DIR/all_R1.fastq.gz" && gzip -t "$OUT_DIR/all_R2.fastq.gz" \
  && echo "[$(date)] Gzip integrity OK" >> "$LOGFILE"

# record file sizes
echo "[$(date)] Combined file sizes:" >> "$LOGFILE"
ls -lh "$OUT_DIR"/all_R*.fastq.gz >> "$LOGFILE"

echo "[$(date)] Done: All trimmed reads combined" >> "$LOGFILE"

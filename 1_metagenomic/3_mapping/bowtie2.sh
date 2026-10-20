#!/usr/bin/env bash
# Requires: Bowtie2 2.5.1, SAMtools 1.18
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

set -euo pipefail

# ===== Paths & parameters =====
CONTIGS="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
READ_DIR="/path/to/base/1_metagenomic/1_qc/1_fastp_out"
OUT_DIR="/path/to/base/1_metagenomic/3_mapping"
INDEX_PREFIX="$OUT_DIR/megahit_contigs_index"   # required: Bowtie2 index prefix
THREADS=16

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ===== 1) Build Bowtie2 index (skip if already exists) =====
echo "[$(date)] Building Bowtie2 index if needed..."
if ls "${INDEX_PREFIX}".*.bt2 >/dev/null 2>&1 || ls "${INDEX_PREFIX}".*.bt2l >/dev/null 2>&1; then
  echo "[$(date)] Index already exists at ${INDEX_PREFIX}*"
else
  bowtie2-build "$CONTIGS" "$INDEX_PREFIX"
  echo "[$(date)] Indexing complete."
fi

# ===== 2) Map each sample based on R1 file pattern =====
shopt -s nullglob
R1_LIST=( "$READ_DIR"/zr15207*_clean_R1.fastq.gz )

if [[ ${#R1_LIST[@]} -eq 0 ]]; then
  echo "[$(date)] ERROR: No R1 files found in $READ_DIR (pattern: zr15207*_clean_R1.fastq.gz)"
  exit 1
fi

echo "[$(date)] Found ${#R1_LIST[@]} R1 files. Starting mapping..."

for R1 in "${R1_LIST[@]}"; do
  base=$(basename "$R1")
  sample=${base%_clean_R1.fastq.gz}
  R2="$READ_DIR/${sample}_clean_R2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "[$(date)] SKIP: $sample (missing R2: $R2)"
    continue
  fi

  echo "[$(date)] Mapping $sample ..."
  BAM="${sample}.sorted.bam"
  BWTLOG="${sample}.bowtie2.log"
  FLAGSTAT="${sample}.flagstat.txt"

  bowtie2 \
    -x "$INDEX_PREFIX" \
    -1 "$R1" -2 "$R2" \
    -p "$THREADS" \
    --very-sensitive \
    2> "$BWTLOG" \
  | samtools view -@ "$THREADS" -b \
  | samtools sort -@ "$THREADS" -o "$BAM" -

  samtools index "$BAM"
  samtools flagstat -@ "$THREADS" "$BAM" > "$FLAGSTAT"

  echo "[$(date)] Finished $sample  ->  $BAM"
done

echo "[$(date)] ✅ All samples mapped and sorted BAM files generated in: $OUT_DIR"

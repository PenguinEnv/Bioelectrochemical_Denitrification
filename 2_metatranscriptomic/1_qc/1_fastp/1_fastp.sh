#!/bin/bash
# Requires: fastp 1.0.1, MultiQC
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Paths
RAW_DIR="/path/to/base/your_metatranscriptome_raw_reads"
OUT_DIR="/path/to/base/2_metatranscriptomic/1_qc/1_fastp/fastp_out"
THREADS=8

# create output directory
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# fastp loop
for r1 in ${RAW_DIR}/*_R1.fastq.gz
do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${RAW_DIR}/${sample}_R2.fastq.gz"

    echo "▶️ Processing $sample..."

    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${OUT_DIR}/${sample}_clean_R1.fastq.gz" \
        -O "${OUT_DIR}/${sample}_clean_R2.fastq.gz" \
	--cut_tail \
        --trim_poly_g \
	-q 20 -u 30 -n 5 -l 30 -w "$THREADS" \
        -h "${OUT_DIR}/${sample}_fastp.html" \
        -j "${OUT_DIR}/${sample}_fastp.json"
done

echo "📊 Running MultiQC..."
multiqc "$OUT_DIR" -o "$OUT_DIR"

echo "✅ All done. Cleaned data and reports are in: $OUT_DIR"

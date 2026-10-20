#!/bin/bash
# Requires: fastp 1.0.1, MultiQC
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Edit
LOGFILE="log_temp.txt"
# ----------

THREADS=8
INPUT_DIR="/path/to/base/your_metagenome_raw_reads"
FASTP_OUT="./1_fastp_out"
REPORT_DIR="./2_fastp_multiqc_reports"
ADAPTER_SEQ="CTGTCTCTTATACACATCT"  # Nextera adapter
LOG="run_log.txt"

# Create output directories
mkdir -p "$FASTP_OUT" "$REPORT_DIR"

# Sample auto-detection
samples=($(ls $INPUT_DIR/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1.fastq.gz//' | sort))

# Start processing
for sample in "${samples[@]}"; do
    echo "🔬 Processing $sample..." | tee -a "$LOG"

    R1="$INPUT_DIR/${sample}_R1.fastq.gz"
    R2="$INPUT_DIR/${sample}_R2.fastq.gz"
    CLEAN_R1="$FASTP_OUT/${sample}_clean_R1.fastq.gz"
    CLEAN_R2="$FASTP_OUT/${sample}_clean_R2.fastq.gz"

    # 1️⃣ fastp
    fastp \
        -i "$R1" -I "$R2" \
        -o "$CLEAN_R1" -O "$CLEAN_R2" \
        --adapter_sequence "$ADAPTER_SEQ" \
        --adapter_sequence_r2 "$ADAPTER_SEQ" \
        --trim_poly_g \
        --cut_tail \
        -q 20 -l 50 \
	-w "$THREADS" \
        -h "$REPORT_DIR/${sample}_fastp.html" \
	-j "$REPORT_DIR/${sample}_fastp.json" \
	>> "$LOG" 2>&1

    echo "✅ fastp done for $sample" | tee -a "$LOG"
    echo "----------------------------------------------" | tee -a "$LOG"
done

# 2️⃣ multiQC on fastp reports
multiqc "$REPORT_DIR" -o "$REPORT_DIR" >> "$LOG" 2>&1
echo "📊 MultiQC report saved at: $REPORT_DIR/multiqc_report.html"

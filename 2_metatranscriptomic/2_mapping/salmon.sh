#!/bin/bash
# Requires: Salmon 1.10.3
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Set directories
CONTIGS="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
SORTMERNA_DIR="/path/to/base/2_metatranscriptomic/1_qc/2_sortmerna"
WORK_DIR="/path/to/base/2_metatranscriptomic/2_mapping"
INDEX_DIR="${WORK_DIR}/salmon_index"

mkdir -p "$INDEX_DIR"

# Step 1: Build Salmon index
echo "📦 Building Salmon index..."
salmon index -t "$CONTIGS" -i "$INDEX_DIR" -k 31

# Step 2: Run Salmon quantification
cd "$WORK_DIR"
echo "🚀 Starting Salmon quantification..."

for sample in zr15207_2T zr15207_3T zr15207_4T zr15207_6T zr15207_7T
do
    echo "🧬 Running Salmon for $sample..."
    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -r "${SORTMERNA_DIR}/${sample}_mRNA_SILVA.fq.gz" \
        -p 8 \
        --validateMappings \
        -o "${WORK_DIR}/${sample}_salmon_quant"
done

echo "✅ All done!"

#!/bin/bash
#
# Requires: SortMeRNA 4.3.4, SILVA 138.1 LSU/SSU NR99 references
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

samples=("zr15207_2T" "zr15207_3T" "zr15207_4T" "zr15207_6T" "zr15207_7T")

INPUT_DIR="/path/to/base/2_metatranscriptomic/1_qc/1_fastp/fastp_out"
OUTPUT_DIR="/path/to/base/2_metatranscriptomic/1_qc/2_sortmerna"
DB_PATH="/path/to/your_db/SILVA_138.1"

# Loop over samples
for sample in "${samples[@]}"
do
    echo "🔁 Running sortMeRNA for $sample..."

    sortmerna \
        --ref "$DB_PATH/SILVA_138.1_LSURef_NR99_tax_silva.fasta" \
        --ref "$DB_PATH/SILVA_138.1_SSURef_NR99_tax_silva.fasta" \
        --reads "$INPUT_DIR/${sample}_clean_R1.fastq.gz" \
        --reads "$INPUT_DIR/${sample}_clean_R2.fastq.gz" \
        --paired_in \
        --fastx \
        --other "$OUTPUT_DIR/${sample}_mRNA_SILVA" \
        --workdir "$OUTPUT_DIR/${sample}_sortmerna_tmp" \
        --threads 8

    echo "✅ Finished $sample"
done

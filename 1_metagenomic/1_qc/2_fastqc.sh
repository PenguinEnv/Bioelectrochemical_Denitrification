#!/bin/bash
#
# Requires: FastQC, MultiQC
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

FASTPOUT="/path/to/base/1_metagenomic/1_qc/1_fastp_out"
FASTQC_RESULT="/path/to/base/1_metagenomic/1_qc/fastqc_results"

mkdir -p "$FASTQC_RESULT"

fastqc "$FASTPOUT"/*.fastq.gz \
       -t 10 \
       -o "$FASTQC_RESULT"

multiqc "$FASTQC_RESULT" -o "$FASTQC_RESULT"

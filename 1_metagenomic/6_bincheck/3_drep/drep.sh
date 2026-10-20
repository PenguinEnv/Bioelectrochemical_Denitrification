#!/bin/bash
# Requires: dRep 3.2.2, CheckM
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Set paths
INPUT="/path/to/base/1_metagenomic/6_bincheck/2_dastool/DAS_Tool_DASTool_bins"
OUTPUT="/path/to/base/1_metagenomic/6_bincheck/3_drep"
THREADS=16
#16GB will be enough

mkdir -p "$OUTPUT"

# Run dRep with CheckM built-in scoring
dRep dereplicate "$OUTPUT" \
  -g "$INPUT"/*.fa \
  -p "$THREADS" \
  --checkM_method taxonomy_wf \
  --completeness 50 \
  --contamination 10 \
  --S_algorithm ANImf \
  --S_ani 0.95 \
  --cov_thresh 0.25

echo "✅ dRep dereplication complete!"
echo "📂 Results saved in: $OUTPUT"

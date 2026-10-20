#!/bin/bash
# Requires: MaxBin2 2.2.7
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

# Paths and settings
FASTA="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"
DEPTH_TSV="/path/to/base/1_metagenomic/4_coverage/coverm_contig_mean.tsv"
OUTDIR="/path/to/base/1_metagenomic/5_binning/3_maxbin2/maxbin_out"
THREADS=8
MINFASTA=200000
# 4-8 CPU and 8G memory will be enough

# Create output directory
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "[`date`] 🔍 Parsing CoverM depth file: $DEPTH_TSV"

# Read header and generate sample-wise .abund.txt files
HEADER=$(head -n 1 "$DEPTH_TSV")
IFS=$'\t' read -ra COLS <<< "$HEADER"

for i in "${!COLS[@]}"; do
    COL="${COLS[$i]}"

    if [[ "$COL" == "Contig" ]]; then
        continue
    fi

    # Clean sample name
    SAMPLE_NAME=$(echo "$COL" | sed 's/\.sorted//;s/ Mean//g')
    OUTFILE="${SAMPLE_NAME}.abund.txt"
    COL_NUM=$((i + 1))  # awk uses 1-based index

    echo "✔️  Generating $OUTFILE from column $COL_NUM ($COL)"
    awk -v col=$COL_NUM 'NR > 1 {print $1 "\t" $col}' "$DEPTH_TSV" > "$OUTFILE"
done

# Build abundance file list for MaxBin2
ABUND_ARGS=""
for FILE in *.abund.txt; do
    ABUND_ARGS+=" -abund $FILE"
done

echo "[`date`] 🚀 Running MaxBin2..."
run_MaxBin.pl \
    -contig "$FASTA" \
    $ABUND_ARGS \
    -out "$OUTDIR/maxbin2" \
    -thread "$THREADS"

echo "[`date`] ✅ MaxBin2 binning complete!"

#flitering bin size >200,000bp
mkdir -p "$OUTDIR/bins_200k"

for f in "$OUTDIR"/*.fasta; do
  size=$(awk 'BEGIN{sum=0} /^>/ {next} {sum+=length($0)} END{print sum}' "$f")
  if [ "$size" -ge "$MINFASTA" ]; then
    cp "$f" "$OUTDIR/bins_200k/"
  fi
done

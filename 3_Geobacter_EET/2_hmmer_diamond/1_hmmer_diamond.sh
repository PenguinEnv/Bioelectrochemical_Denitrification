#!/bin/bash
#
# Requires: Prodigal 2.6.3, HMMER 3.4, DIAMOND; Pfam-A.hmm
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

BIN_DIR="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"
PFAM="/path/to/your_db/Pfam/Pfam-A.hmm"
GEODB="/path/to/base/3_Geobacter_EET/1_blastdb/geobacter_ref/major_ref"
# ===========================================

bins=(semibin_173 maxbin2_085 concoct_bin_89_sub)

mkdir -p EET_check/{prodigal,pfam,diamond,logs}

# find the FASTA file
find_fasta () {
  local base="$1"
  for ext in fa fasta fna fas; do
    if [[ -f "${BIN_DIR}/${base}.${ext}" ]]; then
      echo "${BIN_DIR}/${base}.${ext}"
      return 0
    fi
  done
  # file without an extension
  if [[ -f "${BIN_DIR}/${base}" ]]; then
    echo "${BIN_DIR}/${base}"
    return 0
  fi
  return 1
}

# check that the databases exist
if [[ ! -f "$PFAM" ]]; then
  echo "[ERROR] Pfam file not found: $PFAM"
  exit 1
fi
if [[ ! -f "${GEODB}.dmnd" ]]; then
  echo "[ERROR] DIAMOND db not found: ${GEODB}.dmnd"
  exit 1
fi

for b in "${bins[@]}"; do
  echo "=============================="
  echo "[INFO] Processing: $b"

  fa=$(find_fasta "$b") || { echo "[ERROR] FASTA not found for $b in $BIN_DIR"; exit 1; }
  echo "[INFO] FASTA: $fa"

  # 1) Prodigal: contig -> proteins
  faa="EET_check/prodigal/${b}.faa"
  fna="EET_check/prodigal/${b}.genes.fna"
  gff="EET_check/prodigal/${b}.gff"

  echo "[INFO] Prodigal..."
  prodigal -i "$fa" -a "$faa" -d "$fna" -o "$gff" -p single \
    > "EET_check/logs/${b}.prodigal.log" 2>&1

  # 2) Pfam (hmmscan)
  domtbl="EET_check/pfam/${b}.pfam.domtblout"
  echo "[INFO] hmmscan (Pfam)..."
  hmmscan --cpu 16 --domtblout "$domtbl" "$PFAM" "$faa" \
    > "EET_check/logs/${b}.pfam.log" 2>&1

  # 3) DIAMOND blastp
  hits="EET_check/diamond/${b}.geobacter_ref.tsv"
  echo "[INFO] DIAMOND blastp..."
  diamond blastp -d "$GEODB" -q "$faa" -o "$hits" \
    -e 1e-10 -k 5 --query-cover 50 --subject-cover 50 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
    > "EET_check/logs/${b}.diamond.log" 2>&1

  echo "[DONE] $b"
done

echo "=============================="
echo "[ALL DONE] Outputs are under: $(pwd)/EET_check/"

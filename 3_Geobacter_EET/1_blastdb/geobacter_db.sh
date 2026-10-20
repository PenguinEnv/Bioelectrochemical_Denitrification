#!/bin/bash
#
# Requires: wget, DIAMOND

set -euo pipefail

mkdir -p geobacter_ref && cd geobacter_ref

download_protein_faa () {
  local acc="$1"      # e.g., GCF_000007985.2
  local asm="$2"      # e.g., ASM798v2
  local digits
  digits=$(echo "$acc" | sed -E 's/^GCF_([0-9]{3})([0-9]{3})([0-9]{3})\..*/\1\/\2\/\3/')

  local prefix="${acc}_${asm}"
  local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${digits}/${prefix}/${prefix}_protein.faa.gz"

  echo "[INFO] $acc ($asm)"
  echo "  -> $url"
  wget -c "$url"
  gunzip -f "${prefix}_protein.faa.gz"
}

# three key Geobacter species
download_protein_faa "GCF_000007985.2" "ASM798v2"     # G. sulfurreducens PCA
download_protein_faa "GCF_000012925.1" "ASM1292v1"    # G. metallireducens GS-15
download_protein_faa "GCF_000816575.1" "ASM81657v1"   # G. anodireducens GSS01

# merge and build the DIAMOND reference
cat *_protein.faa > Geobacter_major_ref_proteins.faa

echo "[INFO] Total proteins fasta:"
ls -lh Geobacter_major_ref_proteins.faa

diamond makedb \
  --in Geobacter_major_ref_proteins.faa \
  -d major_ref

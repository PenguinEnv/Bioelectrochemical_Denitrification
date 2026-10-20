#!/bin/bash
#
# Requires: bash, awk
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

ROOT="/path/to/base/3_Geobacter_EET/2_hmmer_diamond"
OUTDIR="${ROOT}/EET_organize/reports"
PFAM_DIR="${ROOT}/EET_check/pfam"
DIA_DIR="${ROOT}/EET_check/diamond"
PROT_DIR="${ROOT}/EET_check/prodigal"

# Geobacter reference FASTA (contains WP IDs in headers)
REF_FASTA="/path/to/base/3_Geobacter_EET/1_blastdb/geobacter_ref/Geobacter_major_ref_proteins.faa"

BINS=(semibin_173 maxbin2_085 concoct_bin_89_sub)

# DIAMOND base filters
EVALUE_CUTOFF="1e-10"
MIN_QCOV=50
MIN_SCOV=50

# multiheme threshold
MIN_CXXCH=3

mkdir -p "${OUTDIR}"

echo "[INFO] OUTDIR: ${OUTDIR}"
echo "[INFO] REF_FASTA: ${REF_FASTA}"
echo "[INFO] DIAMOND: evalue<${EVALUE_CUTOFF}, qcov>=${MIN_QCOV}, scov>=${MIN_SCOV}"
echo "[INFO] multiheme: CXXCH>=${MIN_CXXCH}"
echo

# 1) Build EET-like subject list from reference FASTA headers (WP IDs whose descriptions look EET-ish)
# You can edit keywords if needed.
EET_SUBJECTS="${OUTDIR}/EET_subjects_from_ref.txt"
grep "^>" "${REF_FASTA}" \
  | grep -iE "cytochrome|multiheme|heme|pilin|type[[:space:]]*iv|exta|extb|extc|extd|imch|cbcl|porin|mtr" \
  | grep -oE "WP_[0-9]+\.[0-9]+" \
  | sort -u > "${EET_SUBJECTS}"

echo "[INFO] EET-like subjects in reference: $(wc -l < "${EET_SUBJECTS}")"
echo

SUMMARY="${OUTDIR}/summary_counts.tsv"
echo -e "bin\tpfam_ids\tdiamond_all_ids\tdiamond_eet_ids\tunion_ids\tmultiheme_CXXCH_ge3" > "${SUMMARY}"

for b in "${BINS[@]}"; do
  echo "=============================="
  echo "[INFO] Processing bin: ${b}"

  PFAM_FILE="${PFAM_DIR}/${b}.pfam.domtblout"
  DIA_FILE="${DIA_DIR}/${b}.geobacter_ref.tsv"
  FAA_FILE="${PROT_DIR}/${b}.faa"

  [[ -f "${PFAM_FILE}" ]] || { echo "[ERROR] Missing: ${PFAM_FILE}"; exit 1; }
  [[ -f "${DIA_FILE}"  ]] || { echo "[ERROR] Missing: ${DIA_FILE}";  exit 1; }
  [[ -f "${FAA_FILE}"  ]] || { echo "[ERROR] Missing: ${FAA_FILE}";  exit 1; }

  # A) PFAM candidates (fix awk to one-liner, robust on $0)
  awk 'BEGIN{IGNORECASE=1} $0!~/^#/ && ($1~/cytochrome/ || $0~/cytochrome/ || $1~/pilin/ || $0~/pilin/) {print $4}' \
    "${PFAM_FILE}" | sort -u > "${OUTDIR}/${b}.EET_from_PFAM.ids"

  # B1) DIAMOND candidates (all strong Geobacter homologs)
  awk -v ev="${EVALUE_CUTOFF}" -v q="${MIN_QCOV}" -v s="${MIN_SCOV}" '
    function toNum(x){ return x+0 }
    (toNum($5) < toNum(ev)) && ($7+0 >= q) && ($8+0 >= s) {print $1}
  ' "${DIA_FILE}" | sort -u > "${OUTDIR}/${b}.DIAMOND_strong.ids"

  # B2) DIAMOND EET-like candidates: strong hits whose SUBJECT is in EET_SUBJECTS
  awk -v ev="${EVALUE_CUTOFF}" -v q="${MIN_QCOV}" -v s="${MIN_SCOV}" '
    function toNum(x){ return x+0 }
    (toNum($5) < toNum(ev)) && ($7+0 >= q) && ($8+0 >= s) {print $1"\t"$2}
  ' "${DIA_FILE}" > "${OUTDIR}/${b}.diamond.filtered.tmp"

  # keep only those with EET-like subjects, then output unique query IDs
  awk 'NR==FNR{a[$1]=1; next} ($2 in a){print $1}' \
    "${EET_SUBJECTS}" "${OUTDIR}/${b}.diamond.filtered.tmp" \
    | sort -u > "${OUTDIR}/${b}.EET_from_DIAMOND.ids"

  rm -f "${OUTDIR}/${b}.diamond.filtered.tmp"

  # C) UNION (PFAM + DIAMOND_EET)
  cat "${OUTDIR}/${b}.EET_from_PFAM.ids" "${OUTDIR}/${b}.EET_from_DIAMOND.ids" \
    | sort -u > "${OUTDIR}/${b}.EET_candidates.UNION.ids"

  # D) multiheme motif counting from proteins
  awk -v min="${MIN_CXXCH}" '
    function count_cxxch(s,   c){
      c=0
      while (match(s, /C..CH/)) { c++; s=substr(s, RSTART+1) }
      return c
    }
    /^>/{
      if (seq!="") {
        c=count_cxxch(seq)
        if (c>=min) print id"\t"c
      }
      id=substr($1,2); seq=""; next
    }
    { gsub(/[ \t\r]/,""); seq=seq$0 }
    END{
      if (seq!="") {
        c=count_cxxch(seq)
        if (c>=min) print id"\t"c
      }
    }
  ' "${FAA_FILE}" | sort -k2,2nr > "${OUTDIR}/${b}.multiheme_CXXCH_ge3.tsv"

  PFAM_N=$(wc -l < "${OUTDIR}/${b}.EET_from_PFAM.ids" | tr -d ' ')
  DIA_ALL_N=$(wc -l < "${OUTDIR}/${b}.DIAMOND_strong.ids" | tr -d ' ')
  DIA_EET_N=$(wc -l < "${OUTDIR}/${b}.EET_from_DIAMOND.ids" | tr -d ' ')
  UNI_N=$(wc -l < "${OUTDIR}/${b}.EET_candidates.UNION.ids" | tr -d ' ')
  MHE_N=$(wc -l < "${OUTDIR}/${b}.multiheme_CXXCH_ge3.tsv" | tr -d ' ')

  echo -e "${b}\t${PFAM_N}\t${DIA_ALL_N}\t${DIA_EET_N}\t${UNI_N}\t${MHE_N}" >> "${SUMMARY}"

  echo "[DONE] ${b}  PFAM=${PFAM_N}  DIAMOND_all=${DIA_ALL_N}  DIAMOND_EET=${DIA_EET_N}  UNION=${UNI_N}  multiheme=${MHE_N}"
done

echo "=============================="
echo "[ALL DONE] Summary: ${SUMMARY}"
column -t "${SUMMARY}"

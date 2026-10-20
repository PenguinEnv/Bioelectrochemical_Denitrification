#!/bin/bash
# Requires: HMMER 3.4; Pfam-A.hmm
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

PFAM_DIR="/path/to/your_db/Pfam"
PFAM_HMM="${PFAM_DIR}/Pfam-A.hmm"

OUTDIR="/path/to/base/3_Geobacter_EET/3_community_EET/Pfam_eet"
mkdir -p "${OUTDIR}"

# 1) keyword set
KW='cytochrome|c-type|multiheme|decaheme|mtr|omc|exta|extb|extc|extd|pilin|pilus|type[[:space:]]*iv|T4P'

NAME_LIST="${OUTDIR}/EET_pfam_names.txt"
SUB_HMM="${OUTDIR}/EET_Pfam_subset.hmm"

echo "[INFO] Scanning Pfam-A.hmm to collect NAMES..."
awk -v IGNORECASE=1 -v kw="$KW" '
  /^HMMER3/ {name=""}
  /^NAME[ \t]+/ {name=$2}
  /^DESC[ \t]+/ { if (name!="" && $0 ~ kw) print name }
' "${PFAM_HMM}" | sort -u > "${NAME_LIST}"

echo "[DONE] $(wc -l < "${NAME_LIST}") NAMES -> ${NAME_LIST}"
echo "[INFO] Example NAMES:"
head -n 20 "${NAME_LIST}"

# 2) build index (note: the index file is written next to PFAM_HMM, i.e. in PFAM_DIR)
echo "[INFO] hmmfetch --index ${PFAM_HMM}"
hmmfetch --index "${PFAM_HMM}"

# 3) extract subset
echo "[INFO] hmmfetch -f to create subset HMM..."
hmmfetch -f "${PFAM_HMM}" "${NAME_LIST}" > "${SUB_HMM}"

# 4) hmmpress for hmmscan (pressed files are written next to SUB_HMM, i.e. in OUTDIR)
echo "[INFO] hmmpress ${SUB_HMM}"
hmmpress "${SUB_HMM}"

echo "[DONE] subset HMM created:"
ls -lh "${SUB_HMM}" "${SUB_HMM}".h3*

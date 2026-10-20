#!/bin/bash
# Requires: bash, awk
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

############################
# INPUTS
############################
ROOT_OUT="/path/to/base/3_Geobacter_EET/3_community_EET"
EET_OUT="${ROOT_OUT}/2_eet_define"
SALMON_OUT="${ROOT_OUT}/3_salmon"
REPORT_OUT="${ROOT_OUT}/4_reports"

MAG_DIR="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"

# target MAGs
MAGS=("concoct_bin_121_sub" "concoct_bin_124" "concoct_bin_164" "concoct_bin_34_sub" "concoct_bin_78_sub" "concoct_bin_89_sub" "concoct_bin_95_sub" "maxbin2_047" "maxbin2_085" "maxbin2_105_sub" "semibin_144" "semibin_151" "semibin_173" "semibin_190" "semibin_83" "vamb_13")

# DNA / RNA sample names
DNA_SAMPLES=("zr15207_10" "zr15207_11" "zr15207_13" "zr15207_14" "zr15207_9")
RNA_SAMPLES=("zr15207_2T" "zr15207_3T" "zr15207_4T" "zr15207_6T" "zr15207_7T")

mkdir -p "${REPORT_OUT}/MAG_contributions"

############################
# HELPERS
############################
die(){ echo "[ERROR] $*" >&2; exit 1; }
check_file(){ [[ -f "$1" ]] || die "Missing file: $1"; }

find_fasta(){
  local base="$1"
  for ext in fa fasta fna fas; do
    if [[ -f "${MAG_DIR}/${base}.${ext}" ]]; then
      echo "${MAG_DIR}/${base}.${ext}"
      return 0
    fi
  done
  if [[ -f "${MAG_DIR}/${base}" ]]; then
    echo "${MAG_DIR}/${base}"
    return 0
  fi
  return 1
}

sum_tpm(){
  local gene_list="$1"
  local quant_sf="$2"
  awk 'NR==FNR{a[$1]=1; next} ($1 in a){sum+=$4} END{printf "%.6f", sum}' \
    "${gene_list}" "${quant_sf}"
}

pct(){
  local num="$1"
  local den="$2"
  awk -v n="$num" -v d="$den" 'BEGIN{ if(d==0) print "NA"; else printf "%.3f", (n/d)*100 }'
}

############################
# EXISTING EET GENE SETS
############################
EET_BROAD="${EET_OUT}/lists/EET_genes.BROAD.ids"
EET_HIGH="${EET_OUT}/lists/EET_genes.HIGHCONF.ids"
EET_BROAD_LOC="${EET_OUT}/localization/EET_genes.BROAD.localized.ids"
EET_HIGH_LOC="${EET_OUT}/localization/EET_genes.HIGHCONF.localized.ids"

check_file "${EET_BROAD}"
check_file "${EET_HIGH}"

############################
# BUILD MAG-SPECIFIC GENE LISTS
############################
for MAG in "${MAGS[@]}"; do
  echo "[INFO] Processing MAG: ${MAG}"

  MAG_FASTA=$(find_fasta "${MAG}") || die "FASTA not found for ${MAG} in ${MAG_DIR}"
  echo "[INFO] FASTA: ${MAG_FASTA}"

  MAG_CONTIGS="${REPORT_OUT}/MAG_contributions/${MAG}.contigs.txt"
  MAG_EET_B="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.genes"
  MAG_EET_H="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.genes"
  MAG_EET_B_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.localized.genes"
  MAG_EET_H_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.localized.genes"

  grep "^>" "${MAG_FASTA}" | sed 's/^>//' | awk '{print $1}' | sort -u > "${MAG_CONTIGS}"

  awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if(c in a) print $0}' \
    "${MAG_CONTIGS}" "${EET_BROAD}" | sort -u > "${MAG_EET_B}"

  awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if(c in a) print $0}' \
    "${MAG_CONTIGS}" "${EET_HIGH}" | sort -u > "${MAG_EET_H}"

  if [[ -s "${EET_BROAD_LOC}" ]]; then
    awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if(c in a) print $0}' \
      "${MAG_CONTIGS}" "${EET_BROAD_LOC}" | sort -u > "${MAG_EET_B_LOC}"
  else
    : > "${MAG_EET_B_LOC}"
  fi

  if [[ -s "${EET_HIGH_LOC}" ]]; then
    awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if(c in a) print $0}' \
      "${MAG_CONTIGS}" "${EET_HIGH_LOC}" | sort -u > "${MAG_EET_H_LOC}"
  else
    : > "${MAG_EET_H_LOC}"
  fi

  echo "[INFO] ${MAG}: BROAD=$(wc -l < "${MAG_EET_B}") HIGH=$(wc -l < "${MAG_EET_H}")"
done

############################
# DNA REPORT
############################
DNA_OUT="${REPORT_OUT}/MAG_contributions/DNA_MAG_contributions.tsv"
echo -e "sample\tMAG\tmag_EET_BROAD\tpct_of_total_BROAD\tmag_EET_HIGH\tpct_of_total_HIGH\tmag_LOC_BROAD\tpct_of_total_LOC_BROAD\tmag_LOC_HIGH\tpct_of_total_LOC_HIGH" > "${DNA_OUT}"

for s in "${DNA_SAMPLES[@]}"; do
  q="${SALMON_OUT}/DNA/${s}_salmon/quant.sf"
  [[ -f "${q}" ]] || { echo "[WARN] Missing DNA quant for ${s}: ${q}"; continue; }

  total_b=$(sum_tpm "${EET_BROAD}" "${q}")
  total_h=$(sum_tpm "${EET_HIGH}" "${q}")

  total_lb="0.000000"
  total_lh="0.000000"
  [[ -s "${EET_BROAD_LOC}" ]] && total_lb=$(sum_tpm "${EET_BROAD_LOC}" "${q}")
  [[ -s "${EET_HIGH_LOC}"  ]] && total_lh=$(sum_tpm "${EET_HIGH_LOC}"  "${q}")

  for MAG in "${MAGS[@]}"; do
    MAG_EET_B="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.genes"
    MAG_EET_H="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.genes"
    MAG_EET_B_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.localized.genes"
    MAG_EET_H_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.localized.genes"

    mag_b=$(sum_tpm "${MAG_EET_B}" "${q}")
    mag_h=$(sum_tpm "${MAG_EET_H}" "${q}")

    mag_lb="0.000000"
    mag_lh="0.000000"
    [[ -s "${MAG_EET_B_LOC}" ]] && mag_lb=$(sum_tpm "${MAG_EET_B_LOC}" "${q}")
    [[ -s "${MAG_EET_H_LOC}" ]] && mag_lh=$(sum_tpm "${MAG_EET_H_LOC}" "${q}")

    pct_b=$(pct "${mag_b}" "${total_b}")
    pct_h=$(pct "${mag_h}" "${total_h}")
    pct_lb=$(pct "${mag_lb}" "${total_lb}")
    pct_lh=$(pct "${mag_lh}" "${total_lh}")

    echo -e "${s}\t${MAG}\t${mag_b}\t${pct_b}\t${mag_h}\t${pct_h}\t${mag_lb}\t${pct_lb}\t${mag_lh}\t${pct_lh}" >> "${DNA_OUT}"
  done
done

############################
# RNA REPORT
############################
RNA_OUT="${REPORT_OUT}/MAG_contributions/RNA_MAG_contributions.tsv"
echo -e "sample\tMAG\tmag_EET_BROAD\tpct_of_total_BROAD\tmag_EET_HIGH\tpct_of_total_HIGH\tmag_LOC_BROAD\tpct_of_total_LOC_BROAD\tmag_LOC_HIGH\tpct_of_total_LOC_HIGH" > "${RNA_OUT}"

for s in "${RNA_SAMPLES[@]}"; do
  q="${SALMON_OUT}/RNA/${s}_salmon/quant.sf"
  [[ -f "${q}" ]] || { echo "[WARN] Missing RNA quant for ${s}: ${q}"; continue; }

  total_b=$(sum_tpm "${EET_BROAD}" "${q}")
  total_h=$(sum_tpm "${EET_HIGH}" "${q}")

  total_lb="0.000000"
  total_lh="0.000000"
  [[ -s "${EET_BROAD_LOC}" ]] && total_lb=$(sum_tpm "${EET_BROAD_LOC}" "${q}")
  [[ -s "${EET_HIGH_LOC}"  ]] && total_lh=$(sum_tpm "${EET_HIGH_LOC}"  "${q}")

  for MAG in "${MAGS[@]}"; do
    MAG_EET_B="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.genes"
    MAG_EET_H="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.genes"
    MAG_EET_B_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_BROAD.localized.genes"
    MAG_EET_H_LOC="${REPORT_OUT}/MAG_contributions/${MAG}.EET_HIGHCONF.localized.genes"

    mag_b=$(sum_tpm "${MAG_EET_B}" "${q}")
    mag_h=$(sum_tpm "${MAG_EET_H}" "${q}")

    mag_lb="0.000000"
    mag_lh="0.000000"
    [[ -s "${MAG_EET_B_LOC}" ]] && mag_lb=$(sum_tpm "${MAG_EET_B_LOC}" "${q}")
    [[ -s "${MAG_EET_H_LOC}" ]] && mag_lh=$(sum_tpm "${MAG_EET_H_LOC}" "${q}")

    pct_b=$(pct "${mag_b}" "${total_b}")
    pct_h=$(pct "${mag_h}" "${total_h}")
    pct_lb=$(pct "${mag_lb}" "${total_lb}")
    pct_lh=$(pct "${mag_lh}" "${total_lh}")

    echo -e "${s}\t${MAG}\t${mag_b}\t${pct_b}\t${mag_h}\t${pct_h}\t${mag_lb}\t${pct_lb}\t${mag_lh}\t${pct_lh}" >> "${RNA_OUT}"
  done
done

echo "[DONE] Output files:"
echo "  - ${DNA_OUT}"
echo "  - ${RNA_OUT}"

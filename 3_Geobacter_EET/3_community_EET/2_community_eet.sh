#!/bin/bash
# Requires: Prodigal 2.6.3, HMMER 3.4, DIAMOND, Salmon 1.10.3, Python 3
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.

###############################################################################
# Community-wide EET pipeline + Geobacter contribution (DNA abundance + RNA TPM)
#
# Steps:
# 1) Prodigal on assembly contigs -> gene catalog (CDS + proteins + GFF)
# 2) Define EET genes:
#    - PFAM subset hmmscan (--cut_ga)  [FAST]
#    - multiheme motif (CXXCH >= MIN_CXXCH)
#    - optional DIAMOND (EET-like)
#    - build EET gene sets: BROAD / HIGHCONF
#    - localization (genomic clustering on contigs) -> localized gene lists
# 3) Salmon index on CDS catalog
# 4) DNA quant -> TPM-like abundance
# 5) RNA quant -> TPM
# 6) Reports:
#    - total EET vs Geobacter EET (%)
#    - localized EET vs Geobacter localized EET (%)
###############################################################################

############################
# INPUTS                   #
############################

ASSEMBLY_CONTIGS="/path/to/base/1_metagenomic/2_assemble/megahit_out/final.contigs.fa"

GEO_BIN_DIR="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"
GEO_BINS=("semibin_173.fa" "maxbin2_085.fa" "concoct_bin_89_sub.fa")

# Pfam subset HMM (233 models) + pressed files already exist
PFAM_HMM="/path/to/base/3_Geobacter_EET/3_community_EET/Pfam_eet/EET_Pfam_subset.hmm"

# Optional DIAMOND reference
USE_DIAMOND=1
GEO_REF_FASTA="/path/to/base/3_Geobacter_EET/1_blastdb/geobacter_ref/Geobacter_major_ref_proteins.faa"
GEO_REF_DMND="/path/to/base/3_Geobacter_EET/1_blastdb/geobacter_ref/Geobacter_major_ref_proteins"

# RNA reads
RNA_READ_DIR="/path/to/base/2_metatranscriptomic/1_qc/2_sortmerna"
RNA_SAMPLES=("zr15207_2T" "zr15207_3T" "zr15207_4T" "zr15207_6T" "zr15207_7T")
RNA_READ_SUFFIX="_mRNA_SILVA.fq.gz"

# DNA reads
DNA_MODE="PE"   # "PE" or "SE"
DNA_READ_DIR="/path/to/base/1_metagenomic/1_qc/1_fastp_out"
DNA_SAMPLES=("zr15207_10" "zr15207_11" "zr15207_13" "zr15207_14" "zr15207_9")
DNA_R1_SUFFIX="_clean_R1.fastq.gz"
DNA_R2_SUFFIX="_clean_R2.fastq.gz"
DNA_SE_SUFFIX="_clean.fastq.gz"

# Threads
THREADS=16

# EET definition params
MIN_CXXCH=3
DIAMOND_EVALUE="1e-10"
DIAMOND_MIN_QCOV=50
DIAMOND_MIN_SCOV=50
EET_REF_KEYWORDS_REGEX="cytochrome|multiheme|heme|pilin|type[[:space:]]*iv|exta|extb|extc|extd|imch|cbcl|porin|mtr|omc"

# Localization (genomic clustering) params
CLUSTER_MAX_BP=10000      # EET genes <= 10 kb apart are grouped into one locus
CLUSTER_MIN_GENES=2       # minimum number of EET genes per locus

############################
# OUTPUT DIRS              #
############################
ROOT_OUT="/path/to/base/3_Geobacter_EET/3_community_EET"
PRODIGAL_OUT="${ROOT_OUT}/1_prodigal"
EET_OUT="${ROOT_OUT}/2_eet_define"
SALMON_OUT="${ROOT_OUT}/3_salmon"
REPORT_OUT="${ROOT_OUT}/4_reports"

mkdir -p "${PRODIGAL_OUT}"/{proteins,cds,gff,logs} \
         "${EET_OUT}"/{pfam,diamond,lists,localization,logs} \
         "${SALMON_OUT}"/{index,DNA,RNA,logs} \
         "${REPORT_OUT}"

############################
# HELPER FUNCTIONS         #
############################
die(){ echo "[ERROR] $*" >&2; exit 1; }
check_file(){ [[ -f "$1" ]] || die "Missing file: $1"; }
gene_to_contig(){ sed 's/_[0-9]\+$//'; }

############################
# STEP 1: Prodigal         #
############################
step1_prodigal(){
  echo "============================================================"
  echo "[STEP 1] Prodigal on assembly contigs -> gene catalog: $(date '+%Y-%m-%d %H:%M:%S')"
  check_file "${ASSEMBLY_CONTIGS}"

  local FAA="${PRODIGAL_OUT}/proteins/assembly.proteins.faa"
  local FNA="${PRODIGAL_OUT}/cds/assembly.genes.fna"
  local GFF="${PRODIGAL_OUT}/gff/assembly.prodigal.gff"

  if [[ -s "${FAA}" && -s "${FNA}" && -s "${GFF}" ]]; then
    echo "[INFO] Prodigal outputs already exist. Skipping."
    echo "[INFO] ${FAA}"
    echo "[INFO] ${FNA}"
    echo "[INFO] ${GFF}"
    return 0
  fi

  prodigal -i "${ASSEMBLY_CONTIGS}" \
    -a "${FAA}" \
    -d "${FNA}" \
    -o "${GFF}" -f gff \
    -p meta \
    > "${PRODIGAL_OUT}/logs/prodigal.log" 2>&1

  echo "[DONE] Prodigal outputs:"
  echo "  - ${FAA}"
  echo "  - ${FNA}"
  echo "  - ${GFF}"
}

############################
# STEP 2: EET define + localization
############################
step2_define_eet(){
  echo "============================================================"
  echo "[STEP 2] Define EET genes + localization: $(date '+%Y-%m-%d %H:%M:%S')"

  local FAA="${PRODIGAL_OUT}/proteins/assembly.proteins.faa"
  local GFF="${PRODIGAL_OUT}/gff/assembly.prodigal.gff"
  check_file "${FAA}"
  check_file "${GFF}"
  check_file "${PFAM_HMM}"

  local DOMTBL="${EET_OUT}/pfam/assembly.EETsubset.domtblout"

  ############################################
  # 2A) PFAM scan (subset HMM) on SCRATCH
  ############################################
  if [[ ! -s "${DOMTBL}" ]]; then
    echo "[INFO] Running hmmscan (subset Pfam) on scratch..."

    local JOBID="${PBS_JOBID:-${SLURM_JOB_ID:-$$}}"
    local SCRATCH_BASE="/scratch/${USER}"
    mkdir -p "${SCRATCH_BASE}"
    local SCRATCH_DIR="${SCRATCH_BASE}/pfam_subset_${JOBID}"
    mkdir -p "${SCRATCH_DIR}"

    local PFAM_DIR
    PFAM_DIR="$(dirname "${PFAM_HMM}")"
    local PFAM_BASE
    PFAM_BASE="$(basename "${PFAM_HMM}")"

    cp "${FAA}" "${SCRATCH_DIR}/" || die "Failed to copy FAA to scratch"

    # subset HMM + pressed files copy
    for f in "${PFAM_BASE}" "${PFAM_BASE}.h3f" "${PFAM_BASE}.h3i" "${PFAM_BASE}.h3m" "${PFAM_BASE}.h3p"; do
      [[ -f "${PFAM_DIR}/${f}" ]] || die "Missing pressed file: ${PFAM_DIR}/${f}"
      cp "${PFAM_DIR}/${f}" "${SCRATCH_DIR}/" || die "Failed to copy: ${f}"
    done

    local ORIG_DIR
    ORIG_DIR=$(pwd)
    cd "${SCRATCH_DIR}"

    # --cut_ga: Pfam curated gathering thresholds
    hmmscan --cpu "${THREADS}" --cut_ga \
      --domtblout assembly.EETsubset.domtblout \
      "${PFAM_BASE}" \
      "$(basename "${FAA}")" \
      > hmmscan.log 2>&1 || die "hmmscan failed (see ${SCRATCH_DIR}/hmmscan.log)"

    [[ -s assembly.EETsubset.domtblout ]] || die "PFAM domtblout missing"
    cp assembly.EETsubset.domtblout "${DOMTBL}" || die "Copy back failed"

    cd "${ORIG_DIR}"
    rm -rf "${SCRATCH_DIR}"
  else
    echo "[INFO] PFAM domtblout exists. Skipping hmmscan."
  fi

  ############################################
  # 2B) PFAM hit genes -> PFAM_IDS
  ############################################
  local PFAM_IDS="${EET_OUT}/lists/EET_from_PFAM.ids"
  awk '$0!~/^#/ {print $4}' "${DOMTBL}" | sort -u > "${PFAM_IDS}"

  ############################################
  # 2C) multiheme detection (CXXCH >= MIN_CXXCH)
  ############################################
  local MHE_TSV="${EET_OUT}/lists/multiheme_CXXCH_ge${MIN_CXXCH}.tsv"
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
  ' "${FAA}" | sort -k2,2nr > "${MHE_TSV}"
  cut -f1 "${MHE_TSV}" | sort -u > "${EET_OUT}/lists/multiheme.ids"

  ############################################
  # 2D) DIAMOND filtering (optional)
  ############################################
  local DIA_EET_IDS="${EET_OUT}/lists/EET_from_DIAMOND.ids"
  local DIA_TSV="${EET_OUT}/diamond/assembly_vs_georef.tsv"
  local EET_SUBJECTS="${EET_OUT}/lists/EET_subjects_from_georef.txt"

  if [[ "${USE_DIAMOND}" -eq 1 ]]; then
    check_file "${GEO_REF_FASTA}"

    awk -v IGNORECASE=1 -v re="${EET_REF_KEYWORDS_REGEX}" '
      /^>/{
        hdr=substr($0,2)
        id=hdr; sub(/[ \t].*$/,"",id)
        if (hdr ~ re) print id
      }
    ' "${GEO_REF_FASTA}" | sort -u > "${EET_SUBJECTS}"

    if [[ ! -f "${GEO_REF_DMND}.dmnd" ]]; then
      diamond makedb --in "${GEO_REF_FASTA}" -d "${GEO_REF_DMND}" \
        > "${EET_OUT}/logs/diamond_makedb.log" 2>&1 || die "diamond makedb failed"
    fi

    if [[ ! -s "${DIA_TSV}" ]]; then
      diamond blastp -d "${GEO_REF_DMND}" -q "${FAA}" -o "${DIA_TSV}" \
        -e "${DIAMOND_EVALUE}" -k 5 \
        --query-cover "${DIAMOND_MIN_QCOV}" \
        --subject-cover "${DIAMOND_MIN_SCOV}" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp \
        > "${EET_OUT}/logs/diamond_blastp.log" 2>&1 || die "DIAMOND failed"
    fi

    awk 'NR==FNR{a[$1]=1; next} ($2 in a){print $1}' \
      "${EET_SUBJECTS}" "${DIA_TSV}" | sort -u > "${DIA_EET_IDS}"
  else
    : > "${DIA_EET_IDS}"
  fi

  ############################################
  # 2E) Final gene sets
  ############################################
  cat "${PFAM_IDS}" "${DIA_EET_IDS}" "${EET_OUT}/lists/multiheme.ids" \
    | sort -u > "${EET_OUT}/lists/EET_genes.BROAD.ids"

  if [[ "${USE_DIAMOND}" -eq 1 ]]; then
    comm -12 <(sort "${PFAM_IDS}") <(sort "${DIA_EET_IDS}") \
      > "${EET_OUT}/lists/EET_genes.PFAM_AND_DIAMOND.ids"
    cat "${EET_OUT}/lists/EET_genes.PFAM_AND_DIAMOND.ids" "${EET_OUT}/lists/multiheme.ids" \
      | sort -u > "${EET_OUT}/lists/EET_genes.HIGHCONF.ids"
  else
    cp "${EET_OUT}/lists/multiheme.ids" "${EET_OUT}/lists/EET_genes.HIGHCONF.ids"
  fi

  ############################################
  # 2F) Localization (genomic clustering on contigs)
  ############################################
  echo "[INFO] Building localization (EET loci by contig clustering)..."

  # Make gene coordinate table from Prodigal GFF
  # output: gene_id  contig  start  end  strand
  local GENE_COORDS="${EET_OUT}/localization/gene_coords.tsv"
  awk -F'\t' '
    BEGIN{OFS="\t"}
    $0 ~ /^#/ {next}
    $3=="CDS" {
      contig=$1; start=$4; end=$5; strand=$7
      id=""
      n=split($9,a,";")
      for(i=1;i<=n;i++){
        if(a[i] ~ /^ID=/){ id=a[i]; sub(/^ID=/,"",id) }
      }
      if(id!="") print id, contig, start, end, strand
    }
  ' "${GFF}" > "${GENE_COORDS}"

  # Define localized genes for BROAD + HIGHCONF
  for SET in BROAD HIGHCONF; do
    local EET_LIST="${EET_OUT}/lists/EET_genes.${SET}.ids"
    local OUT_LOC="${EET_OUT}/localization/EET_genes.${SET}.localized.ids"
    local OUT_LOCI_TSV="${EET_OUT}/localization/EET_loci.${SET}.tsv"

    check_file "${EET_LIST}"

    # Join coords with EET gene list, sort by contig+start
    # fields: contig start end gene strand
    local TMP_SORT="${EET_OUT}/localization/tmp.${SET}.sorted.tsv"
    awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2"\t"$3"\t"$4"\t"$1"\t"$5}' \
      "${EET_LIST}" "${GENE_COORDS}" \
      | sort -k1,1 -k2,2n > "${TMP_SORT}"

    # Cluster within contig by distance between consecutive EET genes
    # locus_id = contig|L<n>
    awk -v maxbp="${CLUSTER_MAX_BP}" -v minGenes="${CLUSTER_MIN_GENES}" '
      BEGIN{OFS="\t"}
      {
        contig=$1; s=$2; e=$3; gene=$4
        if(NR==1){
          cur_contig=contig; locus=1
          last_end=e
          genes=gene
          n=1
          start_min=s; end_max=e
          next
        }
        if(contig!=cur_contig){
          # finalize previous contig locus
          if(n>=minGenes){
            split(genes,garr,",")
            for(i in garr) print garr[i] > out_ids
            print cur_contig"|L"locus, cur_contig, start_min, end_max, n, genes > out_loci
          }
          # reset
          cur_contig=contig; locus=1
          last_end=e
          genes=gene
          n=1
          start_min=s; end_max=e
          next
        }
        # same contig
        if(s - last_end <= maxbp){
          genes=genes","gene
          n++
          if(e>end_max) end_max=e
          last_end=e
        } else {
          # close previous locus
          if(n>=minGenes){
            split(genes,garr,",")
            for(i in garr) print garr[i] > out_ids
            print cur_contig"|L"locus, cur_contig, start_min, end_max, n, genes > out_loci
          }
          locus++
          # start new locus
          last_end=e
          genes=gene
          n=1
          start_min=s; end_max=e
        }
      }
      END{
        if(n>=minGenes){
          split(genes,garr,",")
          for(i in garr) print garr[i] > out_ids
          print cur_contig"|L"locus, cur_contig, start_min, end_max, n, genes > out_loci
        }
      }
    ' out_ids="${OUT_LOC}" out_loci="${OUT_LOCI_TSV}" "${TMP_SORT}"

    # de-dup gene list
    sort -u "${OUT_LOC}" -o "${OUT_LOC}" 2>/dev/null || true

    echo "[INFO] ${SET}: localized genes = $(wc -l < "${OUT_LOC}" 2>/dev/null || echo 0)"
    echo "[INFO] ${SET}: loci table -> ${OUT_LOCI_TSV}"
  done

  ############################################
  # 2G) Summary counts
  ############################################
  {
    echo -e "list\tcount"
    echo -e "PFAM_subset_hits\t$(wc -l < "${PFAM_IDS}")"
    echo -e "DIAMOND_EET\t$(wc -l < "${DIA_EET_IDS}")"
    echo -e "multiheme_CXXCH_ge${MIN_CXXCH}\t$(wc -l < "${EET_OUT}/lists/multiheme.ids")"
    echo -e "EET_BROAD\t$(wc -l < "${EET_OUT}/lists/EET_genes.BROAD.ids")"
    echo -e "EET_HIGHCONF\t$(wc -l < "${EET_OUT}/lists/EET_genes.HIGHCONF.ids")"
    echo -e "EET_BROAD_localized\t$(wc -l < "${EET_OUT}/localization/EET_genes.BROAD.localized.ids" 2>/dev/null || echo 0)"
    echo -e "EET_HIGHCONF_localized\t$(wc -l < "${EET_OUT}/localization/EET_genes.HIGHCONF.localized.ids" 2>/dev/null || echo 0)"
  } > "${EET_OUT}/lists/EET_define_counts.tsv"

  echo "[DONE] STEP2 completed: $(date '+%Y-%m-%d %H:%M:%S')"
}

############################
# STEP 3: Salmon index     #
############################
step_salmon_index(){
  echo "============================================================"
  echo "[STEP 3] Build Salmon index on CDS gene catalog: $(date '+%Y-%m-%d %H:%M:%S')"

  local FNA="${PRODIGAL_OUT}/cds/assembly.genes.fna"
  check_file "${FNA}"

  local IDX="${SALMON_OUT}/index/CDS_index"
  if [[ -d "${IDX}" && -s "${IDX}/hash.bin" ]]; then
    echo "[INFO] Salmon index exists. Skipping."
    return 0
  fi

  salmon index -t "${FNA}" -i "${IDX}" -k 31 \
    > "${SALMON_OUT}/logs/salmon_index.log" 2>&1

  echo "[DONE] Salmon index: ${IDX}"
}

############################
# STEP 4: DNA quant        #
############################
step3_dna_quant(){
  echo "============================================================"
  echo "[STEP 4] DNA quant (metagenome -> TPM-like): $(date '+%Y-%m-%d %H:%M:%S')"

  local IDX="${SALMON_OUT}/index/CDS_index"
  [[ -d "${IDX}" ]] || die "Missing Salmon index."

  mkdir -p "${SALMON_OUT}/DNA"

  for s in "${DNA_SAMPLES[@]}"; do
    echo "[INFO] DNA sample: ${s}"
    local OUT="${SALMON_OUT}/DNA/${s}_salmon"
    if [[ -s "${OUT}/quant.sf" ]]; then
      echo "  [SKIP] quant.sf exists"
      continue
    fi

    if [[ "${DNA_MODE}" == "PE" ]]; then
      local R1="${DNA_READ_DIR}/${s}${DNA_R1_SUFFIX}"
      local R2="${DNA_READ_DIR}/${s}${DNA_R2_SUFFIX}"
      check_file "${R1}"; check_file "${R2}"
      salmon quant -i "${IDX}" -l A -1 "${R1}" -2 "${R2}" \
        -p "${THREADS}" --validateMappings -o "${OUT}" \
        > "${SALMON_OUT}/logs/${s}.DNA.salmon.log" 2>&1
    else
      local R="${DNA_READ_DIR}/${s}${DNA_SE_SUFFIX}"
      check_file "${R}"
      salmon quant -i "${IDX}" -l A -r "${R}" \
        -p "${THREADS}" --validateMappings -o "${OUT}" \
        > "${SALMON_OUT}/logs/${s}.DNA.salmon.log" 2>&1
    fi
  done
  echo "[DONE] DNA quant complete."
}

############################
# STEP 5: RNA quant        #
############################
step4_rna_quant(){
  echo "============================================================"
  echo "[STEP 5] RNA quant (metatranscriptome -> TPM): $(date '+%Y-%m-%d %H:%M:%S')"

  local IDX="${SALMON_OUT}/index/CDS_index"
  [[ -d "${IDX}" ]] || die "Missing Salmon index."

  mkdir -p "${SALMON_OUT}/RNA"

  for s in "${RNA_SAMPLES[@]}"; do
    echo "[INFO] RNA sample: ${s}"
    local OUT="${SALMON_OUT}/RNA/${s}_salmon"
    if [[ -s "${OUT}/quant.sf" ]]; then
      echo "  [SKIP] quant.sf exists"
      continue
    fi
    local R="${RNA_READ_DIR}/${s}${RNA_READ_SUFFIX}"
    check_file "${R}"

    salmon quant -i "${IDX}" -l A -r "${R}" \
      -p "${THREADS}" --validateMappings -o "${OUT}" \
      > "${SALMON_OUT}/logs/${s}.RNA.salmon.log" 2>&1
  done
  echo "[DONE] RNA quant complete."
}

############################
# STEP 6: Reports          #
############################
step5_reports(){
  echo "============================================================"
  echo "[STEP 6] Reports: Geobacter contribution (%): $(date '+%Y-%m-%d %H:%M:%S')"

  local EET_BROAD="${EET_OUT}/lists/EET_genes.BROAD.ids"
  local EET_HIGH="${EET_OUT}/lists/EET_genes.HIGHCONF.ids"
  local EET_BROAD_LOC="${EET_OUT}/localization/EET_genes.BROAD.localized.ids"
  local EET_HIGH_LOC="${EET_OUT}/localization/EET_genes.HIGHCONF.localized.ids"
  check_file "${EET_BROAD}"
  check_file "${EET_HIGH}"

  # Geobacter contigs set
  local GEO_CONTIGS="${REPORT_OUT}/Geobacter_contigs.txt"
  : > "${GEO_CONTIGS}"
  for bf in "${GEO_BINS[@]}"; do
    local f="${GEO_BIN_DIR}/${bf}"
    check_file "${f}"
    grep "^>" "${f}" | sed 's/^>//' | awk '{print $1}' >> "${GEO_CONTIGS}"
  done
  sort -u "${GEO_CONTIGS}" -o "${GEO_CONTIGS}"

  # Geobacter gene lists (broad/high + localized versions)
  local GEO_EET_B="${REPORT_OUT}/Geobacter_EET_BROAD.genes"
  local GEO_EET_H="${REPORT_OUT}/Geobacter_EET_HIGHCONF.genes"
  local GEO_EET_B_LOC="${REPORT_OUT}/Geobacter_EET_BROAD.localized.genes"
  local GEO_EET_H_LOC="${REPORT_OUT}/Geobacter_EET_HIGHCONF.localized.genes"

  awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if (c in a) print $0}' \
    "${GEO_CONTIGS}" "${EET_BROAD}" | sort -u > "${GEO_EET_B}"
  awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if (c in a) print $0}' \
    "${GEO_CONTIGS}" "${EET_HIGH}" | sort -u > "${GEO_EET_H}"

  # localized lists might be empty (OK)
  if [[ -s "${EET_BROAD_LOC}" ]]; then
    awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if (c in a) print $0}' \
      "${GEO_CONTIGS}" "${EET_BROAD_LOC}" | sort -u > "${GEO_EET_B_LOC}"
  else
    : > "${GEO_EET_B_LOC}"
  fi

  if [[ -s "${EET_HIGH_LOC}" ]]; then
    awk 'NR==FNR{a[$1]=1; next} {c=$0; sub(/_[0-9]+$/,"",c); if (c in a) print $0}' \
      "${GEO_CONTIGS}" "${EET_HIGH_LOC}" | sort -u > "${GEO_EET_H_LOC}"
  else
    : > "${GEO_EET_H_LOC}"
  fi

  # sum TPM helper
  sum_tpm(){
    local gene_list="$1"
    local quant_sf="$2"
    awk 'NR==FNR{a[$1]=1; next} ($1 in a){sum+=$4} END{printf "%.6f", sum}' \
      "${gene_list}" "${quant_sf}"
  }

  # DNA + RNA report with BOTH: total EET and localized EET
  local DNA_TSV="${REPORT_OUT}/DNA_EET_contribution.tsv"
  local RNA_TSV="${REPORT_OUT}/RNA_EET_contribution.tsv"

  echo -e "sample\ttotal_EET_BROAD\tgeo_EET_BROAD\tgeo_pct_BROAD\ttotal_EET_HIGH\tgeo_EET_HIGH\tgeo_pct_HIGH\ttotal_LOC_BROAD\tgeo_LOC_BROAD\tgeo_pct_LOC_BROAD\ttotal_LOC_HIGH\tgeo_LOC_HIGH\tgeo_pct_LOC_HIGH" > "${DNA_TSV}"

  for s in "${DNA_SAMPLES[@]}"; do
    local q="${SALMON_OUT}/DNA/${s}_salmon/quant.sf"
    [[ -f "${q}" ]] || { echo "[WARN] Missing DNA quant for ${s}: ${q}"; continue; }

    total_b=$(sum_tpm "${EET_BROAD}" "${q}")
    geo_b=$(sum_tpm "${GEO_EET_B}" "${q}")
    total_h=$(sum_tpm "${EET_HIGH}" "${q}")
    geo_h=$(sum_tpm "${GEO_EET_H}" "${q}")

    total_lb="0.000000"; geo_lb="0.000000"; total_lh="0.000000"; geo_lh="0.000000"
    if [[ -s "${EET_BROAD_LOC}" ]]; then total_lb=$(sum_tpm "${EET_BROAD_LOC}" "${q}"); fi
    if [[ -s "${GEO_EET_B_LOC}" ]]; then geo_lb=$(sum_tpm "${GEO_EET_B_LOC}" "${q}"); fi
    if [[ -s "${EET_HIGH_LOC}" ]]; then total_lh=$(sum_tpm "${EET_HIGH_LOC}" "${q}"); fi
    if [[ -s "${GEO_EET_H_LOC}" ]]; then geo_lh=$(sum_tpm "${GEO_EET_H_LOC}" "${q}"); fi

    geo_pct_b=$(awk -v g="${geo_b}" -v t="${total_b}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')
    geo_pct_h=$(awk -v g="${geo_h}" -v t="${total_h}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')

    geo_pct_lb=$(awk -v g="${geo_lb}" -v t="${total_lb}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')
    geo_pct_lh=$(awk -v g="${geo_lh}" -v t="${total_lh}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')

    echo -e "${s}\t${total_b}\t${geo_b}\t${geo_pct_b}\t${total_h}\t${geo_h}\t${geo_pct_h}\t${total_lb}\t${geo_lb}\t${geo_pct_lb}\t${total_lh}\t${geo_lh}\t${geo_pct_lh}" >> "${DNA_TSV}"
  done

  echo -e "sample\ttotal_EET_BROAD\tgeo_EET_BROAD\tgeo_pct_BROAD\ttotal_EET_HIGH\tgeo_EET_HIGH\tgeo_pct_HIGH\ttotal_LOC_BROAD\tgeo_LOC_BROAD\tgeo_pct_LOC_BROAD\ttotal_LOC_HIGH\tgeo_LOC_HIGH\tgeo_pct_LOC_HIGH" > "${RNA_TSV}"

  for s in "${RNA_SAMPLES[@]}"; do
    local q="${SALMON_OUT}/RNA/${s}_salmon/quant.sf"
    [[ -f "${q}" ]] || { echo "[WARN] Missing RNA quant for ${s}: ${q}"; continue; }

    total_b=$(sum_tpm "${EET_BROAD}" "${q}")
    geo_b=$(sum_tpm "${GEO_EET_B}" "${q}")
    total_h=$(sum_tpm "${EET_HIGH}" "${q}")
    geo_h=$(sum_tpm "${GEO_EET_H}" "${q}")

    total_lb="0.000000"; geo_lb="0.000000"; total_lh="0.000000"; geo_lh="0.000000"
    if [[ -s "${EET_BROAD_LOC}" ]]; then total_lb=$(sum_tpm "${EET_BROAD_LOC}" "${q}"); fi
    if [[ -s "${GEO_EET_B_LOC}" ]]; then geo_lb=$(sum_tpm "${GEO_EET_B_LOC}" "${q}"); fi
    if [[ -s "${EET_HIGH_LOC}" ]]; then total_lh=$(sum_tpm "${EET_HIGH_LOC}" "${q}"); fi
    if [[ -s "${GEO_EET_H_LOC}" ]]; then geo_lh=$(sum_tpm "${GEO_EET_H_LOC}" "${q}"); fi

    geo_pct_b=$(awk -v g="${geo_b}" -v t="${total_b}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')
    geo_pct_h=$(awk -v g="${geo_h}" -v t="${total_h}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')

    geo_pct_lb=$(awk -v g="${geo_lb}" -v t="${total_lb}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')
    geo_pct_lh=$(awk -v g="${geo_lh}" -v t="${total_lh}" 'BEGIN{ if(t==0) print "NA"; else printf "%.3f", (g/t)*100 }')

    echo -e "${s}\t${total_b}\t${geo_b}\t${geo_pct_b}\t${total_h}\t${geo_h}\t${geo_pct_h}\t${total_lb}\t${geo_lb}\t${geo_pct_lb}\t${total_lh}\t${geo_lh}\t${geo_pct_lh}" >> "${RNA_TSV}"
  done

  echo "[DONE] Reports:"
  echo "  - ${DNA_TSV}"
  echo "  - ${RNA_TSV}"
  echo "  - loci tables: ${EET_OUT}/localization/EET_loci.BROAD.tsv , EET_loci.HIGHCONF.tsv"
}

############################
# MAIN
############################
main(){
  step1_prodigal
  step2_define_eet
  step_salmon_index
  step3_dna_quant
  step4_rna_quant
  step5_reports

  echo "============================================================"
  echo "[ALL DONE] Root output: ${ROOT_OUT}"
}

main "$@"

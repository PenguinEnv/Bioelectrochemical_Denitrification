# Pipeline steps

Run in this order. Paths are relative to the project folder (`/path/to/base` in the scripts).
DNA samples: `zr15207_9, 10, 11, 13, 14`; RNA samples: `zr15207_2T, 3T, 4T, 6T, 7T` (see README).

## 1. Metagenome

| # | Script | What it does | Main input → output |
|---|---|---|---|
| 1.1 | `1_metagenomic/1_qc/1_fastp.sh` | fastp trimming: Nextera adapter `CTGTCTCTTATACACATCT`, `--trim_poly_g --cut_tail -q 20 -l 50`; MultiQC of fastp reports | raw `*_R1/R2.fastq.gz` → `1_fastp_out/<sample>_clean_R{1,2}.fastq.gz` |
| 1.2 | `1_metagenomic/1_qc/2_fastqc.sh` | FastQC + MultiQC on trimmed reads | clean reads → `fastqc_results/` |
| 1.3 | `1_metagenomic/2_assemble/1_combine_reads.sh` | concatenate all samples' clean reads (gzip integrity check) | clean reads → `all_R1.fastq.gz`, `all_R2.fastq.gz` |
| 1.4 | `1_metagenomic/2_assemble/2_megahit.sh` | MEGAHIT co-assembly, `--min-contig-len 1000` | combined reads → `megahit_out/final.contigs.fa` |
| 1.5 | `1_metagenomic/3_mapping/bowtie2.sh` | Bowtie2 index of contigs; `--very-sensitive` mapping per sample; sort, index, flagstat | contigs + clean reads → `<sample>.sorted.bam` |
| 1.6 | `1_metagenomic/4_coverage/1_coverm_contig.sh` | CoverM contig TPM and mean coverage (`--min-read-percent-identity 0.95 --min-read-aligned-percent 0.75`) | BAMs → `coverm_contig_tpm.tsv`, `coverm_contig_mean.tsv` |
| 1.7 | `1_metagenomic/5_binning/1_vamb/vamb.sh` | VAMB `bin default`, `--minfasta 200000` | contigs + contig mean coverage → `vamb_out/bins/` |
| 1.8 | `1_metagenomic/5_binning/2_concoct/1_concoct.sh` | cut contigs into 10 kb pieces, coverage table, CONCOCT clustering (`-c 200`) | contigs + BAMs → `concoct_out/clustering_gt1000.csv` |
| 1.9 | `1_metagenomic/5_binning/2_concoct/2_concoct_extract_200k.sh` | merge cut-up clustering, extract bin FASTAs, keep bins >= 200 kb | → `concoct_out/clustering_merged.csv`, `bins_200k/` |
| 1.10 | `1_metagenomic/5_binning/3_maxbin2/1_maxbin2.sh` | per-sample abundance files from CoverM mean, MaxBin2, keep bins >= 200 kb | contigs + contig mean coverage → `maxbin_out/` |
| 1.11 | `1_metagenomic/5_binning/3_maxbin2/2_maxbin2_filter_200k.sh` | re-run only the >= 200 kb filter on MaxBin2 bins | → `maxbin_out/bins_200k/` |
| 1.12 | `1_metagenomic/5_binning/4_semibin/semibin.sh` | SemiBin2: features from 5 BAMs, self-supervised training, `bin_short`; keep bins >= 200 kb | contigs + BAMs → `4_semibin/3_bin/` |
| 1.13 | `1_metagenomic/6_bincheck/1_tsv/contig2bin_tsv.sh` | contig-to-bin tables for each binner (`vamb_`, `semibin_`, `concoct_bin_` prefixes) | bins → `vamb.tsv`, `maxbin.tsv`, `semibin.tsv`, `concoct_fixed.tsv` |
| 1.14 | `1_metagenomic/6_bincheck/2_dastool/dastool.sh` | DAS Tool integration of the 4 binners (`--write_bins --write_bin_evals`) | tables + contigs → `DAS_Tool_DASTool_bins/` |
| 1.15 | `1_metagenomic/6_bincheck/3_drep/drep.sh` | dRep: `--completeness 50 --contamination 10 --S_algorithm ANImf --S_ani 0.95 --cov_thresh 0.25` (CheckM `taxonomy_wf`) | DAS Tool bins → `dereplicated_genomes/` (final MAGs) |
| 1.16 | `1_metagenomic/4_coverage/2_coverm_genome.sh` | CoverM genome TPM and relative abundance of MAGs | BAMs + MAGs → `coverm_genome.tsv` |
| 1.17 | `1_metagenomic/4_coverage/3_coverm_genome_mean.sh` | CoverM genome mean coverage of MAGs | BAMs + MAGs → `coverm_genome_mean.tsv` |
| 1.18 | `1_metagenomic/7_interpretation/1_gtdbtk/1_gtdbtk.sh` | GTDB-Tk `classify_wf` on DAS Tool bins and on dRep MAGs | MAGs → `gtdbtk_out_drep/` (summary, alignment, markers) |
| 1.19 | `.../1_gtdbtk/FastTree/1_fasttree.sh` | FastTree (`-lg -gamma`) on the GTDB-Tk bac120 alignment | MSA → `MAG_bac120_user.nwk` |
| 1.20 | `.../1_gtdbtk/FastTree/2_itol_annotation.sh` | iTOL annotation files (family/genus color strips, species labels) | GTDB-Tk summary → `itol_*.txt` |
| 1.21 | `1_metagenomic/7_interpretation/2_enrichm/split_contig_enrichm.sh` | split contigs into chunks of 10,000 contigs; EnrichM `annotate --ko` (DIAMOND) | contigs → `annot_output/annotations_ko/DIAMOND_search.tsv`, GFFs |
| 1.22 | `1_metagenomic/7_interpretation/3_Contigs_TPM/contig_function_tpm.sh` | best DIAMOND hit per gene; contig TPM split among its genes by alignment length; summed per KO; KEGG descriptions and pathways | contig TPM + DIAMOND → `KO_TPM.clean.tsv` |
| 1.23 | `1_metagenomic/7_interpretation/4_MAG_TPM/function_MAG_tpm.sh` | contig-to-MAG table from dRep MAGs; MAG abundance; MAG × KO and KO TPM (length-weighted) | → `contig_to_bin.tsv`, `MAG_abundance_by_sample.tsv`, `MAG_KO_TPM.tsv`, `KO_TPM.tsv` |
| 1.24 | `1_metagenomic/7_interpretation/5_KO_MAX_TPM/KO_MAX_TPM.sh` | top 3 MAGs per KO in each DNA sample | `MAG_KO_TPM.tsv` → `<sample>_KO_Top3MAG.tsv` |

## 2. Metatranscriptome

| # | Script | What it does | Main input → output |
|---|---|---|---|
| 2.1 | `2_metatranscriptomic/1_qc/1_fastp/1_fastp.sh` | fastp `--cut_tail --trim_poly_g -q 20 -u 30 -n 5 -l 30`; MultiQC | raw RNA reads → `fastp_out/<sample>_clean_R{1,2}.fastq.gz` |
| 2.2 | `2_metatranscriptomic/1_qc/1_fastp/2_fastqc.sh` | FastQC + MultiQC on trimmed reads | → `fastqc_results/` |
| 2.3 | `2_metatranscriptomic/1_qc/2_sortmerna/sortmerna.sh` | SortMeRNA against SILVA 138.1 LSU + SSU NR99 (`--paired_in --fastx`), keep non-rRNA reads | clean reads → `<sample>_mRNA_SILVA.fq.gz` |
| 2.4 | `2_metatranscriptomic/2_mapping/salmon.sh` | Salmon index of metagenome contigs (`-k 31`); `quant -l A --validateMappings` | mRNA reads + contigs (1.4) → `<sample>_salmon_quant/quant.sf` |
| 2.5 | `2_metatranscriptomic/3_TPM/mag_contig_TPM.sh` | contig TPM split among CDS by CDS length (EnrichM GFF); aggregated to KO and MAG × KO; KEGG annotation | quant.sf + 1.21 + 1.23 → `KO_TPM_contig_level_with_annotation.csv`, `MAG_KO_TPM_with_annotation.csv` |
| 2.6 | `2_metatranscriptomic/4_KO_MAX_TPM/KO_MAX_TPM.sh` | top 3 MAGs per KO in each RNA sample | 2.5 → `<sample>_KO_Top3MAG.tsv` |

## 1-2. Transcription normalized by gene copy (DNA + RNA)

| # | Script | What it does | Main input → output |
|---|---|---|---|
| N.1 | `1_metagenomic/7_interpretation/6_Gene_normal/1_normalization.sh` | community genome equivalents from GTDB-Tk bac120 single-copy markers; KO DNA coverage per genome equivalent; KO copies per MAG; RNA TPM / DNA ratios at community and MAG level (DNA-RNA pairs in the built-in sample map) | GTDB-Tk markers, DIAMOND, contig/MAG coverage, 2.5 → `community_KO_activity_RNA_over_DNA_GTDB_marker.tsv`, `MAG_KO_activity_RNA_over_DNAcopy.tsv` (+ filtered versions) |
| N.2 | `1_metagenomic/7_interpretation/6_Gene_normal/2_group_weighted.sh` | sum RNA and DNA per MAG × function group (nitrate/nitrite reduction, denitrification, DNRA, N fixation, hydrogenases, oxidative stress, terminal oxidases; KO list inside the script); weighted activity = ΣRNA / ΣDNA copy | N.1 → `grouped_MAG_activity/` |
| N.3 | `1_metagenomic/7_interpretation/6_Gene_normal/3_pathway.sh` | same summary per KEGG pathway (broad pathways optionally excluded) | N.1 → `pathway_MAG_activity/` |

## 3. EET genes

| # | Script | What it does | Main input → output |
|---|---|---|---|
| 3.1 | `3_Geobacter_EET/1_blastdb/geobacter_db.sh` | download proteomes of *G. sulfurreducens* PCA, *G. metallireducens* GS-15, *G. anodireducens* GSS01 from NCBI; DIAMOND database | → `Geobacter_major_ref_proteins.faa`, `.dmnd` |
| 3.2 | `3_Geobacter_EET/2_hmmer_diamond/1_hmmer_diamond.sh` | for the three Geobacter MAGs: Prodigal, Pfam hmmscan, DIAMOND blastp (`-e 1e-10 --query-cover 50 --subject-cover 50`) | Geobacter MAGs → `EET_check/` |
| 3.3 | `3_Geobacter_EET/2_hmmer_diamond/2_EET_organize.sh` | per MAG: Pfam EET candidates, DIAMOND hits to EET-like Geobacter proteins, multiheme cytochromes (>= 3 CXXCH) | 3.2 → `EET_organize/reports/` |
| 3.4 | `3_Geobacter_EET/3_community_EET/1_pfam_eet_subset.sh` | subset Pfam-A to EET-related families (cytochrome, multiheme, Mtr/Omc, Ext, type IV pili keywords) | Pfam-A.hmm → `Pfam_eet/EET_Pfam_subset.hmm` |
| 3.5 | `3_Geobacter_EET/3_community_EET/2_community_eet.sh` | Prodigal `-p meta` on all contigs; Pfam subset `hmmscan --cut_ga`; CXXCH >= 3; DIAMOND vs Geobacter EET-like proteins. Gene sets: BROAD (union) and HIGHCONF (Pfam ∩ DIAMOND, plus multiheme); localized loci (>= 2 EET genes within 10 kb). Salmon on the CDS catalog with DNA and RNA reads; total vs Geobacter EET | contigs, Geobacter MAGs, clean DNA reads, mRNA reads → `1_prodigal/` … `4_reports/` |
| 3.6 | `3_Geobacter_EET/3_community_EET/3_EET_MAGs.sh` | EET gene abundance and transcription contributed by selected MAGs (list in the script) | 3.5 → `4_reports/MAG_contributions/` |

## 4. Correlation

| # | Script | What it does | Main input → output |
|---|---|---|---|
| 4.1 | `4_correlation/1_correlation.sh` | Pearson r and p-values between all variables | `data.csv` → `correlation_r.csv`, `correlation_p.csv`, `correlation_matrix.csv` |
| 4.2 | `4_correlation/2_figure.sh` | lower-triangle heatmap with significance stars (* p < 0.05, ** p < 0.01) | `data.csv` → `correlation_heatmap_WR_style.png/.pdf` |

`data.csv` is not included. `data_template.csv` gives the column layout used in the paper:
reactor performance (`NO3`, `e/N`) and nitrate, nitrite, H2, ROS, O2 and EET function values in four groups:

| Prefix | Meaning |
|---|---|
| `CA` | community-wide abundance (DNA) |
| `CE` | community-wide expression (RNA) |
| `DA` | dominant taxa abundance (DNA) |
| `DE` | dominant taxa expression (RNA) |

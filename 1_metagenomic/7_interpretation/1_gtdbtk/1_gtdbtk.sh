#!/bin/bash
# Requires: GTDB-Tk 2.4.1 (GTDBTK_DATA_PATH must point to the GTDB reference data)
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.
THREADS=16

# Paths
INPUT_das="/path/to/base/1_metagenomic/6_bincheck/2_dastool/DAS_Tool_DASTool_bins"
OUTDIR_das="/path/to/base/1_metagenomic/7_interpretation/1_gtdbtk/gtdbtk_out_das"

INPUT_drep="/path/to/base/1_metagenomic/6_bincheck/3_drep/dereplicated_genomes"
OUTDIR_drep="/path/to/base/1_metagenomic/7_interpretation/1_gtdbtk/gtdbtk_out_drep"

mkdir -p "$OUTDIR_das"
mkdir -p "$OUTDIR_drep"

# Run GTDB-Tk classify workflow
gtdbtk classify_wf \
  --genome_dir "$INPUT_das" \
  --out_dir "$OUTDIR_das" \
  --extension fa \
  --mash_db "$OUTDIR_das/gtdb_ref_sketch.msh" \
  --cpus 8

# Run GTDB-Tk classify workflow
gtdbtk classify_wf \
  --genome_dir "$INPUT_drep" \
  --out_dir "$OUTDIR_drep" \
  --extension fa \
  --mash_db "$OUTDIR_drep/gtdb_ref_sketch.msh" \
  --cpus 8

#!/bin/bash
#
# Requires: FastTree 2.1.10
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.
export OMP_NUM_THREADS=2

zcat /path/to/base/1_metagenomic/7_interpretation/1_gtdbtk/gtdbtk_out_drep/align/gtdbtk.bac120.user_msa.fasta.gz \
  | FastTree -lg -gamma \
  > MAG_bac120_user.nwk

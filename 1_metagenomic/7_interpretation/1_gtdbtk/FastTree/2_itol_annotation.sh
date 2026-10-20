#!/bin/bash
# Requires: bash, awk
# Paths below are placeholders (/path/to/base = your project folder, /path/to/your_db = reference databases); edit before running.
set -euo pipefail

IN="/path/to/base/1_metagenomic/7_interpretation/1_gtdbtk/gtdbtk_out_drep/gtdbtk.bac120.summary.tsv"

awk -F'\t' '
function hashcolor(s,   i, c, hex) {
  # label -> simple hash so the same label always gets the same color
  c = 0
  for (i=1; i<=length(s); i++) c = (c*33 + ord(substr(s,i,1))) % 16777215
  hex = sprintf("#%06X", c)
  return hex
}
function ord(ch){
  return index("\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1A\x1B\x1C\x1D\x1E\x1F !\"#$%&'\''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", ch)
}

BEGIN {
  # Family color strip
  print "DATASET_COLORSTRIP" > "itol_family.txt"
  print "SEPARATOR\tTAB" >> "itol_family.txt"
  print "DATASET_LABEL\tFamily" >> "itol_family.txt"
  print "STRIP_WIDTH\t18" >> "itol_family.txt"
  print "SHOW_INTERNAL\t0" >> "itol_family.txt"
  print "DATA" >> "itol_family.txt"

  # Genus color strip
  print "DATASET_COLORSTRIP" > "itol_genus.txt"
  print "SEPARATOR\tTAB" >> "itol_genus.txt"
  print "DATASET_LABEL\tGenus" >> "itol_genus.txt"
  print "STRIP_WIDTH\t18" >> "itol_genus.txt"
  print "SHOW_INTERNAL\t0" >> "itol_genus.txt"
  print "DATA" >> "itol_genus.txt"

  # Species text labels
  print "DATASET_TEXT" > "itol_species.txt"
  print "SEPARATOR\tTAB" >> "itol_species.txt"
  print "DATASET_LABEL\tSpecies" >> "itol_species.txt"
  print "COLOR\t#000000" >> "itol_species.txt"
  print "SIZE\t2" >> "itol_species.txt"
  print "DATA" >> "itol_species.txt"
}

NR>1 {
  # GTDB-Tk summary.tsv: col1=user_genome, col2=classification (usual layout)
  split($2, a, ";")
  fam="Unclassified"; gen="Unclassified"; sp=""

  for (i in a) {
    if (a[i] ~ /^f__/) {fam=a[i]; sub(/^f__/, "", fam)}
    if (a[i] ~ /^g__/) {gen=a[i]; sub(/^g__/, "", gen)}
    if (a[i] ~ /^s__/) {sp=a[i]; sub(/^s__/, "", sp)}
  }

  # clean empty/odd values
  if (fam == "" || fam == "uncultured") fam="Unclassified"
  if (gen == "" || gen == "uncultured") gen="Unclassified"

  # if species is missing, use genus + " sp."
  if (sp == "" || sp == "uncultured") {
    if (gen != "Unclassified") sp = gen " sp."
    else sp = "Unclassified sp."
  } else {
    # GTDB placeholder species (e.g. sp905339345) are kept as they are;
    # a name like "Geobacter sp905339345" is safest left unchanged
  }

  # Colorstrip: ID  COLOR  LABEL
  print $1 "\t" hashcolor(fam) "\t" fam >> "itol_family.txt"
  print $1 "\t" hashcolor(gen) "\t" gen >> "itol_genus.txt"

  # Text: ID  TEXT
  print $1 "\t" sp >> "itol_species.txt"
}

' "$IN"

echo "✅ Done:"
echo " - itol_family.txt (Color strip)"
echo " - itol_genus.txt  (Color strip)"
echo " - itol_species.txt (Text labels; filled with genus sp. if missing)"

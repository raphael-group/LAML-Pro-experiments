#!/bin/bash

BASE="/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/BaseMEM_Magic/Magic/Figures3_5_and_SupFigures2_4-6/lineage_analysis/BEAST_analysis/Trial4"
OUTBASE="/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/inputs/trees/"

files=(
  "$BASE/bigmem_test_all_trees_plus_space_geo-p2_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p3_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p4_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p5_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p6_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p78_barcs.trees.tree"
  "$BASE/bigmem_test_all_trees_plus_space_geo-p910_barcs.trees.tree"
)

for fname in "${files[@]}"; do
  echo "Processing: $fname"
  basefname="${fname##*/}"
  stem="${basefname%.trees.tree}"

  python3 nexus_to_newick.py "$fname" > "${OUTBASE}/${stem}.both.newick"

  infile="${OUTBASE}/${stem}.both.newick"
  full="${OUTBASE}/${stem}.full.newick"
  short="${OUTBASE}/${stem}.newick"

  sed -n '1p' "$infile" > "$full"
  sed -n '2p' "$infile" > "$short"

  out=$(nw_stats $short)
  echo $out

done


#!/bin/bash

outfile="tmp"

round1_chars="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/PEtracer-2025/barcoded_tracing/data/round1_annotations.txt"
round2_chars="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/PEtracer-2025/barcoded_tracing/data/round2_annotations.txt"

laml_dir="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs"
petracer_tree="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/barcode/clone4_exp3.published_petracer.binary.nwk"

laml_trees=(
  $laml_dir/barcode_experiment3/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.01_tree.newick
  $laml_dir/barcode_experiment3/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.03_tree.newick
  $laml_dir/barcode_experiment3/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.05_tree.newick
  $laml_dir/barcode_experiment3/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.07_tree.newick
  $laml_dir/barcode_experiment3/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.10_tree.newick
  # plus the pt2 topo tree 
  $laml_dir/barcode_experiment3_pt2/fastlaml_barcode_fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.10_tree_topo_binary_pt2_tree.newick
)

echo "Round 1"
char_file="$round1_chars"

echo "  LAML-Pro"
for tree_file in "${laml_trees[@]}"; do
  echo "    $(basename "$tree_file")"
  python helpers/small_MP.py "${char_file}" "${tree_file}" "${outfile}"
done

echo "  PEtracer"
python helpers/small_MP.py "${char_file}" "${petracer_tree}" "${outfile}"

echo
echo "Round 2"
char_file="$round2_chars"

echo "  LAML-Pro"
for tree_file in "${laml_trees[@]}"; do
  echo "    $(basename "$tree_file")"
  python helpers/small_MP.py "${char_file}" "${tree_file}" "${outfile}"
done

echo "  PEtracer"
python helpers/small_MP.py "${char_file}" "${petracer_tree}" "${outfile}"


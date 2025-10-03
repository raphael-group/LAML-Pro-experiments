#!/bin/bash

inputs_basename="/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/inputs/"

colony=2
lamlpro_basename="/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/runjobs/outputs_baseMemoir_colony${colony}/fastlaml_baseMemoir.colony2.published"

# run spatial analysis
bm_treefile="${inputs_basename}/trees/baseMemoir.colony${colony}.published.scaled.newick"
lp_treefile="${lamlpro_basename}_tree.scaled.newick"
centroids_fname="${inputs_basename}/colony${colony}_centroids.txt"
#(ete3-py312) gc3045@PU-T0264KMR6R scripts % python euclidean_solver.py
echo "LAML-Pro:"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $lp_treefile --leaf_df ${centroids_fname} --output "colony${colony}_lp_ancestral_labeling.txt" --use_branch_length
echo "baseMemoir:"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $bm_treefile --leaf_df ${centroids_fname} --output "colony${colony}_bm__ancestral_labeling.txt" --use_branch_length

########################################################
colony=5
lamlpro_basename="/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/runjobs/outputs_baseMemoir_colony${colony}/fastlaml_baseMemoir.colony5.stepwise_addition_03"

# run spatial analysis
bm_treefile="${inputs_basename}/trees/baseMemoir.colony${colony}.published.scaled.newick"
lp_treefile="${lamlpro_basename}_tree.scaled.newick"
centroids_fname="${inputs_basename}/colony${colony}_centroids.txt"

#(ete3-py312) gc3045@PU-T0264KMR6R scripts % python euclidean_solver.py
echo "LAML-Pro:"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $lp_treefile --leaf_df ${centroids_fname} --output "colony${colony}_lp_ancestral_labeling.txt" --use_branch_length
echo "baseMemoir:"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $bm_treefile --leaf_df ${centroids_fname} --output "colony${colony}_bm__ancestral_labeling.txt" --use_branch_length


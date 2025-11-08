#!/bin/bash


colony=3
centroids_fname="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/colonies_clone3_centroids_subset.txt"

# run spatial analysis
pet_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/colonies_clone3.petracer.neighbor_joining.scaled.newick"
#pet_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/colonies_clone3.petracer.neighbor_joining.nwk"
#pet_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/colonies/clone3.petracer.neighbor_joining.scaled.newick"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1/fastlaml_colonies_fastlaml_colonies_clone3.petracer.neighbor_joining_tree.round2_tree.scaled.newick"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1/fastlaml_colonies_fastlaml_colonies_clone3_cmat.petracer.neighbor_joining_tree.round3_tree_tree.scaled.newick"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1/fastlaml_colonies_clone3.petracer.neighbor_joining.binary_tree.scaled.newick"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1/fastlaml_colonies_fastlaml_colonies_clone3_cmat.petracer.neighbor_joining_tree.round3_tree_run2_tree.scaled.newick"
lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1_pt2/fastlaml_colonies_fastlaml_colonies_clone3.petracer.neighbor_joining.binary_run2_tree_topo_binary_pt2_tree.scaled.newick"
laml_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1_cmat/fastlaml_colonies_cmat_fastlaml_colonies_clone3.petracer.neighbor_joining_tree.round2_tree.scaled.newick"
# Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1_cmat/fastlaml_colonies_cmat_clone3.petracer.neighbor_joining.binary_tree.scaled.newick"

echo "Experiment 1"
echo "LAML-Pro: ${lp_treefile}"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $lp_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp1_lp_ancestral_labeling.txt" --use_branch_length #--epsilon 1346 # 0.06 * scale of tree #0.02 #--epsilon 2500 #0.02

echo "PETracer: ${pet_treefile}"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $pet_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp1_pet_ancestral_labeling.txt"  --use_branch_length #--epsilon 1346 #0.02 #--epsilon 2500

echo "LAML: ${laml_treefile}"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $laml_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp1_laml_ancestral_labeling.txt"  --use_branch_length #--epsilon 1346 #0.02 #--epsilon 2500
##################################################################################################
echo "Experiment 2"
lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment2/fastlaml_colonies_clone3_exp2.neighbor_joining_tree.scaled.newick"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment2/outputs_petracer_colonies_clone3/fastlaml_colonies_clone3.neighbor_joining_tree.subset.scaled.newick"

centroids_fname="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/colonies_clone3_centroids_filtered.txt"
#echo "LAML-Pro induced on 695 cells: ${lp_treefile}"
python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $lp_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp2_lp_ancestral_labeling.txt" --use_branch_length #--epsilon 1346

#echo "PETracer: ${pet_treefile}"
#python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $pet_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp1_pet_ancestral_labeling.txt"  --use_branch_length

##################################################################################################

#centroids_fname="${inputs_basename}/colonies_clone${colony}_centroids.txt"
#lp_treefile="/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment2/outputs_petracer_colonies_clone3/fastlaml_colonies_clone3.neighbor_joining_tree.scaled.newick"

#echo "LAML-Pro on 832 cells: ${lp_treefile}"
#python /Users/gc3045/git/fast-laml/scripts/euclidean_solver.py --tree $lp_treefile --leaf_df ${centroids_fname} --output "colonies_clone${colony}_exp2_lp_ancestral_labeling.txt"  --use_branch_length


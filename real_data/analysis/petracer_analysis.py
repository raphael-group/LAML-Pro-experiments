import pandas as pd
import numpy as np
import subprocess
from helpers.utils import build_summary_df, plot_genotype_confidence, clustermap_genos, plot_genotypecall_summary

from helpers.utils import distdict_to_df, leaf_pairs, get_geno_dict
from helpers.utils import im_ehd, empirical_site_dists, ehd, sm_ehd, pair_metrics
from helpers.utils import plot_concordance_scatterplot, plot_concordance_distribution, plot_state_counts
from helpers.utils import report_genotype_call_stats, save_df_to_pdf, plot_correlation, branch_table, plot_bl_variance, add_internal_labels
from helpers.utils import plot_tree_3d, edge_ratio_table

import treeswift
import sys
import matplotlib.pyplot as plt

plt.ioff()

lamlpro_basename = sys.argv[1]
colony = sys.argv[2]

####################################################################
plotdir = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/plots/"
inputs_basename = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/"
analysis_basename = "/Users/gc3045/git/laml2-experiments/real_data/analysis/"

#bm_input_geno = inputs_basename + f"baseMemoir_colony{colony}_baseMemoir_genotypes.csv"
#lp_map_geno = lamlpro_basename + "_posterior_argmax.csv"
pe_tree_fname = inputs_basename + f"trees/colonies/clone{colony}.petracer.neighbor_joining.binary.nwk"
lp_tree_fname = lamlpro_basename + "_tree.newick"

#bm_input_geno_df = pd.read_csv(bm_input_geno)
#lp_map_geno_df = pd.read_csv(lp_map_geno, skiprows=2, index_col=0)

#summary_df = build_summary_df(bm_input_geno_df, lp_map_geno_df)

bm_tree = treeswift.read_tree_newick(bm_tree_fname)
lp_tree = treeswift.read_tree_newick(lp_tree_fname)

#centroids_fname = inputs_basename + f"colony{colony}_centroids.txt"
#centroids_df = pd.read_csv(centroids_fname, sep="\t", header=None, index_col=0)

tau = 72 # hours
short_branch_thresh = 0.0005

####### Metrics
N, K = lp_map_geno_df.loc[[n for n in lp_map_geno_df.index if not n.startswith('internal')],].shape
baseMemoir_impute_per = np.sum(bm_input_geno_df['bM_pmax'] > 0.70) / (N * 396)
print("BaseMemoir Impute %:", baseMemoir_impute_per, "over", N, "cells")

# robinson-foulds
output = subprocess.run([
    "python", "/Users/gc3045/git/fast-laml/scripts/compare_two_trees.py",
    "-t1", bm_tree_fname,
    "-t2", lp_tree_fname
])
print("Robinson-Foulds distance between published baseMemoir tree and LAML-Pro tree:", output)

####### Genotype plausibility 
"""
plot_genotype_confidence(summary_df, title=f"Colony {colony}: baseMemoir probabilities\nLAML-Pro vs baseMemoir genotypes",
                         outfile=f"{plotdir}/colony{colony}_genotype_probabilities.pdf",
                         bins=100)

bm_geno_df = bm_input_geno_df[['cell_name', 'target_site', 'bM_geno']]
bm_geno_df = bm_geno_df.pivot(
    index="cell_name",
    columns="target_site",
    values="bM_geno"
).sort_index(axis=1).fillna(-1)

lp_map_geno_df.columns = bm_geno_df.columns
lp_map_geno_df = lp_map_geno_df.loc[[x for x in lp_map_geno_df.index if not x.startswith('internal')], :]
lp_map_geno_df.index = [int(x) for x in lp_map_geno_df.index]
bm_geno_df.index = [int(x) for x in bm_geno_df.index]

clustermap_genos(bm_geno_df, lp_map_geno_df, metric='hamming', method='complete', title=f"Colony {colony}", outfile=f"{plotdir}/colony{colony}_genotype_comparison.pdf")

counts = plot_state_counts(bm_geno_df, lp_map_geno_df, title=f"Colony {colony}: BaseMemoir vs. LAML-Pro",
                           outfile=f"{plotdir}/colony{colony}_genotype_counts.pdf")

counts_trim = (
    counts[counts.index != '-1']
    .drop(columns='-1'))
unedited = ['AA']            # baseMemoir (rows): AA = unedited
edited   = ['GG','GA','AG']  # GG/GA/AG = edited

uu = counts_trim.loc[unedited, unedited].to_numpy().sum()  # unedited → unedited
ue = counts_trim.loc[unedited, edited].to_numpy().sum()    # unedited → edited
eu = counts_trim.loc[edited, unedited].to_numpy().sum()    # edited → unedited
ee = counts_trim.loc[edited, edited].to_numpy().sum()      # edited → edited

summary_2x2 = pd.DataFrame(
    [[uu, ue],
     [eu, ee]],
    index=pd.Index(['Unedited','Edited'], name=''),
    columns=pd.Index(['Unedited','Edited'])
)

shared_sites = [c for c in bm_geno_df.columns if c in lp_map_geno_df.columns]

# Flatten both matrices into tidy long format
bm_long = bm_geno_df[shared_sites].stack().reset_index()
bm_long.columns = ["cell_name", "target_site", "bm_state"]

lp_long = lp_map_geno_df[shared_sites].stack().reset_index()
lp_long.columns = ["cell_name", "target_site", "lp_state"]

merged = (
    bm_input_geno_df.merge(bm_long, on=["cell_name", "target_site"])
                    .merge(lp_long, on=["cell_name", "target_site"])
)
fig, ax = plot_genotypecall_summary(
    counts_trim, merged, y_max=1.0, outfile=f"{plotdir}/colony{colony}_genotypecall_summary.pdf"
)

stats_df = report_genotype_call_stats(counts)
save_df_to_pdf(stats_df, filename=f"{plotdir}/colony{colony}_genotype_call_stats.pdf", title=None, floatfmt="{:.4f}", fontsize=10)

####### Genotype concordance with tree topology

bm_distmat = distdict_to_df(bm_tree.distance_matrix(leaf_labels=True))
lp_distmat = distdict_to_df(lp_tree.distance_matrix(leaf_labels=True))

bm_leaf_pairs = leaf_pairs(bm_distmat)
lp_leaf_pairs = leaf_pairs(lp_distmat)

bm_genodict = get_geno_dict(bm_input_geno_df)
hat = empirical_site_dists(bm_genodict)

bm_concordance = pair_metrics(bm_leaf_pairs, bm_genodict, hat)
lp_concordance = pair_metrics(lp_leaf_pairs, bm_genodict, hat)

plot_concordance_scatterplot(bm_concordance, lp_concordance,
                          outfile=f"{plotdir}/colony{colony}_genotype_scatterplot.pdf")

plot_concordance_distribution(bm_concordance, lp_concordance, bm_genodict, hat,
                          outfile=f"{plotdir}/colony{colony}_genotype_distribution.pdf")

####### Spatial annotation

# scale trees first
lp_tree.collapse_short_branches(short_branch_thresh)
lp_tree.root.set_edge_length(0.0)

bm_tree.scale_edges(tau/bm_tree.height())
lp_tree.scale_edges(tau/lp_tree.height())

var_x, var_y = centroids_df.var(axis=0)
print("Variance", var_x, var_y)
diffusion_scale = (var_x + var_y)/(2*tau)

bm_tree.scale_edges(diffusion_scale)
lp_tree.scale_edges(diffusion_scale)

add_internal_labels(bm_tree)
add_internal_labels(lp_tree)

print("Tree heights", bm_tree.height(), lp_tree.height())

bm_tree.write_tree_newick(bm_tree_fname[:-4] + ".scaled.newick")
lp_tree.write_tree_newick(lp_tree_fname[:-7] + ".scaled.newick")

# compute correlation

lp_tree = treeswift.read_tree_newick(lp_tree_fname[:-7] + ".scaled.newick")
bm_tree = treeswift.read_tree_newick(bm_tree_fname[:-4] + ".scaled.newick")

bm_distmat = distdict_to_df(bm_tree.distance_matrix(leaf_labels=True))
lp_distmat = distdict_to_df(lp_tree.distance_matrix(leaf_labels=True))

# plot phylogenetic distance against spatial distance and compute pearson correlation
fig, axes = plot_correlation(centroids_df, bm_distmat, lp_distmat, title=f"Colony {colony} Leaf-pair distances", outfile=f"{plotdir}/colony{colony}_spatial_correlation.pdf")

# plot branch length variance
bm_tree_df = branch_table(bm_tree)
lp_tree_df = branch_table(lp_tree)

plot_bl_variance(bm_tree_df, lp_tree_df, title=None, 
                 outfile=f"{plotdir}/colony{colony}_branchlen_variance.pdf")

# generate the 3D interactive plot

lp_tree = treeswift.read_tree_newick(lp_tree_fname)
bm_tree = treeswift.read_tree_newick(bm_tree_fname)

lp_tree.collapse_short_branches(short_branch_thresh)
lp_tree.root.set_edge_length(0.0)

bm_tree.scale_edges(tau/bm_tree.height())
lp_tree.scale_edges(tau/lp_tree.height())

centroids_df.columns = ["X", "Y"]
ancestral_labeling_df = pd.read_csv(f"{analysis_basename}/colony{colony}_bm_ancestral_labeling.txt", index_col=0)
ancestral_labeling_df.columns = centroids_df.columns
bm_merged_df = pd.concat([ancestral_labeling_df, centroids_df], axis=0)

centroids_df.columns = ["X", "Y"]
ancestral_labeling_df = pd.read_csv(f"{analysis_basename}/colony{colony}_lp_ancestral_labeling.txt", index_col=0)
ancestral_labeling_df.columns = centroids_df.columns
lp_merged_df = pd.concat([ancestral_labeling_df, centroids_df], axis=0)

fig = plot_tree_3d(lp_tree, lp_merged_df, x_col="X", y_col="Y",
                   title="LAML-Pro: Cells in 3D",
                   outfile=f"{plotdir}/colony{colony}_lamlpro_tree.html")

fig = plot_tree_3d(bm_tree, bm_merged_df, x_col="X", y_col="Y",
                   title="baseMemoir: Cells in 3D",
                   outfile=f"{plotdir}/colony{colony}_baseMemoir_tree.html")
"""

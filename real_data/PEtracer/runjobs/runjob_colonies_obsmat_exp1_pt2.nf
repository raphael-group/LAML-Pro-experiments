clone=3
params.matrix = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/colonies_clone${clone}_kde_scores_subset.csv"
params.matrix_type = "observation-matrix"
params.fastlaml_bin = "/Users/gc3045/git/fast-laml/build/src/fastlaml"
params.seed = 0
params.max_iter = 10000 
params.tree_dir = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/colonies_experiment1"
params.nw_topology = "nw_topology"   
params.python = "/Users/gc3045/miniconda3/bin/python3" 

Channel
    .fromPath("${params.tree_dir}/*run2*.newick")
    .set { tree_files }

// ---- topology extraction ----
process make_topology {
    tag "${tree.getBaseName()}"

    input:
    path tree

    output:
    path "${tree.getBaseName()}_topo_binary.newick", emit: topo

    publishDir "colonies_experiment1_pt2", mode: 'copy'

    script:
    def basename = tree.getBaseName()
    """
    echo "Extracting topology from: ${tree}"
    ${params.nw_topology} -I ${tree} > ${basename}_topo.newick

    echo "Resolving polytomies with TreeSwift to binarize tree"
    ${params.python} - <<'PY'
    import sys
    import treeswift as ts

    inp = "${basename}_topo.newick"
    outp = "${basename}_topo_binary.newick"

    with open(inp, "r") as f:
        newick_str = f.read()
    t = ts.read_tree_newick(newick_str)

    t.resolve_polytomies()

    with open(outp, "w") as f:
        f.write(t.newick())

    print(outp)
    PY

    echo "Wrote: ${basename}_topo_binary.newick"
    """
}

process run_fastlaml {
    tag "${tree.getBaseName()}"

    input:
    path tree

    output:
    path "fastlaml_*", emit: results
    publishDir "colonies_experiment1_pt2", mode: 'copy'

    script:
    def basename = tree.getBaseName()
    def outdir = "fastlaml_colonies_${basename}_pt2"

    """
    echo "Running fastlaml on tree topology: ${tree}"
    ${params.fastlaml_bin} \\
        --matrix ${params.matrix} \\
        --tree ${tree} \\
        --output ${outdir} \\
        --mode search \\
        --seed ${params.seed} \\
        --max-iterations ${params.max_iter} \\
        --ultrametric \\
        --min-branch-length 0.01 \\
        -d ${params.matrix_type} \\
        -v
    echo "Finished: ${tree} -> ${outdir}"
    """
}

workflow {
    make_topology(tree_files)
    run_fastlaml(make_topology.out.topo)
}

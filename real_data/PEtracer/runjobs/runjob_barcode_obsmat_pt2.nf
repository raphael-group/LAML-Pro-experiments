params.matrix = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/barcode_tracing/state_probs_eta_0.10.csv"
params.matrix_type = "observation-matrix"
params.fastlaml_bin = "/Users/gc3045/git/fast-laml/build/src/fastlaml"
params.seed = 0
params.max_iter = 20000 
params.tree_dir = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/runjobs/barcode_experiment3"
params.nw_topology = "nw_topology"   
params.python = "/Users/gc3045/miniconda3/bin/python3" 

Channel
    .fromPath("${params.tree_dir}/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.01_tree.newick")
    //.fromPath("${params.tree_dir}/fastlaml_barcode_clone4_exp3.published_petracer.binary_state_probs_eta_0.10_tree.newick")
    .set { tree_files }

// ---- topology extraction ----
process make_topology {
    tag "${tree.getBaseName()}"

    input:
    path tree

    output:
    path "${tree.getBaseName()}_topo_binary.newick", emit: topo

    publishDir "barcode_experiment3", mode: 'copy'

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
    publishDir "barcode_experiment3_pt2", mode: 'copy'

    script:
    def basename = tree.getBaseName()
    def outdir = "fastlaml_barcode_${basename}_pt2"

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

clone=3
params.matrix = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/colonies_clone${clone}_petracer_character_matrix.csv"
params.matrix_type = "character-matrix"
params.fastlaml_bin = "/Users/gc3045/git/fast-laml/build/src/fastlaml"
params.seed = 0
params.max_iter = 5000
params.tree_dir = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/colonies_inputpet/"

Channel
    .fromPath("${params.tree_dir}/*.nwk")
    .set { tree_files }

process run_fastlaml {
    tag "${tree.getBaseName()}"

    input:
    path tree

    output:
    path "fastlaml_*", emit: results
    publishDir "colonies_experiment1_cmat", mode: 'copy'

    script:
    def basename = tree.getBaseName()
    def outdir = "fastlaml_colonies_cmat_${basename}"

    """
    echo "Running fastlaml on tree: ${tree}"
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
    run_fastlaml(tree_files)
}

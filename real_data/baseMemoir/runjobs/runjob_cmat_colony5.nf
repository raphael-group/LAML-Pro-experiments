colony = 5
params.matrix = "/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/inputs/baseMemoir_colony${colony}_kde_character_matrix.csv"
params.matrix_type = "character-matrix"
params.fastlaml_bin = "/Users/gc3045/git/fast-laml/build/src/fastlaml"
params.seed = 0
params.max_iter = 5000
params.tree_dir = "/Users/gc3045/git/laml2-experiments/real_data/baseMemoir/inputs/trees/baseMemoir.colony${colony}"

Channel
    .fromPath("${params.tree_dir}/*.nwk")
    .set { tree_files }

process run_fastlaml {
    tag "${tree.getBaseName()}"

    input:
    path tree

    output:
    path "fastlaml_*", emit: results
    publishDir "outputs_baseMemoir_colony${colony}_cmat", mode: 'copy'

    script:
    def basename = tree.getBaseName()
    def outdir = "fastlaml_${basename}"

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
        -d ${params.matrix_type} \\
        -v
    echo "Finished: ${tree} -> ${outdir}"
    """
}

workflow {
    run_fastlaml(tree_files)
}

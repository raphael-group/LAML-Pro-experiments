params.matrix_dir = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/barcode_tracing/"
params.matrix_type = "observation-matrix"
params.fastlaml_bin = "/Users/gc3045/git/fast-laml/build/src/fastlaml"
params.seed = 0
params.max_iter = 5000
params.tree = "/Users/gc3045/git/laml2-experiments/real_data/PEtracer/inputs/trees/barcode/clone4_exp3.published_petracer.binary.nwk"
params.python = "/Users/gc3045/miniconda3/bin/python3" 

Channel
    .fromPath("${params.matrix_dir}/state_probs_eta_*")
    .set { matrix_files }

process run_fastlaml {
    tag "${matrix.baseName}"

    input:
    path matrix

    output:
    path "fastlaml_*", emit: results
    publishDir "barcode_experiment3", mode: 'copy'

    script:
    def tree_file = file(params.tree)
    def tree_base = tree_file.baseName
    def matrix_base = matrix.baseName
    def outdir = "fastlaml_barcode_${tree_base}_${matrix_base}"

    """
    echo "Running fastlaml on tree topology: ${tree_file} and matrix: ${matrix}"
    ${params.fastlaml_bin} \\
        --matrix ${matrix} \\
        --tree ${tree_file} \\
        --output ${outdir} \\
        --mode search \\
        --seed ${params.seed} \\
        --max-iterations ${params.max_iter} \\
        --ultrametric \\
        --min-branch-length 0.01 \\
        -d ${params.matrix_type} \\
        -v
    echo "Finished: ${tree_file} and ${matrix} -> ${outdir}"
    """
}

workflow {
    //run_fastlaml(tree_files, matrix_files)
    run_fastlaml(matrix_files)
}

process QUAST {
    tag "$meta"
    label 'process_medium'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    // found on https://github.com/Arcadia-Science/hifi2genome/blob/a969c2b78a519814c712c3d77fdd598cd5655c8f/modules/local/nf-core-modified/quast/main.nf
    // modified from nf-core module to copy the report.tsv file

    input:
    tuple val(meta), path("*")

    output:
    path "quast/*", type: 'dir'     , emit: qc
    tuple val(meta), path('*.tsv')  , emit: results
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    quast.py \\
        --output-dir quast \\
        *.fasta \\
        --threads $task.cpus \\
        $args

    ln -s quast/report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}

process SEQKIT_SEQ {
    tag "$meta.id"
    label 'process_low'


    conda "bioconda::seqkit=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0':
        'biocontainers/seqkit:2.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(sequence)
    val min_len

    output:
    path("len*")  , emit: filter
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        seq \\
        --threads $task.cpus \\
        ${sequence} \\
        --min-len $min_len >> lf_${sequence} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}

process JELLYFISH_HIST {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jellyfish:2.2.6--0' :
        'biocontainers/jellyfish:2.2.6--0' }"


    input:
    tuple val(meta), path(shortkmer)
    val kmernum

    output:
    tuple val(meta), path("*.histo")               , emit: shortkmer_hist

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    jellyfish histo -o ${kmernum}_mer_out_${prefix}.histo $shortkmer
    """
}

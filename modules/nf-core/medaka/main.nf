process MEDAKA {
    tag "$meta"
    label 'process_high_memory', 'error_ignore'

    conda "bioconda::medaka=1.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.8.0--py39h771796b_0' :
        'biocontainers/medaka:1.8.0--py39h771796b_0' }"

    input:
    tuple val(meta2), path(reads)
    tuple val(meta), path(assembly)
    val model

    output:
    path("medaka/*.fasta")              , emit: assembly
    path("medaka/versions.yml")            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    medaka_consensus \\
        -t $task.cpus \\
        -i $reads \\
        -d $assembly \\
        -m $model
        $args
    
    cd medaka
    mv consensus.fasta medaka_${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}

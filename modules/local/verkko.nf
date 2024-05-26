process VERKKO {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko:2.0--h45dadce_0' :
        'quay.io/biocontainers/verkko:2.0--h45dadce_0' }"

    input:
    tuple val(meta), path(fasta)
    path reference

    output:
    path "ragtag*/*.fasta"                           , emit: scaffolded_assembly
    tuple val(meta), path("ragtag*/*.stats")                 , emit: summary
    path  "ragtag*/versions.yml"                           , emit: versions

    script:
    def VERSION = '2.0'

    """
    ragtag.py scaffold $reference $fasta -f 1000 -o ragtag_${fasta.baseName}

    cd ragtag_${fasta.baseName}
    mv ragtag.scaffold.fasta ragtag.scaffold.${fasta.baseName}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RagTag: $VERSION
    END_VERSIONS
    """
}

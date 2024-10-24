process RECENTRIFUGE_C {
tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/recentrifuge:1.9.1--pyhdfd78af_0' :
        'biocontainers/recentrifuge:1.9.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(centrifuge_output) //path_to/centrifuge/output/files
    path rcf_db
    
    output:
    tuple val(meta), path("*.html")     , emit: html
    path "versions.yml"                , emit: versions

    script:
    def VERSION = '0.28.8'
    """
    rcf -n $rcf_db -f $centrifuge_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Recentrifuge: $VERSION
    END_VERSIONS
    """
}
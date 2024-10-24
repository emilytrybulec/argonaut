process PYCOQC {
    tag "$summary"
    label 'process_medium'

    conda "bioconda::pycoqc=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pycoqc:2.5.2--py_0' :
        'biocontainers/pycoqc:2.5.2--py_0' }"

    input:
    tuple val(meta), path(summary)
    tuple val(meta), path(inbam)
    tuple val(meta), path(inbai)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    pycoQC \\
        -f ${summary} \\
        -o ${meta.id}.html \\
        -j ${meta.id}.json \\
        -a ${inbam} \\
        --skip_coverage_plot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version 2>&1 | sed 's/^.*pycoQC v//; s/ .*\$//')
    END_VERSIONS
    """
}
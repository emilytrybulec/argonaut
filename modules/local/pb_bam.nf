process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pbtk==3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.0--h9ee0642_0':
        'quay.io/biocontainers/pbtk' }"

    input:
    tuple val(meta), path(bam, stageAs: "?/*")

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.pbi"), emit: pbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmerge \\
        -o ${prefix}.bam \\
        $args \\
        */*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbbam: \$( pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//' )
    END_VERSIONS
    """
}

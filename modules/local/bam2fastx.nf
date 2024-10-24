process BAM2FASTX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbtk==3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.0--h9ee0642_0':
        'quay.io/biocontainers/pbtk' }"

    input:
    tuple val(meta), path(bam, stageAs: '??.bam'), path(index, stageAs: '??.bam.pbi')

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam \\
        > ${prefix}.bam2fastx.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq: \$(bam2fastq --version | sed 's/bam2fastq //g')
    END_VERSIONS
    """
}
